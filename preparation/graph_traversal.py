import logging
import sys
import os

import graph_tool.all as gt
from typing import Tuple

from src.log import log_command, logger
from src.utils import executeCmd, prepare_tmp_file
from preparation.homoseq_region import HOMOSEQ_REGION


def inspect_cnode_along_route(graph, 
                              verts, 
                              edges, 
                              avg_frag_size = 500, 
                              std_frag_size = 150,
                              logger = logger):
    qnode = HOMOSEQ_REGION(verts[0], graph)
    for i in range(0, len(edges)):
        logger.debug(f"For qnode {qnode}: relative start is {qnode.rela_start} and relative end is {qnode.rela_end}")
        vertex = verts[i + 1]
        edge = edges[i]
        cnode = HOMOSEQ_REGION(vertex, graph)
        logger.debug(f"The cnode is {cnode} at the other end of the {i}th edge. The relative_start is {cnode.rela_start} and relative end is {cnode.rela_end}")
        if graph.ep["type"][edge] == "segmental_duplication":
            if qnode.strand == cnode.strand:
                cnode.ups_rela_start = max(qnode.rela_start, 0)
                cnode.ups_rela_end = min(qnode.rela_end, cnode.size)
            else:
                cnode.ups_rela_end = min(qnode.size - qnode.rela_start, cnode.size)
                cnode.ups_rela_start = max(qnode.size - qnode.rela_end, 0)
            cnode.rela_start = cnode.ups_rela_start
            cnode.rela_end = cnode.ups_rela_end
            if not cnode.rela_end > cnode.rela_start:
                logger.error(f"The relative end is smaller than the relative start for {cnode} at the end of {i}th edge, rela_start is {cnode.rela_start} and rela_end is {cnode.rela_end}. The upstream qnode is {qnode}, rela_start is {qnode.rela_start} and rela_end is {qnode.rela_end}")
                sys.exit(1)
            cnode.traverse_route = tuple(list(qnode.traverse_route) + [(qnode, "segmental_duplication")])
        elif graph.ep["overlap"][edge] == "True":
            qnode_rela_start = cnode.start - qnode.start
            qnode_rela_end = min(qnode_rela_start + cnode.size, qnode.size)
            qnode_rela_start = max(qnode_rela_start, 0)
            
            qnode.down_rela_start = qnode_rela_start
            qnode.down_rela_end = qnode_rela_end
            qnode_overlap_start = max(qnode.down_rela_start, qnode.ups_rela_start)
            qnode_overlap_end = min(qnode.down_rela_end, qnode.ups_rela_end)
            qnode_overlap_size = qnode_overlap_end - qnode_overlap_start
            qnode.rela_start = qnode_overlap_start
            qnode.rela_end = qnode_overlap_end
            
            if qnode_overlap_size <= max(avg_frag_size - std_frag_size, 140):
                # 0.4 quantile in insert size distribution
                logger.debug(f"The overlapping region for (cnode is {cnode}) query node {qnode}, ({qnode.chrom}:{qnode.start + qnode_overlap_start}-{qnode.start + qnode_overlap_end}) is too small, the traverse route ends here: {list(qnode.traverse_route) + [(qnode, 'overlap')]}")
                return None
            else:
                qnode_overlap_abs_start = qnode.start + qnode_overlap_start
                cnode_rela_start = qnode_overlap_abs_start - cnode.start
                cnode_rela_end = min(cnode.size, cnode_rela_start + qnode_overlap_size)
                cnode_rela_start = max(0, cnode_rela_start)     
                cnode.ups_rela_start = cnode_rela_start
                cnode.ups_rela_end = cnode_rela_end
                cnode.rela_start = cnode_rela_start
                cnode.rela_end = cnode_rela_end

                if cnode_rela_end - cnode_rela_start <= max(avg_frag_size - std_frag_size, 140):
                    logger.debug(f"The cnode segment ({cnode}) is too small. Skipping ... ")
                    return None
                    
                cnode.traverse_route = tuple(list(qnode.traverse_route) + [(qnode, "overlap")])

        logger.debug(f"The traverse route for {cnode} at the other end of {i}th edge is {cnode.traverse_route}")
        qnode = cnode
    return cnode


def compare_homologous_sequences(
    ori_qnode: tuple,
    cnode: HOMOSEQ_REGION,
    reference_fasta: str,
    tmp_dir: str = "/tmp",
    min_similarity: float = 0.9,
    logger: logging.Logger = logger) -> Tuple[bool, float]:
    """
    Compare genomic sequences between original query node and counterpart node.
    
    Args:
        ori_qnode: Original query node
        cnode: Counterpart node
        reference_fasta: Path to reference genome fasta
        min_similarity: Minimum sequence similarity threshold (default: 0.9)
        tmp_dir: Directory for temporary files
        logger: Logger instance
    
    Returns:
        Tuple of (is_similar: bool, similarity_score: float)
    """
    # Get relative coordinates for comparison
    qnode_rela_start, qnode_rela_end = cnode.qnode_relative_region(ori_qnode.data, logger=logger)
    qnode_region = f"{ori_qnode[0]}:{ori_qnode[1] + qnode_rela_start}-{ori_qnode[1] + qnode_rela_end}"
    qnode_strand = ori_qnode[3]
    logger.debug(f"The query node {ori_qnode} has relative coordinates {qnode_rela_start} and {qnode_rela_end} and the region is {qnode_region} at strand {qnode_strand}")

    cnode_region_tup = tuple([str(field) for field in cnode])
    cnode_region = f"{cnode_region_tup[0]}:{cnode_region_tup[1]}-{cnode_region_tup[2]}"
    cnode_strand = cnode_region_tup[3]
    logger.debug(f"The counterpart node {cnode.data} has relative coordinates {cnode.ups_rela_start} and {cnode.ups_rela_end} and the region is {cnode_region} at strand {cnode_strand}")

    if qnode_rela_start == "NaN" or qnode_rela_end == "NaN":
        logger.warning(f"Invalid relative coordinates for {cnode.data}")
        return False, 0.0

    # Create temporary FASTA files for the regions
    tmp_q = prepare_tmp_file(suffix=".fasta", tmp_dir=tmp_dir)
    tmp_c = prepare_tmp_file(suffix=".fasta", tmp_dir=tmp_dir)
    
    try:
        # Extract query sequence
        cmd_q = f"samtools faidx {reference_fasta} "
        cmd_q += f"{qnode_region}"
        cmd_q += f" > {tmp_q.name}"
        executeCmd(cmd_q, logger=logger)

        # Extract counterpart sequence
        cmd_c = f"samtools faidx {reference_fasta} "
        cmd_c += f"{cnode_region}"
        # Handle different strands
        if cnode_strand != qnode_strand:
            cmd_c += " -i"
        cmd_c += f" > {tmp_c.name}"
        executeCmd(cmd_c, logger=logger)

        # Close files to ensure writing is complete
        tmp_q.close()
        tmp_c.close()

        # Run minimap2 for alignment
        tmp_paf = prepare_tmp_file(suffix=".paf", tmp_dir=tmp_dir)
        cmd_map = f"minimap2 -x asm10 -t 1 --eqx --cs -c {tmp_q.name} {tmp_c.name} > {tmp_paf.name} && rm {tmp_q.name} {tmp_c.name}"
        executeCmd(cmd_map, logger=logger)

        # Parse PAF file to calculate similarity
        with open(tmp_paf.name) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                    
                # Extract alignment statistics
                query_len = int(fields[1])
                matches = int(fields[9])
                aln_len = int(fields[10])
                
                # Calculate similarity score
                similarity = matches / aln_len if aln_len > 0 else 0
                
                # Check if similarity meets threshold
                is_similar = similarity > min_similarity
                
                logger.debug(f"Sequence similarity between qnode {qnode_region} and cnode {cnode_region}: {similarity:.3f}")
                return is_similar, similarity

        # If no alignment found
        logger.warning(f"No alignment found between qnode {qnode_region} and cnode {cnode_region}")
        return False, 0.0

    finally:
        # Cleanup temporary files
        for tmp_file in [tmp_q, tmp_c, tmp_paf]:
            try:
                os.unlink(tmp_file.name)
            except:
                pass
    return False, 0.0



def get_overlap_neighbors(g, vertex):
    """
    Directly inspect edges for the overlap property
    
    Args:
        g: graph-tool Graph object
        vertex: vertex object or index
        
    Returns:
        List of vertex objects connected via "True" overlap edges
    """
    # Convert vertex index to vertex object if needed
    if isinstance(vertex, int):
        vertex = g.vertex(vertex)
    
    data_prop = g.vertex_properties["node_index"]
    vertex_data = data_prop[vertex]  # Tuple with 4 items (chrom, start, end, strand)

    overlap_neighbors = []
    
    # Check out-edges
    for neighbor in g.get_all_neighbors(vertex):
        neighbor = g.vertex(neighbor)
        neighbor_data = data_prop[neighbor]
        overlap_neighbors.append(neighbor_data)

    return overlap_neighbors



def summarize_shortest_paths_per_subgraph(ori_qnode, 
                                          subgraph, 
                                          graph, 
                                          subgraph_label,
                                          all_qnode_vertices,
                                          reference_fasta,
                                          tmp_dir = "/tmp",
                                          avg_frag_size = 500, 
                                          std_frag_size = 150,
                                          logger = logger):
    '''
    ori_qnode: HOMOSEQ_REGION object
    subgraph is a GraphView object containing only one subgraph connected by overlap edges
    '''
    from functools import reduce
    from operator import mul

    # Create a boolean vertex property to mark query nodes
    is_query_node = subgraph.new_vertex_property("bool")
    
    # Set the property to True for vertices that are query nodes
    for v in subgraph.vertices():
        v_int = int(v)
        is_query_node[v] = v_int in all_qnode_vertices
    
    # Create a filtered view of the subgraph containing only query nodes
    query_nodes_view = gt.GraphView(subgraph, vfilt=is_query_node, efilt = lambda e: subgraph.ep["overlap"][e] == "True")
    
    # Store the property in the subgraph for future use
    subgraph.vertex_properties["is_query_node"] = is_query_node

    counter_nodes = []
    cnode_similarities = []
    counter_qnodes = []

    n = 0
    for v in subgraph.vertices():
        n += 1
        # Find the shortest path in original graph
        if int(v) == ori_qnode.vertex:
            shortest_path_verts = [graph.vertex(ori_qnode.vertex)]
            continue
        
        # all_shortest_path() function is too time consuming because it uses bfs search by default to find all shortest paths between a pair of vertices.
        shortest_path_verts, shortest_path_edges = gt.shortest_path(graph, 
                                                                    source=graph.vertex(ori_qnode.vertex), 
                                                                    target=v,
                                                                    weights=graph.ep["weight"])

        # Skip the current path if the last edge is a PO edge
        if graph.ep["overlap"][shortest_path_edges[-1]] == "True":
            continue
        # Skip the current path in case of adjacent PO edges
        if len([i for i in range(0, len(shortest_path_edges) - 1) if graph.ep["overlap"][shortest_path_edges[i]] == "True" and graph.ep["overlap"][shortest_path_edges[i+1]] == "True"]) > 1:
            continue
        if reduce(mul, [1 - graph.ep["weight"][e] for e in shortest_path_edges if graph.ep["type"][e] == "segmental_duplication"]) <= 0.9:
            continue
        
        # If the last edge is segmental duplication, then check if the region is considerably large
        cnode = inspect_cnode_along_route(graph, 
                                          shortest_path_verts, 
                                          shortest_path_edges, 
                                          avg_frag_size, 
                                          std_frag_size, 
                                          logger = logger)
        
        if cnode:
            is_similar, similarity = compare_homologous_sequences(ori_qnode, 
                                                                  cnode, 
                                                                  reference_fasta, 
                                                                  tmp_dir = tmp_dir,
                                                                  min_similarity = 0.9, 
                                                                  logger = logger )
        else:
            is_similar = False
            similarity = 0.0
        
        if is_similar:
            logger.debug(f"Found a new counterparts node {cnode} for query node {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes. The similarity is {similarity}. The traverse route is {cnode.traverse_route} \n")
            counter_nodes.append(cnode)
            cnode_similarities.append(similarity)
            if cnode.vertex in all_qnode_vertices:
                # Grab the overlapping query nodes with the current cnode and add them to the counter_qnodes list
                counter_qnodes.append(cnode.data)
                counter_qnodes += get_overlap_neighbors(query_nodes_view, cnode.vertex)
                
    if len(counter_nodes) == 0:
        logger.debug(f"No cnodes found for query node {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes, (one of the node is {HOMOSEQ_REGION(shortest_path_verts[-1], graph)}).")
        return [], []
    else:
        logger.debug(f"Found {len(counter_nodes)} new cnodes for qnode {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes, (one of the node is {HOMOSEQ_REGION(shortest_path_verts[-1], graph)})")
        return zip(counter_nodes, cnode_similarities), counter_qnodes


@log_command
def traverse_network_to_get_homology_counterparts(qnode, 
                                                  component_graph, 
                                                  all_qnode_vertices,
                                                  reference_fasta,
                                                  tmp_dir = "/tmp",
                                                  avg_frag_size=500, 
                                                  std_frag_size=150,
                                                  logger = logger):
    '''
    The input directed graph is only the subgraph that containing the qnode (which is a tuple instead of an integer index)
    '''
    import numpy as np

    all_qnode_vertices = set(all_qnode_vertices)
    prop_map = component_graph.vertex_properties["node_index"]
    node_to_vertex = {prop_map[v]: v for v in component_graph.vertices()}
    
    if isinstance(qnode, tuple):
        qnode = HOMOSEQ_REGION(node_to_vertex[qnode], component_graph)

    # First identify subgraphs only connected by overlap edges
    overlap_graph = gt.GraphView(component_graph, efilt = lambda e: component_graph.ep["overlap"][e] == "True")
    overlap_components, hist = gt.label_components(overlap_graph, directed = False)
    
    uniq_comp_labels = np.unique(overlap_components.a)
    subgraphs = [gt.GraphView(component_graph, directed = False, vfilt = lambda v: overlap_components[v] == comp_label) for comp_label in uniq_comp_labels]
    
    counterparts_nodes = []
    total_query_counter_nodes = []
    logger.debug(f"Traversing the network to find homology counterparts for {qnode} in {len(uniq_comp_labels)} components")
    for i in range(len(uniq_comp_labels)):
        cnodes, query_counter_nodes = summarize_shortest_paths_per_subgraph(qnode, 
                                                                            subgraphs[i], 
                                                                            component_graph, 
                                                                            uniq_comp_labels[i],
                                                                            all_qnode_vertices, 
                                                                            reference_fasta,
                                                                            tmp_dir = tmp_dir,
                                                                            avg_frag_size = avg_frag_size, 
                                                                            std_frag_size = std_frag_size,
                                                                            logger = logger)
        counterparts_nodes += [t[0] for t in cnodes if t[1] >= 0.95]
        total_query_counter_nodes += query_counter_nodes  # Use loosen similarity threshold to include more query nodes, improving graph coloring for better partitioning the realign groups
    
    if len(counterparts_nodes) == 0:
        logger.warning(f"Query node {qnode} is not associated with any cnodes.")
    else:
        if not all([type(cnode) == HOMOSEQ_REGION for cnode in counterparts_nodes]):
            logger.critical("cnodes are not of type HOMOSEQ_REGION")
            sys.exit(2)
        counter_bed_str = "\n".join([str(cnode) for cnode in counterparts_nodes])
        precise_counter_bed_str = "\n".join(["\t".join([str(field) for field in cnode]) for cnode in counterparts_nodes])
        logger.debug(f"This query node {qnode} has {len(counterparts_nodes)} counterpart nodes and {len(total_query_counter_nodes)} of them are query nodes:\n{counter_bed_str}\nPrecise counterpart bed coordinates are:\n{precise_counter_bed_str}\n\n")

    return (qnode.data, tuple(counterparts_nodes), tuple(total_query_counter_nodes))  