import logging
import sys

import networkx as nx
import graph_tool.all as gt

from src.log import log_command
from preparation.homoseq_region import HOMOSEQ_REGION

logger = logging.getLogger('SDrecall')

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

def summarize_shortest_paths_per_subgraph(ori_qnode, 
                                          subgraph, 
                                          graph, 
                                          subgraph_label,
                                          all_qnode_vertices,
                                          avg_frag_size = 500, 
                                          std_frag_size = 150,
                                          logger = logger):
    '''
    ori_qnode: HOMOSEQ_REGION object
    subgraph is a GraphView object containing only one subgraph connected by overlap edges
    '''
    from functools import reduce
    from operator import mul

    counter_nodes = []

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
        if reduce(mul, [1 - graph.ep["weight"][e] for e in shortest_path_edges if graph.ep["type"][e] == "segmental_duplication"]) < 0.6:
            continue
        
        # If the last edge is segmental duplication, then check if the region is considerably large
        cnode = inspect_cnode_along_route(graph, 
                                          shortest_path_verts, 
                                          shortest_path_edges, 
                                          avg_frag_size, 
                                          std_frag_size, 
                                          logger = logger)
        if cnode:
            logger.debug(f"Found a new counterparts node {cnode} for query node {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes. The traversal route is {cnode.traverse_route}")
            counter_nodes.append(cnode)

    if len(counter_nodes) == 0:
        logger.debug(f"No cnodes found for query node {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes, (one of the node is {HOMOSEQ_REGION(shortest_path_verts[-1], graph)}).")
        return [], []
    else:
        logger.debug(f"Found {len(counter_nodes)} new cnodes for qnode {ori_qnode} in the subgraph {subgraph_label} containing {n} nodes, (one of the node is {HOMOSEQ_REGION(shortest_path_verts[-1], graph)})")
        return counter_nodes, [c for c in counter_nodes if c.vertex in all_qnode_vertices]

@log_command
def traverse_network_to_get_homology_counterparts(qnode, 
                                                  directed_graph, 
                                                  all_qnode_vertices,
                                                  avg_frag_size=500, 
                                                  std_frag_size=150,
                                                  logger = logger):
    '''
    The input directed graph is only the subgraph that containing the qnode (which is a tuple instead of an integer index)
    '''
    import numpy as np

    all_qnode_vertices = set(all_qnode_vertices)
    prop_map = directed_graph.vertex_properties["node_index"]
    node_to_vertex = {prop_map[v]: v for v in directed_graph.vertices()}
    
    if isinstance(qnode, tuple):
        qnode = HOMOSEQ_REGION(node_to_vertex[qnode], directed_graph)

    # First identify subgraphs only connected by overlap edges
    overlap_graph = gt.GraphView(directed_graph, efilt = lambda e: directed_graph.ep["overlap"][e] == "True")
    overlap_components, hist = gt.label_components(overlap_graph, directed = False)
    
    uniq_comp_labels = np.unique(overlap_components.a)
    subgraphs = [gt.GraphView(directed_graph, directed = False, vfilt = lambda v: overlap_components[v] == comp_label) for comp_label in uniq_comp_labels]
    
    counterparts_nodes = []
    total_query_counter_nodes = []
    for i in range(len(uniq_comp_labels)):
        cnodes, query_counter_nodes = summarize_shortest_paths_per_subgraph(qnode, 
                                                                            subgraphs[i], 
                                                                            directed_graph, 
                                                                            uniq_comp_labels[i],
                                                                            all_qnode_vertices, 
                                                                            avg_frag_size, 
                                                                            std_frag_size,
                                                                            logger = logger)
        counterparts_nodes += cnodes
        total_query_counter_nodes += query_counter_nodes
    
    if len(counterparts_nodes) == 0:
        logger.warning(f"Query node {qnode} is not associated with any cnodes.")
    else:
        if not all([type(cnode) == HOMOSEQ_REGION for cnode in counterparts_nodes]):
            logger.critical("cnodes are not of type HOMOSEQ_REGION")
            sys.exit(2)
        counter_bed_str = "\n".join([str(cnode) for cnode in counterparts_nodes])
        precise_counter_bed_str = "\n".join(["\t".join([str(field) for field in cnode]) for cnode in counterparts_nodes])
        logger.debug(f"This query node {qnode} has {len(counterparts_nodes)} counterpart nodes:\n{counter_bed_str}\nPrecise counterpart bed coordinates are:\n{precise_counter_bed_str}\n\n")

    return (qnode.data, tuple(counterparts_nodes), tuple(total_query_counter_nodes))  