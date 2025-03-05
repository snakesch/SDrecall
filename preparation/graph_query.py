import networkx as nx
import graph_tool.all as gt

from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys
import logging
from itertools import repeat
from collections import defaultdict

from preparation.graph_traversal import traverse_network_to_get_homology_counterparts
from preparation.homoseq_region import HOMOSEQ_REGION

from src.log import logger

def convert_networkx_to_graphtool(nx_graph):
    # Create a new graph-tool graph
    gt_graph = gt.Graph(directed=nx_graph.is_directed())

    # Create a mapping between node objects and their indices
    node_mapping = {node: i for i, node in enumerate(nx_graph.nodes())}

    # Add nodes to the graph-tool graph
    gt_graph.add_vertex(n=nx_graph.number_of_nodes())

    # Create a vertex property map to store the original node objects
    vprop_node_index = gt_graph.new_vertex_property("object")
    for node, index in node_mapping.items():
        vprop_node_index[gt_graph.vertex(index)] = node

    # Determine node attribute categories and data types
    node_attrs = set()
    for _, data in nx_graph.nodes(data=True):
        node_attrs.update(data.keys())

    vprop = {}
    for attr in node_attrs:
        attr_types = set(type(data[attr]) for _, data in nx_graph.nodes(data=True) if attr in data)
        if len(attr_types) == 1:
            attr_type = next(iter(attr_types))
            if attr_type == int:
                vprop[attr] = gt_graph.new_vertex_property("int")
            elif attr_type == float:
                vprop[attr] = gt_graph.new_vertex_property("double")
            elif attr_type == str:
                vprop[attr] = gt_graph.new_vertex_property("string")
            elif attr_type == bool:
                vprop[attr] = gt_graph.new_vertex_property("bool")
            else:
                vprop[attr] = gt_graph.new_vertex_property("object")
        else:
            vprop[attr] = gt_graph.new_vertex_property("object")

    # Set node attributes
    for v, data in nx_graph.nodes(data=True):
        gt_v = gt_graph.vertex(node_mapping[v])
        for attr, value in data.items():
            vprop[attr][gt_v] = value

    # Determine edge attribute categories and data types
    edge_attrs = set()
    for _, _, data in nx_graph.edges(data=True):
        edge_attrs.update(data.keys())

    eprop = {}
    for attr in edge_attrs:
        attr_types = set(type(data[attr]) for _, _, data in nx_graph.edges(data=True) if attr in data)
        if len(attr_types) == 1:
            attr_type = next(iter(attr_types))
            if attr_type == int:
                eprop[attr] = gt_graph.new_edge_property("int")
            elif attr_type == float:
                eprop[attr] = gt_graph.new_edge_property("double")
            elif attr_type == str:
                eprop[attr] = gt_graph.new_edge_property("string")
            elif attr_type == bool:
                eprop[attr] = gt_graph.new_edge_property("bool")
            else:
                eprop[attr] = gt_graph.new_edge_property("object")
        else:
            eprop[attr] = gt_graph.new_edge_property("object")
    
    # Make sure "overlap" is added as an edge attribute
    if "overlap" not in eprop:
        eprop["overlap"] = gt_graph.new_edge_property("string")

    # Set edge attributes
    for u, v, data in nx_graph.edges(data=True):
        gt_u = gt_graph.vertex(node_mapping[u])
        gt_v = gt_graph.vertex(node_mapping[v])
        e = gt_graph.add_edge(gt_u, gt_v)
        for attr, value in data.items():
            eprop[attr][e] = value

    # Assign vertex and edge property maps to the graph-tool graph
    gt_graph.vertex_properties["node_index"] = vprop_node_index
    for attr, prop in vprop.items():
        gt_graph.vertex_properties[attr] = prop
    for attr, prop in eprop.items():
        gt_graph.edge_properties[attr] = prop

    return gt_graph

def sort_query_nodes(query_nodes, graph):
    """
    Sort the query_nodes for load balancing
    """
    # Create a dictionary to store the edge count and subgraph size for each query node
    node_info = {}
    
    all_graph_nodes = list(graph.nodes)
    
    # Iterate over each query node
    for node in query_nodes:
        if node not in all_graph_nodes:
            raise ValueError(f"Node {node} not found in graph.")
        # Get the edge count for the current node
        edge_count = graph.degree(node)
        # Find the connected component (subgraph) containing the current node
        subgraph = nx.node_connected_component(graph, node)
        # Calculate the size of the subgraph
        subgraph_size = len(subgraph)
        # Store the edge count and subgraph size in the dictionary
        node_info[node] = edge_count * subgraph_size
    # Sort the query nodes based on the edge count and subgraph size
    sorted_query_nodes = sorted(query_nodes, key=lambda x: node_info[x], reverse=True)
    return sorted_query_nodes


def imap_traverse(tup_args):
    return traverse_network_to_get_homology_counterparts(*tup_args)


def extract_SD_paralog_pairs_from_graph(query_nodes, 
                                        directed_graph, 
                                        graph_path = "", 
                                        reference_fasta = "",
                                        tmp_dir = "/tmp",
                                        avg_frag_size = 500, 
                                        std_frag_size = 150, 
                                        threads = 12,
                                        logger = logger):
    '''
    Returns:
    sd_paralog_pairs: dict = {(qnode): (cnode1, cnode2, ...), ...}
    connected_qnodes_gt: graph-tool Graph object
    '''
    # query_nodes: pd.DataFrame (columns = chr, start, end, strand)
    query_nodes = list(zip(query_nodes["chr"], query_nodes["start"], query_nodes["end"], query_nodes["strand"]))
    logger.info(f"Input graph has {len(directed_graph.nodes())} nodes and {len(directed_graph.edges())} edges")
    
    for node in directed_graph.nodes():
        directed_graph.nodes[node]["query_node"] = bool(node in query_nodes)
    logger.debug("Query nodes tagged")

    unfiltered_graph = directed_graph.copy()
    # Filter out small SDs
    cutoff = max(140, avg_frag_size - 1 * std_frag_size)
    tobe_removed_nodes = []
    for node in directed_graph.nodes(data=True):
        size = float(node[0][2]) - float(node[0][1])
        directed_graph.nodes[node[0]]["size"] = int(size)
        if size <= cutoff:
            # Dont modify collections during the iteration
            tobe_removed_nodes.append(node[0])
    tobe_removed_nodes = list(dict.fromkeys(tobe_removed_nodes))
    directed_graph.remove_nodes_from(tobe_removed_nodes)
    logger.info("After removing small target SDs, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))
    
    # Remove PO edges with minimal overlap
    tobe_removed_edges = []
    for edge in directed_graph.edges(data=True):
        if edge[2].get("overlap", False):
            smaller_node = edge[1]
            smaller_node_size = smaller_node[2] - smaller_node[1]
            weight = float(edge[2]["weight"])
            overlap_size = smaller_node_size * (1/weight)
            if overlap_size < cutoff and 1/weight < 0.5:
                tobe_removed_edges.append(edge[:2])
        if edge[0] == edge[1]:
            # Remove self-to-self edges
            tobe_removed_edges.append(edge[:2])
    tobe_removed_edges = list(dict.fromkeys(tobe_removed_edges))
    directed_graph.remove_edges_from(tobe_removed_edges)  # This directed_graph also needs to be returned for further saving as a file
    logger.info("After removing edges with minimal overlap, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))
    
    """
    Then we can handle the BFS search to another function and perform this on each query node in parallel
    We need to sort the query_nodes to do the load balancing
    """

    assert len(query_nodes) > 0, "No query nodes are found in the graph"
    undirected_graph = directed_graph.to_undirected()
    gt_graph = convert_networkx_to_graphtool(undirected_graph)
    sorted_query_nodes = sort_query_nodes(query_nodes, undirected_graph)
    
    # Identify all subgraphs associated with each query_node
    prop_map = gt_graph.vertex_properties["node_index"]
    node_to_vertex = {prop_map[v]: v for v in gt_graph.vertices()}
    components, _ = gt.label_components(gt_graph, directed = False) # Only use the returned components as a VertexPropertyMap object (mapping from vertices to properties, here the property is the label of the components)
    qnode_vertices = [node_to_vertex[qnode] for qnode in sorted_query_nodes]
    qnode_components = [gt.GraphView(gt_graph, directed = False, vfilt = lambda v: components[v] == components[vertex]) for vertex in qnode_vertices ]
    gt.openmp_set_num_threads(threads)

    with Pool(threads) as pool:
        hom_results = pool.imap_unordered(imap_traverse,
                                          zip(sorted_query_nodes,
                                              qnode_components,
                                              repeat(tuple([int(v) for v in qnode_vertices])),
                                              repeat(reference_fasta),
                                              repeat(tmp_dir),
                                              repeat(avg_frag_size), 
                                              repeat(std_frag_size)))

        # Instead of creating a NetworkX graph, create a graph-tool graph directly
        connected_qnodes_gt = gt.Graph(directed=False)
        
        # Create vertex property to store node data
        v_prop = connected_qnodes_gt.new_vertex_property("object")
        connected_qnodes_gt.vertex_properties["node_data"] = v_prop
        
        # Create a mapping between qnodes and vertices in the graph-tool graph
        node_to_vertex = {}
        
        # Add all qnodes as vertices
        for qnode in sorted_query_nodes:
            v = connected_qnodes_gt.add_vertex()
            v_prop[v] = qnode
            node_to_vertex[qnode] = v
        
        i = 0
        sd_paralog_pairs = {}
        for success, tup, logs in hom_results:
            logger.debug(f"*********************************** {i}_subprocess_start_for_traverse_network ***************************************")
            if not success:
                error_mes, tb_str = tup
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
                sys.exit(1)
            else:
                if len(tup) == 0:
                    logger.warning(f"No SD-paralog pairs are found for query node {tup[0]}, so we will skip this query node")
                    print(logs, file = sys.stderr)
                    print(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************\n\n", file = sys.stderr)
                    i+=1
                    continue
                qnode, counterparts_nodes, query_counter_nodes = tup
                if len(counterparts_nodes) > 0:
                    assert isinstance(qnode, tuple)
                    assert isinstance(counterparts_nodes[0], HOMOSEQ_REGION)
                    sd_paralog_pairs[tup[0]] = tup[1]
                else:
                    logger.warning(f"No counterparts nodes are found for query node {qnode}, skipping current query node")
                    print(logs, file = sys.stderr)
                    print(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************\n\n", file = sys.stderr)
                    i+=1
                    continue
                
                # Add edges between qnode and its counterparts in the graph-tool graph
                for cnode in query_counter_nodes:
                    # cnode: HOMOSEQ_REGION object
                    cnode_data = cnode.data
                    
                    # Add vertex for cnode if it doesn't exist yet
                    if cnode_data not in node_to_vertex:
                        v = connected_qnodes_gt.add_vertex()
                        v_prop[v] = cnode_data
                        node_to_vertex[cnode_data] = v
                    
                    # Add edge between qnode and cnode
                    connected_qnodes_gt.add_edge(node_to_vertex[qnode], node_to_vertex[cnode_data])
                
            logger.debug(logs)
            logger.debug(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************")
            i+=1
        
        # Save the annotated graph to a graphml file
        if graph_path:
            nx.write_graphml(unfiltered_graph, graph_path)
        return sd_paralog_pairs, connected_qnodes_gt


def pick_from_each_group(qnode_groups):
    '''
    Identify non-homologous groups of qnodes for building masked genomes
    '''
    from itertools import zip_longest

    unique_qnodes = list(zip_longest(*qnode_groups, fillvalue=None))
    unique_qnodes = list(dict.fromkeys([tuple([ e for e in sublist if e is not None]) for sublist in unique_qnodes]))
    return unique_qnodes


def query_connected_nodes(sd_paralog_pairs, 
                          connected_qnodes_gt):
    """
    Inputs:
    sd_paralog_pairs: dict = {(qnode): (cnode1, cnode2, ...), ...}
    connected_qnodes_gt: graph-tool Graph

    Note: A query node can also be a cnode of another qnode.
    """
    # Now connected_qnodes_gt is already a graph-tool graph
    # Identify qnodes that can be put into the same masked genome
    unique_qnodes = optimal_node_grouping(connected_qnodes_gt)
    
    # Get the original node data from graph-tool vertices
    v_prop = connected_qnodes_gt.vertex_properties["node_data"]
    node_groups = []
    for group in unique_qnodes:
        # Convert vertex indices back to node data
        node_group = [v_prop[connected_qnodes_gt.vertex(v_idx)] for v_idx in group]
        node_groups.append(node_group)
    
    logger.debug("{} groups of qnodes (splitted via graph coloring) can be put into the same masked genome: \n{}".format(
        len(node_groups), "\n".join(str(g) for g in node_groups)))
    
    grouped_results = []
    for group in node_groups:
        sd_counterparts = [cnode for qnode in group for cnode in sd_paralog_pairs.get(qnode, [])]
        if len(group) == 0:
            logger.warning(f"Empty group in the unique_qnodes: {node_groups}")
        elif len(sd_counterparts) == 0:
            logger.warning(f"No SD counterparts in group {group}. Please check sd_paralog_pairs: {sd_paralog_pairs}.")
        else:
            assert isinstance(sd_counterparts[0], HOMOSEQ_REGION)
            grouped_result = { "SD_qnodes": group, "SD_counterparts": sd_counterparts }
            grouped_results.append(grouped_result)

    logger.info(f"Total {len(grouped_results) + 1} SD-paralog pairs for the preparation of realignment bed file")
    return grouped_results


def optimal_node_grouping(g, min_distance=6):
    """
    Group nodes from a graph such that:
    1. Nodes from the same component must be at least min_distance edges apart
    2. Minimize the number of groups
    
    Args:
        g: A graph-tool Graph
        min_distance: Minimum distance between nodes in same component (default: 6)
        
    Returns:
        List of node groups
    """
    # Step 1: Identify components
    component_labels, _ = gt.label_components(g)
    components = defaultdict(list)
    
    # Fix: Access component_labels with vertex object directly
    for v in g.vertices():
        comp_id = component_labels[v]  # Use vertex object directly
        components[comp_id].append(int(v))
    
    # Step 2: Build conflict graph (nodes that can't be in the same group)
    conflict_graph = gt.Graph(directed=False)
    vertex_map = {}  # Maps original vertices to conflict graph vertices
    
    # Add all vertices to conflict graph
    for comp_id, vertices in components.items():
        for v in vertices:
            vertex_map[v] = conflict_graph.add_vertex()
    
    # Step 3: Add edges for conflicts (distance < min_distance)
    for comp_id, vertices in components.items():
        if len(vertices) <= 1:
            continue
            
        # For each component, compute all-pairs shortest paths
        for i, v1 in enumerate(vertices):
            for v2 in vertices[i+1:]:
                # Find shortest path distance
                dist_map = gt.shortest_distance(g, g.vertex(v1), g.vertex(v2))
                dist = dist_map[g.vertex(v2)]
                
                # If distance is less than minimum, add conflict edge
                if 0 < dist < min_distance:
                    conflict_graph.add_edge(vertex_map[v1], vertex_map[v2])
    
    # Step 4: Color the conflict graph
    colors = gt.sequential_vertex_coloring(conflict_graph)
    
    # Step 5: Group vertices by color
    groups = defaultdict(list)
    for v in g.vertices():
        v_int = int(v)
        conflict_v = vertex_map[v_int]
        color = colors[conflict_v]
        groups[color].append(v_int)
    
    return list(groups.values())