import networkx as nx
import graph_tool.all as gt

from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys
import logging
from itertools import repeat

from graph_traversal import traverse_network_to_get_homology_counterparts
from homoseq_region import HOMOSEQ_REGION

logger = logging.getLogger('SDrecall')

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
    # Iterate over each query node
    for node in query_nodes:
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
                                    avg_frag_size = 500, 
                                    std_frag_size = 150, 
                                    threads = 12,
                                    logger = logger):
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
                                              repeat(avg_frag_size), 
                                              repeat(std_frag_size)))

        # Identify the paralogous regions (counterparts) associated with each SD (qnode) and annotate them on the directed graph
        i = 0
        filtered_results = set([])
        connected_qnodes_graph = nx.Graph()
        for success, tup, logs in hom_results:
            logger.debug(f"*********************************** {i}_subprocess_start_for_traverse_network ***************************************")
            if not success:
                error_mes, tb_str = tup
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
                sys.exit(1)
            else:
                if len(tup) == 0:
                    logger.debug(f"No SD-paralog pairs are found for query node {tup[0]}, so we will skip this query node")
                    logger.debug(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************")
                    i+=1
                    continue
                filtered_results.add(tup[:2])
                qnode, counterparts_nodes, query_counter_nodes = tup
                # Add edges between qnode and its counterparts in the connected qnodes graph
                connected_qnodes_graph.add_node(qnode)
                for cnode in query_counter_nodes:
                    # cnode: HOMOSEQ_REGION object
                    connected_qnodes_graph.add_edge(qnode, cnode.data)
                # qnode: tuple = (tuple_node, [list of nodes in self-defined class HOMOSEQ_REG])
                assert isinstance(qnode, tuple)
                try:
                    assert isinstance(counterparts_nodes[0], HOMOSEQ_REGION)
                except IndexError:
                    # WGAC contains one interval (chrX, 70902050, 71018311) which only corresponds with itself, and is thus filtered out
                    logger.debug(f"No counterparts nodes are found for query node {qnode}, skipping current query node")
                    logger.debug(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************")
                    i+=1
                    continue
                
                # First annotate the query node
                unfiltered_graph.nodes[qnode]["FRA_dest"] = str(i)
                # Then annotate the counterparts nodes
                for cnode in counterparts_nodes:
                    unfiltered_graph.nodes[cnode.data]["FRA_source"] = (unfiltered_graph.nodes[cnode.data].get("FRA_source", "") + "," + str(i)).lstrip(",")
            logger.debug(logs)
            logger.debug(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************")
            i+=1
   
    # Save the annotated graph to a graphml file
    if graph_path:
        nx.write_graphml(unfiltered_graph, graph_path)
    # Identify the qnodes that are in the same graph component in the new connected_qnodes_graph
    qnode_components = list(nx.connected_components(connected_qnodes_graph))
    # The qnode components should be like this:
    # [ [(chr,start,end,strand), (), ..., ()],
    #   [], 
    #   [], ... ,]
    logger.debug(f"The connected qnode graph has {len(qnode_components)} subgraphs. One component has {len(qnode_components[0])} nodes:\n{qnode_components[0]}")
    return filtered_results, qnode_components

def query_connected_nodes(sd_paralog_list, 
                          connected_qnode_components):
    """
    Input fc_nfc_list: list of tuples (fc_node, nfc_nodes)
    Each tuple contains:
        fc_node: (chr, start, end, strand)
        nfc_nodes: list of HOMOSEQ_REGION objects
    Note that an FC node can also be an NFC node.
    """

    def pick_from_each_group(groups):
        '''
        Input groups should be a list of list
        '''
        import numpy as np
        from itertools import zip_longest

        # Each group in the groups is a list of tuples. Each tuple represents a query node
        transposed = list(zip_longest(*groups))
        # remove dummy None values from each sublist
        result = list(dict.fromkeys([tuple([ e for e in sublist if e not in [np.nan, None]]) for sublist in transposed]))
        return result
    
    print(sd_paralog_list)
    sd_paralog_pairs = { n:ns for n, ns in sd_paralog_list }
    
    # The the FC nodes that can be put into the same masked genome
    component_delimited_fcnodes = pick_from_each_group(connected_qnode_components)
    logger.info("Here are the {} groups of FC nodes that can be put into the same masked genome: \n{}\n".format(len(component_delimited_fcnodes),
                                                                                                                "\n".join(str(g) for g in component_delimited_fcnodes)))
    final_fcnode_groups = component_delimited_fcnodes
    
    grouped_results = []
    for group in final_fcnode_groups:
        sd_counterparts = [nfc_node for fcnode in group for nfc_node in sd_paralog_pairs[fcnode]]
        assert isinstance(sd_counterparts[0], HOMOSEQ_REGION)
        grouped_result = { "PCs":[ n for n in group ], "SD_counterparts": sd_counterparts }
        if len(group) == 0:
            logger.warning(f"There is an empty group in the component_delimited_fcnodes, please check the component_delimited_fc_nodes: \n{component_delimited_fcnodes}\n\n")
            continue
        if len(sd_counterparts) == 0:
            logger.warning(f"There is an empty set of SD counterparts in the connected_result, please check the input sd_paralog_pairs: \n{sd_paralog_pairs}\nAnd the fc nodes group: \n{group}\n\n")
            continue
        grouped_results.append(grouped_result)
        
    logger.info(f"There are {len(grouped_results) + 1} SD-paralog pairs for the preparation of realignment bed file\n\n")
    return grouped_results, sd_paralog_pairs