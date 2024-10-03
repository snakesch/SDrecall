import networkx as nx
import graph_tool.all as gt

from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys
import logging
from itertools import repeat

from homoseq_region import HOMOSEQ_REGION

logger = logging.getLogger('SDrecall')

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

def extract_FC_NFC_pairs_from_graph(query_nodes, 
                                    directed_graph, 
                                    graph_path = "", 
                                    avg_frag_size = 500, 
                                    std_frag_size = 150, 
                                    threads = 12,
                                    logger = logger):
    # First we need to annotate the graph about which nodes are query_nodes
    # query nodes are dataframe with chr, start, end ,three columns
    query_nodes = list(zip(query_nodes["chr"], query_nodes["start"], query_nodes["end"], query_nodes["strand"]))
    logger.info(f"How many query nodes we have ? {len(query_nodes)}, stored as a list of {type(query_nodes[0])}")
    logger.info(f"We got a graph with {len(directed_graph.nodes())} nodes and {len(directed_graph.edges())} edges")
    
    for node in directed_graph.nodes():
        if node in query_nodes:
            directed_graph.nodes[node]["query_node"] = True
            if len(directed_graph.nodes[node].get("FC_overlap", "")) <= 1:
                logger.warning(f"The query node {node} has not been annotated to overlap with FC region, but why it is recorded as an query_node region")
        else:
            directed_graph.nodes[node]["query_node"] = False
    logger.info("Now we have marked the query nodes in the graph")

    unfiltered_graph = directed_graph.copy()
    # Then we need to filter out nodes that are too small
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
    logger.info("Now we have removed nodes that are too small, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))
    
    # Then we filter out the edges that the overlap size is too small
    tobe_removed_edges = []
    for edge in directed_graph.edges(data=True):
        overlap = edge[2].get("overlap", False)
        if overlap:
            smaller_node = edge[1]
            smaller_node_size = smaller_node[2] - smaller_node[1]
            weight = float(edge[2]["weight"])
            overlap_size = smaller_node_size * (1/weight)
            if overlap_size < cutoff and 1/weight < 0.5:
                tobe_removed_edges.append(edge[:2])
        if edge[0] == edge[1]:
            # Remove self loop edges
            tobe_removed_edges.append(edge[:2])
    tobe_removed_edges = list(dict.fromkeys(tobe_removed_edges))
    directed_graph.remove_edges_from(tobe_removed_edges)  # This directed_graph also needs to be returned for further saving as a file
    logger.info("Now we have removed edges that implies an overlap with a size too small, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))
    
    """
    Then we can handle the BFS search to another function and perform this on each query node in parallel
    We need to sort the query_nodes to do the load balancing
    """

    assert len(query_nodes) > 0, "No query nodes are found in the graph"
    undirected_graph = directed_graph.to_undirected()
    gt_graph = convert_networkx_to_graphtool(undirected_graph)
    sorted_query_nodes = sort_query_nodes(query_nodes, undirected_graph)
    
    # Identify the subgraphs each query_node is in
    prop_map = gt_graph.vertex_properties["node_index"]
    node_to_vertex = {prop_map[v]: v for v in gt_graph.vertices()}
    components, hist = gt.label_components(gt_graph, directed = False) # Only use the returned components as a VertexPropertyMap object (mapping from vertices to properties, here the property is the label of the components)
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
        # Now we need to identify the FC-NFC pairs and annotate them on the directed graph
        i = 0
        filtered_results = set([])
        connected_qnodes_graph = nx.Graph()
        for success, tup, logs in hom_results:
            print(f"\n\n*********************************** {i}_subprocess_start_for_traverse_network ***************************************", file = sys.stderr)
            if not success:
                error_mes, tb_str = tup
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
            else:
                if len(tup) == 0:
                    logger.warning(f"No FC-NFC pairs are found for query node {tup[0]}, so we will skip this query node")
                    print(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************\n\n", file = sys.stderr)
                    i+=1
                    continue
                filtered_results.add(tup[:2])
                qnode, counterparts_nodes, query_counter_nodes = tup
                # Add edges between qnode and its query counterparts in the connected qnodes graph
                connected_qnodes_graph.add_node(qnode)
                for cnode in query_counter_nodes:
                    # Note cnode here is HOMOSEQ_REGION object
                    connected_qnodes_graph.add_edge(qnode, cnode.data)
                # Here the tuple should contain (tuple_node, [list of nodes in self-defined class HOMOSEQ_REG])
                assert isinstance(qnode, tuple)
                try:
                    assert isinstance(counterparts_nodes[0], HOMOSEQ_REGION)
                except IndexError:
                    # There is one interval (chrX, 70902050, 71018311) in the WGAC database, that it only corresponds with itself, so it will be filtered out
                    logger.warning(f"No counterparts nodes are found for query node {qnode}, so we will skip this query node")
                    print(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************\n\n", file = sys.stderr)
                    i+=1
                    continue
                
                # First we need to annotate the query node
                unfiltered_graph.nodes[qnode]["FRA_dest"] = str(i)
                # Then we need to annotate the counterparts nodes
                for cnode in counterparts_nodes:
                    unfiltered_graph.nodes[cnode.data]["FRA_source"] = (unfiltered_graph.nodes[cnode.data].get("FRA_source", "") + "," + str(i)).lstrip(",")
            print(logs, file = sys.stderr)
            print(f"*********************************** {i}_subprocess_end_for_traverse_network ***************************************\n\n", file = sys.stderr)
            i+=1
   
    # After the annotation with the FC-NFC pairs in the graph, we need to save it to a graphml file
    if graph_path:
        nx.write_graphml(unfiltered_graph, graph_path)
    # Now we need to identify the qnodes that are in the same components in the new raph connected_qnodes_graph
    qnode_components = list(nx.connected_components(connected_qnodes_graph))
    # The qnode components should be like this:
    # [ [(chr,start,end,strand), (), ..., ()],
    #   [], 
    #   [], ... ,]
    logger.info(f"The connected qnodes graph has {len(qnode_components)} subgraphs. One component has {len(qnode_components[0])} nodes look like this :\n{qnode_components[0]}\n")
    return filtered_results, qnode_components

def query_connected_nodes(fc_nfc_list, 
                          connected_qnode_components,
                          multiplex_graph = nx.Graph(), 
                          threads = 12, 
                          avg_frag_size = 500, 
                          std_frag_size = 150, 
                          graph_path = ""):
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
        transposed = list(itertools.zip_longest(*groups))
        # remove None values from each sublist (which are fill values for exhausted groups)
        result = list(dict.fromkeys([tuple([ e for e in sublist if e not in [np.nan, None]]) for sublist in transposed]))
        return result

    fc_nfc_dict = { n:ns for n, ns in fc_nfc_list }
    # We also need to build a dict that records which FC node each NFC node maps to
    fc_nodes = [ n for n, ns in fc_nfc_list]
    nfc_nodes = [ nfc for n, ns in fc_nfc_list for nfc in ns ]
    
    # The the FC nodes that can be put into the same masked genome
    component_delimited_fcnodes = pick_from_each_group(connected_qnode_components)
    logger.info("Here are the {} groups of FC nodes that can be put into the same masked genome: \n{}\n".format(len(component_delimited_fcnodes),
                                                                                                                "\n".join(str(g) for g in component_delimited_fcnodes)))
    final_fcnode_groups = component_delimited_fcnodes
    
    grouped_results = []
    for group in final_fcnode_groups:
        sd_counterparts = [nfc_node for fcnode in group for nfc_node in fc_nfc_dict[fcnode]]
        assert isinstance(sd_counterparts[0], HOMOSEQ_REGION)
        grouped_result = { "PCs":[ n for n in group ], "SD_counterparts": sd_counterparts }
        if len(group) == 0:
            logger.warning(f"There is an empty group in the component_delimited_fcnodes, please check the component_delimited_fc_nodes: \n{component_delimited_fcnodes}\n\n")
            continue
        if len(sd_counterparts) == 0:
            logger.warning(f"There is an empty set of SD counterparts in the connected_result, please check the input fc_nfc_dict: \n{fc_nfc_dict}\nAnd the fc nodes group: \n{group}\n\n")
            continue
        grouped_results.append(grouped_result)
        
    logger.info(f"There are {len(grouped_results) + 1} FC-NFC pairs for the preparation of realignment bed file\n\n")
    return grouped_results, fc_nfc_dict