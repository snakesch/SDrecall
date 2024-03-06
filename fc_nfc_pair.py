#!/usr/bin/env python3

import pandas as pd
import networkx as nx
import logging
logger = logging.getLogger("root")

class HOMOSEQ_REGION(object):
    def __init__(self, chrom, start, end, strand, ref_fasta=""):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.ref_fasta = ref_fasta
        
        # start and end coordinates relative to 0
        self.rela_start = 0
        self.rela_end = end - start
    
    @classmethod
    def from_tuple(self, tup):
        return HOMOSEQ_REGION(tup[0], tup[1], tup[2], tup[3])
    
    def to_tuple(self):
        return (self.chrom, self.start, self.end, self.strand)
    
    def to_bedtools_object(self):
        import pybedtools as pb
        import pandas as pd
        
        return pb.BedTool([self.to_tuple()])
        
    def fix_coord(self):
        return tuple([self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand])
    
    @property
    def size(self):
        return self.end - self.start
        
    def __hash__(self):
        return hash(self.to_tuple())
    
    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end}, {self.strand}"
    
    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.chrom == other.chrom and self.start == other.start and self.end == other.end and self.strand == other.strand
        elif isinstance(other, tuple):
            if len(other) != 2:
                return self.chrom == other[0] and self.start == other[1] and self.end == other[2] and self.strand == other[3]
            else:
                return self.__eq__(other[0])
        else:
            raise ValueError(f"Unrecognized comparison type ({type(other)}). ")
    
    ## This magic is critical for pybedtools
    def __iter__(self):
        for item in [self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand]:
             yield item
    
    ## Each homologous region is counted once
    def __len__(self):
        return 1
                
    def __getitem__(self):
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self.data))
            return [self.data[i] for i in range(start, stop, step)]
        else:
            return self.data[index]
        
def search_neighbor_nodes(qnode, undirected_graph, avg_frag_size, std_frag_size, visited_nodes: list[HOMOSEQ_REGION|None], extend_type: str|None=None) -> HOMOSEQ_REGION | None:
        '''This function takes a query node and returns all associated counterpart nodes from the undirected graph. '''
        
        import networkx as nx
        
        counterpart_edges = nx.bfs_edges(undirected_graph, source = qnode, depth_limit=1)
             
        counter_nodes = []
        extend_nodes = []
        
        for node1, node2 in counterpart_edges:
            cnode = HOMOSEQ_REGION.from_tuple(node2) if node1 == qnode else HOMOSEQ_REGION.from_tuple(node1)
            
            # Skip the current iteration if cnode is already visited
            if cnode in visited_nodes:
                continue
            
            # Graph attributes
            edge_attrs = undirected_graph.get_edge_data(cnode, qnode, None)
            qnode_attrs, cnode_attrs = undirected_graph.nodes[qnode], undirected_graph.nodes[cnode]
            logger.debug(f"Node attributes for query and counterpart nodes are {qnode_attrs} and {cnode_attrs} respectively. ")
            
            if edge_attrs.get("type", None) == "segmental duplication":
                ## coordinates of counterparts w.r.t. query node
                cnode.rela_start = qnode.rela_start
                cnode.rela_end = qnode.rela_end
                counter_nodes.append(cnode)
                extend_nodes.append((cnode, "SD"))
            # Case: cnode overlaps with qnode
            elif edge_attrs.get("overlap", None):
                # qnode is large
                if cnode_attrs.get("size") <= qnode_attrs.get("size"):
                    if edge_attrs.get("weight") * cnode_attrs.get("size") >= avg_frag_size - .675 * std_frag_size:
                        if cnode.start >= qnode.start and cnode.end <= qnode.end: 
                            ## qnode completely encapsulates cnode
                            cnode.rela_start = 0
                            cnode.rela_end = cnode.size
                        elif cnode.start >= qnode.start:
                            ## qnode is larger but cnode end is larger than the qnode end
                            cnode.rela_start = 0
                            cnode.rela_end = qnode.end - cnode.start
                        else:
                            ## qnode is larger but cnode start is smaller than the qnode start
                            cnode.rela_start = qnode.start - cnode.start
                            cnode.rela_end = qnode.end - cnode.start

                        if extend_type == "overlap":
                            visited_nodes.append(cnode)
                        else:
                            extend_nodes.append((cnode, "overlap"))
                    else:
                        visited_nodes.append(cnode)
                else:
                    # counterpart node is larger; qnode and cnode physically overlap with one another
                    # not necessarily that cnode completely enclosed qnode
                    # first estiamte the relative position of cnode
                    if edge_attrs.get("weight") * qnode_attrs.get("size") >= avg_frag_size - .675 * std_frag_size:
                        if cnode.start <= qnode.start and cnode.end >= qnode.end:
                            cnode.rela_start = qnode.start - cnode.start
                            cnode.rela_end = cnode.rela_start + qnode_attrs.get("size")
                        elif cnode.start <= qnode.start:
                            cnode.rela_start = qnode.start - cnode.start
                            cnode.rela_end = cnode_attrs.get("size")
                        else:
                            cnode.rela_start = 0
                            cnode.rela_end = qnode.end - cnode.start
                        if extend_type == "overlap":
                            visited_nodes.append(cnode)
                        else:
                            extend_nodes.append((cnode, "overlap"))
                    else:
                        visited_nodes.append(cnode)
         
        return counter_nodes, extend_nodes, visited_nodes
    
def get_homologous_counterparts_by_graph_traversion(qnode, undirected_graph, avg_frag_size, std_frag_size):

    import pybedtools as pb
    
    qnode = HOMOSEQ_REGION.from_tuple(qnode)

    counterparts_nodes, extend_nodes, visited_nodes = search_neighbor_nodes(qnode, undirected_graph, avg_frag_size=avg_frag_size, std_frag_size=std_frag_size, visited_nodes=[qnode], extend_type=None)

    visited_nodes = list(dict.fromkeys(visited_nodes + counterparts_nodes + extend_nodes))
    extend_nodes = {1: extend_nodes}
    n_bfs = 1

    while len([ v for vs in extend_nodes.values() for v in vs ]):
        '''
        Only case to break is if:
        1. All counterparts nodes are exhausted
        2. All extend nodes are checked
        '''

        ## BFS counter
        n_bfs += 1

        ## Before extending search depth, we sort the current node list and identify the node from which to extend
        extend_nodes = { k:sorted(vs, key = lambda t: abs(t[0].end - t[0].start - qnode.size)) for k, vs in extend_nodes.items() }
        closest_layer = min([k for k in extend_nodes.keys() if len(extend_nodes[k]) > 0])
        closest_layer_enodes = extend_nodes[closest_layer]
        enode, enode_type = closest_layer_enodes.pop(0)

        ## Search depth extension
        iter_counter_nodes, iter_extend_nodes, visited_nodes = search_neighbor_nodes(enode, undirected_graph, avg_frag_size=avg_frag_size, std_frag_size=std_frag_size, visited_nodes=visited_nodes, extend_type=enode_type)
        counterparts_nodes = counterparts_nodes + iter_counter_nodes
        extend_nodes[closest_layer + 1] = extend_nodes.get(closest_layer + 1, []) + iter_extend_nodes
        visited_nodes = list(dict.fromkeys(visited_nodes + iter_counter_nodes + iter_extend_nodes))

        assert all(isinstance(cnode, HOMOSEQ_REGION) for cnode in counterparts_nodes), "Counterparts nodes are not of type HOMOSEQ_REGION. "

        ## Compute total coverage of identified counterparts
        counter_bed_size = pb.BedTool.from_dataframe(cnode_df := pd.DataFrame(counterparts_nodes)).sort().total_coverage()
        logger.debug(f"{n_bfs}th BFS for node {enode}: acquired {len(iter_counter_nodes)} counterpart nodes and number of enodes by the end of current BFS is {len(extend_nodes.values())}. ")
        if len(counterparts_nodes) > 9 or counter_bed_size/qnode.size > 9 or closest_layer > 3:
            logger.debug(f"Stop BFS for node {qnode} because the number of counterpart nodes has reached 10. We might have already collected enough reads from other regions to the current target region. ")
            break

    # Now both qnode and its counterpart nodes ready for return
    if len(counterparts_nodes) == 0:
        logger.warning(f"This query node {qnode} cannot identify any region in the graph that shares homologous sequences. ")
    else:
        assert all(isinstance(cnode, HOMOSEQ_REGION) for cnode in counterparts_nodes), "Counterparts nodes are not of type HOMOSEQ_REGION. "
        logger.info(f"Query node {qnode} has {len(counterparts_nodes)} counterpart nodes sharing homologous sequences, including {cnode_df}")

    return (qnode, counterparts_nodes)

def extract_FC_NFC_pairs_from_graph(query_nodes, merged_directed_graph, avg_frag_size = 500, std_frag_size = 150, threads = 12):
    '''
    This function extracts homologous regions associated with a given query 
    region (qnode) using an undirected graph. Returned nodes with the same 
    FRA_* labels represent one group of homologous regions. 
    
    Annotated (Order-annotated) graph is written to multiplexed_homologous_sequences.trim.annoPC.graphml.
    '''
        
    import multiprocessing as mp
    from functools import partial
    import networkx as nx
    
    # Annotate the merged graph by the existence of query nodes
    query_nodes = query_nodes.to_records(index=False).tolist()
    for node in merged_directed_graph.nodes():
        merged_directed_graph.nodes[node]["query_node"] = True if node in query_nodes else False
    logger.debug(f"Number of query nodes: {len(query_nodes)}")
    logger.debug(f"Number of nodes: {merged_directed_graph.number_of_nodes()}; Number of edges: {merged_directed_graph.number_of_edges()}")
                     
    # Create a copy of the unannotated graph where FC-NFC pairs are not recognized
    original_graph = merged_directed_graph.copy()
    
    # Filter out nodes that are too small
    size_cutoff = max(140, avg_frag_size - .675 * std_frag_size)
    tobe_removed_nodes = []
    for node in merged_directed_graph.nodes(data=True):
        size = float(node[0][2]) - float(node[0][1])
        merged_directed_graph.nodes[node[0]]["size"] = int(size)
        if size <= size_cutoff:
            # Dont modify collections during the iteration
            tobe_removed_nodes.append(node[0])
    tobe_removed_nodes = list(dict.fromkeys(tobe_removed_nodes))
    merged_directed_graph.remove_nodes_from(tobe_removed_nodes)  
    
    # Filter out edges with small overlap size
    tobe_removed_edges = []
    for edge in merged_directed_graph.edges(data=True):
        overlap = edge[2].get("overlap", False)
        if overlap:
            smaller_node = edge[1]
            smaller_node_size = smaller_node[2] - smaller_node[1]
            weight = float(edge[2]["weight"])
            overlap_size = smaller_node_size * weight
            if overlap_size < size_cutoff and weight < 0.6:
                tobe_removed_edges.append(edge[:2])
        if edge[0] == edge[1]:
            # Remove self loop edges
            tobe_removed_edges.append(edge[:2])
    tobe_removed_edges = list(dict.fromkeys(tobe_removed_edges))
    merged_directed_graph.remove_edges_from(tobe_removed_edges)  # This directed_graph also needs to be returned for further saving as a file
    logger.debug(f"Removed {len(tobe_removed_nodes)} nodes and {len(tobe_removed_edges)} edges.")
    
    # For each query node, search for its counterpart nodes using BFS algorithm
    # The directed graph input here should be handled in advance to remove nodes that are too small or edges with low weight
    #
    # For each FC-NFC pair, FC is a single genomic interval and NFC can be a list of genomic intervals.
    # 1. Perform BFS layer by layer, each BFS only performs at depth equals 1
    # 2. Make sure the following BFS search does not traverse back to the visited noedes
    # 3. Filter out nodes based on the connecting edge type and attributes
    #
    # Notice that only nodes from the directed graph were filtered by node size, but not qnode]
    # 
    assert len(query_nodes) > 0, "No query nodes were found from the graph. "
    undirected_graph = merged_directed_graph.to_undirected()
    
    pool = mp.Pool(threads)
    partial_func = partial(get_homologous_counterparts_by_graph_traversion, undirected_graph=undirected_graph, avg_frag_size=avg_frag_size, std_frag_size=std_frag_size)
    homologous_nodes = pool.imap_unordered(partial_func, query_nodes)
    homologous_nodes = list(homologous_nodes)
    pool.close()
    pool.join()

    # After acquiring a list of homologous nodes w.r.t. a given qnode, annotate this info on the original undirected graph
    for i, (qnode, cnodes) in enumerate(homologous_nodes):
        original_graph.nodes[qnode]["FRA_dest"] = i
        
        for cnode in cnodes:
            original_graph.nodes[cnode]["FRA_source"] = i
    
    ordered_graph = original_graph
    
    nx.write_graphml(ordered_graph, "multiplexed_homologous_sequences.trim.annoPC.graphml")
    logger.info(f"Extracted {len(homologous_nodes)} groups of homologous sequences. ")
     
    return homologous_nodes

def get_qnode_combinations(qnode_groups):
    '''
    Given a qnode_group: [(qnode1, qnode2, qnode3, ...), (qnodeA, qnodeB, ...)],
    returns: [(qnode1, qnodeA), (qnode2, qnodeB), ...]
    so that qnodes in each tuple are non-overlapping
    '''
    from itertools import zip_longest
    
    transposed = list(zip_longest(*qnode_groups))
    results = list(dict.fromkeys([tuple([ e for e in sublist if e is not None]) for sublist in transposed]))
    
    return results

def bin_on_ploidy_fold_change(non_overlapping_qnodes, qnode_ploidy_dict):
        '''
        non_overlapping_qnodes: [(qnode1, qnode2, qnode3), (qnode1, qnode2), ...]
        qnode_ploidy_dict: {qnode1: 4.0, qnode2: 9.8, ...}
        
        Groups query nodes by the estimated ploidies. 
        '''
        
        updated_qnode_list = []
        for qnode_group in non_overlapping_qnodes:
            low_ploidy_group, mid_ploidy_group, high_ploidy_group, ultra_high_ploidy_group, extreme_ploidy_group = [], [], [], [], []
            for qnode in qnode_group:
                estimated_ploidy = qnode_ploidy_dict[qnode]
                if estimated_ploidy <= 3:
                    low_ploidy_group.append(qnode)
                elif estimated_ploidy <= 5:
                    mid_ploidy_group.append(qnode)
                elif estimated_ploidy <= 8:
                    high_ploidy_group.append(qnode)
                elif estimated_ploidy <= 10:
                    ultra_high_ploidy_group.append(qnode)
                else:
                    extreme_ploidy_group.append(qnode)
            logger.debug(f"Number of query nodes in low ploidy group: {len(low_ploidy_group)}")
            logger.debug(f"Number of query nodes in mid ploidy group: {len(mid_ploidy_group)}")
            logger.debug(f"Number of query nodes in high ploidy group: {len(high_ploidy_group)}")
            logger.debug(f"Number of query nodes in ultra-high ploidy group: {len(ultra_high_ploidy_group)}")
            logger.debug(f"Number of query nodes in extreme ploidy group: {len(extreme_ploidy_group)}")
            
            all_groups = [ low_ploidy_group, mid_ploidy_group, high_ploidy_group, ultra_high_ploidy_group, extreme_ploidy_group]
            group = [ g for g in all_groups if len(g) > 0]
            updated_qnode_list.extend(group)
        
        return updated_qnode_list
    
def query_connected_nodes(homologous_nodes: list[tuple[HOMOSEQ_REGION]] | None, G_merged: nx.Graph, avg_frag_size, std_frag_size):
    '''
    homologous_nodes: [(qnode, [cnode1, cnode2, cnode3]), ...]
    
    Note that a qnode can be a cnode for another qnode. 
    '''
    import pybedtools as pb
    
    # Estimate ploidy fold changes for sequence remapping
    qnode_ploidy_change = {}
    for qnode, cnodes in homologous_nodes:
        fc_bed_size = qnode.to_bedtools_object().sort().total_coverage()
        nfc_bed_size = pb.BedTool.from_dataframe(pd.DataFrame(cnodes)).sort().total_coverage()
        ploidy_fold_change = (nfc_bed_size + fc_bed_size) / fc_bed_size
        qnode_ploidy_change[qnode] = 2 * ploidy_fold_change
    
    undirected_graph = G_merged.to_undirected()
    connected_components = list(nx.connected_components(undirected_graph))
    
    # For each component, identify qnodes that are also one of the cnode
    qnode_components = []
    disconnected_result = { "PCs": set(), "SD_counterparts": set() }
    all_qnodes = list(qnode_ploidy_change.keys())
    for component in connected_components:
        coexisting_qnodes_per_component = list( set(all_qnodes) & set(component) )
        if not isinstance(coexisting_qnodes_per_component, HOMOSEQ_REGION):
            coexisting_qnodes_per_component = list( map(lambda x: HOMOSEQ_REGION.from_tuple(tuple(x)), coexisting_qnodes_per_component) )
        qnode_components.append(coexisting_qnodes_per_component)
        ## We define disconnected results (isolated query nodes) as query nodes with one and exactly one occurrence as cnode
        if len(coexisting_qnodes_per_component) == 1:
            disconnected_result["PCs"].add(coexisting_qnodes_per_component[0])
            counterparts = list(map(lambda x: HOMOSEQ_REGION.from_tuple(x), component))
            disconnected_result["SD_counterparts"] = disconnected_result["SD_counterparts"].union(counterparts)

    logger.info(f"A total of {len(disconnected_result['PCs'])} query nodes disconnected from other query nodes. ")
    
    # Identify overlapping query nodes (i.e. non-isolated query nodes)
    overlap_qnodes = list( filter(lambda node: node not in disconnected_result["PCs"], all_qnodes) )
    logger.info(f"A total of {len(overlap_qnodes)} query nodes overlaps with homologous sequences of other query nodes. ")

    # Identify query node set (more than one qnode) with multiple occurrence as cnode
    same_component_qnodes = [ qnode for qnode in qnode_components if len(qnode) > 1 ]

    # For each group of qnodes, group non-overlapping qnodes
    non_overlapping_qnode_group = get_qnode_combinations(same_component_qnodes)
    logger.info(f"Number of groups of non-overlapping qnodes: {len(non_overlapping_qnode_group)}")
    
    final_qnode_groups = bin_on_ploidy_fold_change(non_overlapping_qnode_group, qnode_ploidy_change)
    logger.info(f"Number of qnode groups by ploidy change similarity: {len(final_qnode_groups)}")
    
    # Now we write the connected nodes
    connected_results = []
    for group in final_qnode_groups:
        connected_result = {"PCs": set(), "SD_counterparts": set()}
        for g in group:
            for qnode, cnodes in homologous_nodes:
                if g == qnode:
                    connected_result["PCs"].add(qnode)
                    connected_result["SD_counterparts"] = connected_result["SD_counterparts"].union(cnodes)
                    break
    
        # Checking for connected_results
        if len(connected_result["PCs"]) == 0:
            logger.warning(f"Empty group in connected_results, please check same_component_qnodes: {same_component_qnodes} \nand non_overlapping_qnode_group: {non_overlapping_qnode_group}\n")
        if len(connected_result["SD_counterparts"]) == 0:
            logger.warning(f"Empty SD counterparts in connected_results, please check the provided homologous nodes: {homologous_nodes}\n and the related PCs/qnodes: {connected_result['PCs']}")
        connected_results.append(connected_result)
    
    return disconnected_result, connected_results

def annotate_FC_NFC_graph(graph: nx.Graph, disconnected_result, connected_results):
    
    all_results = [disconnected_result] + connected_results
    
    final_graph = graph.copy()
    for i, pair in enumerate(all_results):
        PC_tag = "PC" + str(i)
        for node in pair["PCs"]:
            ## Add attribute "FC" to the nodes
            original_tags = final_graph.nodes[node].get("FC", None)
            final_graph.nodes[node]["FC"] = ",".join(original_tags.split(",") + [PC_tag]) if original_tags else PC_tag
        
        for node in pair["SD_counterparts"]:
            ## Add attribute "NFC" to the nodes
            original_tags = final_graph.nodes[node].get("NFC", None)
            final_graph.nodes[node]["NFC"] = ",".join(original_tags.split(",") + [PC_tag]) if original_tags else PC_tag
    
    nx.write_graphml(final_graph, "final_graph.graphml")
    
    return 