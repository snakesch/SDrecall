import logging
import gc
import numpy as np
import graph_tool.all as gt
from numba import get_num_threads
from collections import defaultdict


from bk_algorithm import bk_algorithm
from numba_operators import numba_isin, \
							numba_and, \
							numba_sum, \
							apply_index_mask



"""
Haplotype Phasing Module for Read Clustering

This module implements the steps of haplotype phasing by grouping reads
into distinct haplotypes using graph-based methods. It works in conjunction with
the Bron-Kerbosch algorithm implementation (bk_algorithm.py) and follows the graph building
module (graph_build.py).

Key features:
- Utilizes graph-tool library for efficient graph operations
- Implements a clique-based approach for identifying haplotype clusters
- Perform twice BK algorithm to balance accuracy and sensitivity

Main functions:
- clique_generator_per_component: Identifies cliques within a given component
- find_components_inside_filtered_cliques: Removing zero weight edges and identifying components inside the cliques

Dependencies:
- graph_tool
- numpy
- bk_algorithm (custom module for Bron-Kerbosch algorithm)
- numba_operators (custom module for optimized operations)

This module contains the high-level framework of haplotype phasing pipeline, taking the
phasing graph and weight matrix as input and producing the final haplotype
assignments for each read.
"""



logger = logging.getLogger("SDrecall")



def graph_vertex_iter(vertex_indices, graph):
    for vid in vertex_indices:
        yield graph.vertex(vid)



def clique_generator_per_component(graph, weight_matrix, ew_cutoff = 0.101, logger = logger):
    '''
    This function is to find the largest cliques in each component of the graph

    - One critical consideration is that when extending cliques, if we allow two read pairs with low edge weight to be grouped into the same clique, there is a significant chance that this call is a false positive.
    - Low edge weight means the overlapping region between two read pairs does not contain sufficient variants to confidently determine the phase.
    - The clique identification is performed by the BK algorithm. Which starts with a node and extend cliques as much as possible.
    - When applying an edge weight cutoff during the clique extension process, we make sure the cliques identified are all high confidence.
    - This ensures that phasing accuracy is high. But for read pairs containing less or no variants, they are no connected by edges with sufficiently high edge weight, leading to exclusion from any clique extensions.
    - This also means for haplotypes closer to the reference genome, their reads are less likely to be grouped into the clique, thus stay fragmented.

    - To address this issue, we perform twice BK algorithm with two levels of edge weight cutoffs.
    - In the first round, we apply an edge weight cutoff to ensure two read pairs share at least two SNVs or 1 indels.
    - In the second round, we apply an edge weight cutoff to ensure two read pairs share at least one SNV.
    - After the first round, during the second round, the remaining fragmented read pairs are then grouped into cliques, to increase the phasing sensitivity a little bit for read pairs carrying less variants.

    As to the read pairs carrying no variants, they are just gonna stay fragmented. This limitation is rooted from the NGS technology.

    Since misaligned haplotypes tend to have more variants than correctly aligned haplotypes, read pairs from misaligned haplotypes are more likely to be grouped into the same clique accurately because overlapping regions among them tend to have more variants, thus higher edge weight.
    '''


    # Find all components for this graph
    components, _ = gt.label_components(graph)
    component_dict = defaultdict(set)
    for v in graph.vertices():
        cid = components[v]
        component_dict[cid].add(np.int32(v))

    # Find row indices where the entire row is below or equal to 0.1
    small_row_mask = np.all(weight_matrix <= 0.1, axis=1)
    big_row_mask = np.logical_not(small_row_mask)
    small_row_indices = set(np.where(small_row_mask)[0])

    # Ensure the weight matrix is contiguous in memory, which consumes less memory and is faster to process
    if not weight_matrix.flags['C_CONTIGUOUS']:
        weight_matrix = np.ascontiguousarray(weight_matrix)

    logger.info(f"Found {len(component_dict)} components in the graph. The big weight matrix has {numba_sum(big_row_mask)} rows and columns. The original weight_matrix has {weight_matrix.shape[0]} rows and columns and its contiguity in memory is {weight_matrix.flags['C_CONTIGUOUS']}. The small row indices are {small_row_indices}")

    for component_id, component_verts in component_dict.items():
        # Iterate through each component
        comp_index_mask = numba_isin(np.arange(weight_matrix.shape[0], dtype=np.int32), component_verts)
        comp_index_mask = numba_and(comp_index_mask, big_row_mask)
        selected_indices = apply_index_mask(weight_matrix.shape[0], comp_index_mask)
        # Extract the subgraph of the component (extraction performed on the original weight matrix)
        big_weight_matrix = weight_matrix[np.ix_(comp_index_mask, comp_index_mask)]

        logger.info(f"Start to find the largest clique in the component {component_id}, which contains {len(component_verts)} vertices. Contiguous? {big_weight_matrix.flags['C_CONTIGUOUS']}. Does the current matrix a view of the original one ? {big_weight_matrix.base is weight_matrix}")
        # logger.debug(f"The selected indices are {selected_indices.tolist()}, here are the component vertices: {component_verts}")
        if len(component_verts) <= 3:
            # If the component is too small, then we add all the vertices in the component to a collection of indices which will be used for clique identification in the next round
            small_row_indices.update(component_verts)
            logger.info(f"Adding the {len(component_verts)} vertices in the component {component_id} to the small row indices")
            continue
        else:
            cliques_iter = bk_algorithm(selected_indices, big_weight_matrix, cutoff = ew_cutoff, logger = logger)
            # logger.debug(f"Found {len(cliques)} cliques in the component {component_id}\n")
            for clique in cliques_iter:
                logger.info(f"Receiving a clique containing {[graph.vertex_properties['qname'][qid] for qid in clique]} in the component {component_id}")
                if len(clique) <= 3:
                    logger.info(f"Found {[graph.vertex_properties['qname'][qid] for qid in clique]} read pairs in a very small clique.")
                    small_row_indices.update(clique)
                    logger.info(f"Adding the {len(clique)} vertices in the clique to the small row indices")
                    continue
                else:
                    logger.info(f"The clique is big enough to be directly yielded")
                    yield clique

    logger.info(f"Remaining {len(small_row_indices)} vertices that are not included in the cliques. Here we find cliques again among them:\n{small_row_indices}")

    # Now this is the second round of clique finding, assembling the remaining read pairs to increase the phasing sensitivity a little bit
    if len(small_row_indices) > 0:
        small_row_mask = numba_isin(np.arange(weight_matrix.shape[0], dtype=np.int32), small_row_indices)
        selected_indices = apply_index_mask(weight_matrix.shape[0], small_row_mask)
        small_weight_matrix = weight_matrix[np.ix_(small_row_mask, small_row_mask)]

        logger.info(f"Start to find the largest clique in the small weight matrix, which contains {numba_sum(small_row_mask)} rows and columns. Contiguous? {small_weight_matrix.flags['C_CONTIGUOUS']}")
        cliques_iter = bk_algorithm(selected_indices, small_weight_matrix, cutoff = ew_cutoff/2, logger = logger)
        for clique in cliques_iter:
            yield clique



def find_components_inside_filtered_cliques(final_cliques,
                                            graph,
                                            final_components,
                                            weight_matrix,
                                            ew_cutoff = 0.101,
                                            logger=logger):
    '''
    Drop the zero weight edges for the graph before performing this function
    For each clique, use vfilter to extract the subgraph, then use efilter to only view the weight > 0 edges
    Then use gt.label_components to find the connected components

    These components are the final groups of read pairs belonged to the same haplotype
    '''

    haplotype_idx = 0
    for clique in final_cliques:
        # Prepare a boolean v property map used for vfilt
        vertex_filter = graph.new_vertex_property("bool", val=False)
        # Set the filter to True for vertices in the list
        for vid in clique:
            vertex_filter[graph.vertex(vid)] = True

        subgraph = gt.GraphView(graph, vfilt = vertex_filter)
        clique = list(clique)
        # Supplement the small weight edges to the graph
        for i in range(len(clique)):
            for j in range(i+1, len(clique)):
                    if 0.05 < weight_matrix[clique[i], clique[j]] <= ew_cutoff:
                        logger.info(f"Adding an edge between {graph.vertex_properties['qname'][graph.vertex(clique[i])]}, {clique[i]} and {graph.vertex_properties['qname'][graph.vertex(clique[j])]}, {clique[j]} because their weight in matrix at {clique[i]} row and {clique[j]} column is {weight_matrix[clique[i], clique[j]]}")
                        v1 = subgraph.vertex(clique[i])
                        v2 = subgraph.vertex(clique[j])
                        subgraph.add_edge(v1, v2)

        component_labels, hist = gt.label_components(subgraph)
        components_dict = {}
        components_dict = defaultdict(set)
        for v in graph_vertex_iter(clique, subgraph):
            c_index = component_labels[v]
            components_dict[c_index].add(int(v))

        logger.info(f"Found {len(components_dict)} components in the clique {clique}")
        for c_index, vs in components_dict.items():
            for v in vs:
                final_components[v] = haplotype_idx
            qnames = "\n".join([graph.vertex_properties['qname'][graph.vertex(v)] for v in vs])
            logger.info(f"Assigning \n{qnames}\nto the component group {haplotype_idx}")
            haplotype_idx += 1

    # logger.info(f"This is the final components: {final_components}")
    return final_components



def phasing_realigned_reads(phased_graph, weight_matrix, edge_weight_cutoff, logger = logger):
    logger.info(f"Now start finding haplotypes in the setup weight matrix, the numba parallel threads are set to {get_num_threads()}")
    total_cliques = clique_generator_per_component( phased_graph,
                                                    weight_matrix,
                                                    ew_cutoff = edge_weight_cutoff,
                                                    logger = logger )

    total_cliques = list(total_cliques)

    # clique_sep_component_idx = 0
    qname_hap_info = defaultdict(int)
    # qname_hap_info is a map (qname -> haplotype index)
    qname_hap_info = find_components_inside_filtered_cliques( total_cliques,
                                                              phased_graph,
                                                              qname_hap_info,
                                                              weight_matrix,
                                                              edge_weight_cutoff,
                                                              logger = logger )
    gc.collect()

    logger.info(f"The final components are {qname_hap_info}")
    hap_qname_info = defaultdict(set)
    # hap_qname_info is a map (haplotype index -> set of qnames)
    for vid, hid in qname_hap_info.items():
        qname = phased_graph.vp.qname[vid]
        hap_qname_info[hid].add(qname)
    logger.info(f"The final haplotype clusters are {hap_qname_info}")

    # Return the maps for both directions
    return qname_hap_info, hap_qname_info