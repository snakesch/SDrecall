#!/usr/bin/env python3

# This module contains all functions related to graphs.

import networkx as nx
import logging
logger = logging.getLogger("root")

def populate_graph_by_sd(sd_map, G):
    '''
    This function populates an empty graph G with SD relationships obtained from sd_map.
    
    Note: This step is independent of involved chromosomes.
    '''
    for _, row in sd_map.iterrows():
        G.add_edge(
            (row["chr_1"], row["start_1"], row["end_1"], row["strand1"]), # Node 1
            (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]), # Node 2
            type = "segmental duplication"
        )
    return G
    
def create_overlap_digraph(sd_map_by_chr):
    '''
    This function creates DiGraphs for SD pairs on each chromosome. Returns a dict of graphs.
    '''
    from intervaltree import Interval, IntervalTree
    
    graph_dict = {}
    # We first aggregate SD coordinates of the same contig, then build a directed graph based on overlap fraction for each contig
    for chrom, sd_map in sd_map_by_chr.items():
        directed_overlap_graph = nx.DiGraph()
        sd_interval_tree = IntervalTree()

        # For each pair of SDs, we create two nodes representing each SD component and define an interval spanning both SDs
        for _, sd_pair in sd_map.iterrows():
            node1, node2 = (sd_pair["chr_1"], sd_pair["start_1"], sd_pair["end_1"], sd_pair["strand1"]), \
                           (sd_pair["chr_2"], sd_pair["start_2"], sd_pair["end_2"], sd_pair["strand2"])

            ## Add node1 and node2 to the directed graph
            directed_overlap_graph.add_node(node1, size = sd_pair["end_1"] - sd_pair["start_1"])        
            directed_overlap_graph.add_node(node2, size = sd_pair["end_2"] - sd_pair["start_2"])        

            ## Compute the overlapped fragments by a growing tree (only consider SD on the same contig)
            interval_1, interval_2 = Interval(sd_pair["start_1"], sd_pair["end_1"], node1), \
                                     Interval(sd_pair["start_2"], sd_pair["end_2"], node2)
            intervals = filter(lambda x: x.data[0] == chrom, [interval_1, interval_2])

            ## Add intervals to the growing tree sequentially to identify overlap fragments
            for interval in intervals:
                ### Grow the tree if no overlap is found
                overlaps = sd_interval_tree.overlap(interval.begin, interval.end)
                if interval not in overlaps:
                    sd_interval_tree.add(interval)
                    for overlap in overlaps:
                        interval_size, overlap_size = interval.end - interval.begin, overlap.end - overlap.begin
                        overlap_fragment = min([interval.end, overlap.end]) - max([interval.begin, overlap.begin])
                        weight = max( overlap_fragment / overlap_size, overlap_fragment / interval_size )
                        uncovered_size = min(overlap_size, interval_size) - overlap_fragment
                        ### Here we define edge directions: from large SD to small SD
                        if interval_size >= overlap_size:
                            directed_overlap_graph.add_edge(interval.data, overlap.data, weight = weight, overlap = "True", uncovered_smaller_interval = uncovered_size)
                        else:
                            directed_overlap_graph.add_edge(overlap.data, interval.data, weight = weight, overlap = "True", uncovered_smaller_interval = uncovered_size)
                else:
                    continue
        graph_dict[chrom] = directed_overlap_graph
    
    return graph_dict
    
def create_multiplex_graph(sd_map, target_region = "", 
                           overlap_frac_min = 0.8, 
                           avg_insert_size = 350, std_insert_size = 120,
                           resolution = 0.1, 
                           thread = 10):
    '''
    This function creates two graphs, one representing SDs and one representing overlap fractions.
    A multiplex graph is returned which contains both graphs.
    '''
    
    import os
    import pandas as pd
    
    from intervaltree import Interval, IntervalTree
    import pybedtools as pb
    
    from graphs import populate_graph_by_sd
    
    ## Initialize an output directory
    os.makedirs(os.path.basename(bam_file).split(".")[0] + "_graph", exist_ok=True)
    
    ## Initialize an empty graph G
    G = nx.Graph()
    
    ## Extract unique contigs in SD map
    all_chrs = set(sd_map["chr_1"].unique().tolist()) | set(sd_map["chr_2"].unique().tolist())
    if os.path.exists(target_region):
        target_chrs = pd.read_csv(target_region, sep="\t", header=None).iloc[:, 0].unique().tolist()
        all_chrs = all_chrs & set(target_chrs)
      
    ## Populate the graph G by SD relationship
    G_sd = populate_graph_by_sd(sd_map, G)
    
    ## Initialize a dictionary of SD maps for each chromosome
    sd_map_by_chr = {}
    for chrom in all_chrs:
        sd_map_by_chr[chrom] = sd_map[(sd_map["chr_1"] == chrom) | (sd_map["chr_2"] == chrom)]
    
    ## Create a DiGraph based on overlap information
    Gdict_overlap = create_overlap_digraph(sd_map_by_chr)
    
    ## Merge all directed graphs
    G_merged = nx.DiGraph()
    for dg in Gdict_overlap.values():
        G_merged = nx.compose(G_merged, dg)
    
    ## Overlap G_overlap with SD relationship
    for node1, node2, dat in G.edges(data=True):
        if dat.get("type", None) == "segmental duplication":
            if node1 in G_merged.nodes() and node2 in G_merged.nodes():
                if G_merged.has_edge(node1, node2):
                    G_merged[node1][node2]["type"] = "segmental duplication"
                else:
                    G_merged.add_edge(node1, node2, **dat)
            else:
                logger.warning(f"These two genomic regions share sequences but are not included in the overlap graph: \n{node1}\n{node2}\n")
    
    ## Compute overlap fraction with gene regions provided
    dg_bed = pb.BedTool.from_dataframe(pd.DataFrame(G_merged.nodes()))
    if os.path.exists(target_region):
        target_bed = pb.BedTool(target_region)
        dg_bed = dg_bed.intersect(target_bed, wo=True, F=0.3, f=0.3, e=True).to_dataframe(disable_auto_names=True, names=["chr_graph", "start_graph", "end_graph", "strand_graph", "chr_feature", "start_feature", "end_feature", "overlap_len"])
    else:
        dg_bed.columns = ["chr_graph", "start_graph", "end_graph", "strand_graph"]
    logger.debug("Overlap dataframe looks like: ")
    logger.debug(dg_bed.head(5))
    
    ## Add overlap span as an attribute
    for _, region in dg_bed.iterrows():
        node = (region["chr_graph"], region["start_graph"], region["end_graph"], region["strand_graph"])
        ## node attributes are accessed by get_node_attributes method
        if node in G_merged.nodes:
            nx.set_node_attributes(G_overlap, {node: region["overlap_len"]}, name="FC_overlap")

    sd_edge_count = [e.get("type", None) for u,v,e in G_merged.edges(data=True)].count("segmental duplication")

    logger.info(f"Number of nodes representing overlap and homology: {G_merged.number_of_nodes()}")
    logger.info(f"Number of edges representing segmental duplications: {sd_edge_count}")

    nx.write_graphml(G_merged, "multiplexed_homologous_sequences.graphml")
    
    return G_merged