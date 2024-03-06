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
        dg_bed = dg_bed.to_dataframe(disable_auto_names=True, names = ["chr_graph", "start_graph", "end_graph", "strand_graph"])
    logger.debug("Overlap dataframe looks like: ")
    logger.debug(dg_bed.head(5))
    
    ## Add overlap span as an attribute
    for _, region in dg_bed.iterrows():
        node = (region["chr_graph"], region["start_graph"], region["end_graph"], region["strand_graph"])
        ## node attributes are accessed by get_node_attributes method
        if node in G_merged.nodes and os.path.exists(target_region):
            nx.set_node_attributes(G_merged, {node: region["overlap_len"]}, name="FC_overlap")

    sd_edge_count = [e.get("type", None) for u,v,e in G_merged.edges(data=True)].count("segmental duplication")

    logger.info(f"Number of nodes representing overlap and homology: {G_merged.number_of_nodes()}")
    logger.info(f"Number of edges representing segmental duplications: {sd_edge_count}")

    nx.write_graphml(G_merged, "multiplexed_homologous_sequences.graphml")
    
    return G_merged

def overlap_graph_nodes_with_multialign_bed(G_merged, avg_frag_size = 400, std_frag_size = 100, target_regions="selected_regions.bed"):
    '''
    Takes a merged graph (overlap + SD) and intersect with regions selected from alignment files in the previous steps. 
    A BED-like dataframe is returned, containing all regions of interest.
    A dataframe of query nodes is also returned containing ambiguously mapped regions relevant to the sample.
    '''
    from itertools import combinations
    import pandas as pd
    
    import pybedtools as pb
    
    ## Extract overlapping / homologous regions deduced from merged graph
    graph_bed = pb.BedTool.from_dataframe(pd.DataFrame(G_merged.nodes)).sort()

    ## Filter out overlapping / homologous regions with interval length >= 140 (common read length is 150bp)
    graph_bed = graph_bed.filter(lambda span: len(span) >= 140)

    ## Overlap with target regions defined in previous steps
    intersect_df = graph_bed.intersect(pb.BedTool(target_regions).sort().merge(), wo=True) \
                            .to_dataframe(disable_auto_names=True, names=["chr", "start", "end", "strand", "chr_query", "start_query", "end_query", "overlap_len"])
    intersect_df = intersect_df.loc[(intersect_df["end"] > intersect_df["start"]) & (intersect_df["end_query"] > intersect_df["start_query"]), :]

    ## For each query regions (defined by sample alignments and provided target regions),
    ## we try to fetch the graph regions with substantial overlapping (>0.95) 
    selected_rows = pd.DataFrame()
    for qr, gr in intersect_df.groupby(["chr_query", "start_query", "end_query"]):
        ### Take query intervals as *fixed intervals*
        fixed_chr, fixed_start, fixed_end = qr[0], int(qr[1]), int(qr[2])
        ### Fetch graph regions and sort the list by end coordinates
        gr_intervals = gr.iloc[:, :4].sort_values(by=["end"])
        ### Compute the span of overlapped fragments
        gr_with_good_coverages = []
        for _, r in gr_intervals.iterrows():
            overlap_span = min(r[2], fixed_end) - max(r[1], fixed_start)
            r_size, fixed_size = r[2] - r[1], fixed_end - fixed_start
            coverage = max(overlap_span/fixed_size, overlap_span/r_size) if overlap_span >= 0 else 0
            #### Only keep overlapped fragments with coverage > 0.95
            if coverage > 0.95:
                gr_with_good_coverages.append(r)
        ### Pick the shortest interval
        if len(gr_with_good_coverages) > 0:
            selected_interval = sorted(gr_with_good_coverages, key = lambda r: abs(r[2]-r[1]-(avg_frag_size + 1.96*std_frag_size)))[0]
        else:
            ### Deals with regions with minimal graph-query overlap
            logger.warning(f"No intervals were found to significantly overlap with query/fixed interval {fixed_chr}:{fixed_start}-{fixed_end} , attempting to find the best graph interval combination.")
            gr_intervals["peripheral_interval_size"] = gr_intervals["end"] - gr_intervals["start"]
            
            ### Compute the overlap size w.r.t. each graph region
            qr_bed_obj = pb.BedTool("\t".join(map(str, qr)), from_string=True)
            gr_bed_obj = pb.BedTool.from_dataframe(gr_intervals)
            gr_intervals_with_overlap = gr_bed_obj.intersect(qr_bed_obj, wo=True).to_dataframe()
            gr_intervals_with_overlap = gr_intervals_with_overlap.iloc[:, [0,1,2,3,8,4]]
            
            gr_intervals_with_overlap.columns = ["chrom", "start", "end", "strand", "central_interval_overlap", "peripheral_interval_size"]
            gr_intervals_with_overlap["central_interval_overlap"] = gr_intervals_with_overlap["central_interval_overlap"] / (fixed_end - fixed_start)
            gr_intervals_with_overlap = gr_intervals_with_overlap.sort_values(by=["central_interval_overlap","peripheral_interval_size"], ascending=[False, True]).reset_index(drop=True)
            
            ### Select the graph region with most overlap against the fixed regions
            if gr_intervals_with_overlap.loc[0, "central_interval_overlap"] >= 0.95:
                selected_interval = gr_intervals_with_overlap.iloc[0, :4]
            ### Select the best combination of graph regions (eg. (interval_1, interval_2) or (interval_1)
            else:
                #### List out all possible combinations
                all_combs = []
                for num in range(1, gr_intervals_with_overlap.shape[0] + 1, 1):
                    combs = list(combinations(gr_intervals_with_overlap.iloc[:, :4].itertuples(), num))
                    all_combs.extend(combs)
                #### Compute graph region-graph region overlap
                overlap_metrics = []
                for i, comb in enumerate(all_combs):
                    # Overlap metrics
                    comb_df = pd.DataFrame(comb).iloc[:, 1:]
                    cumulative_size = (comb_df["end"] - comb_df["start"]).sum()
                    merged_bed_obj = pb.BedTool.from_dataframe(comb_df).sort().merge()
                    peripheral_internal_size = merged_bed_obj.total_coverage()
                    peripheral_internal_overlap = cumulative_size - peripheral_internal_size
                    central_interval_overlap = qr_bed_obj.intersect(merged_bed_obj, wo=True).to_dataframe().iloc[:, -1].sum()
                    overlap_tuple = (i, central_interval_overlap, peripheral_internal_size, peripheral_internal_overlap)
                    overlap_metrics.append(overlap_tuple)
                overlap_df = pd.DataFrame(overlap_metrics)
                overlap_df.columns = ["peripheral_interval", "central_interval_overlap", "peripheral_interval_size", "peripheral_internal_overlap"]
                overlap_df = overlap_df.sort_values(by=["central_interval_overlap", "peripheral_interval_size", "peripheral_internal_overlap"], 
                                                   ascending=[False, True, True]).reset_index(drop=True)

                overlap_df["final_sort_index"] = overlap_df["central_interval_overlap"].astype(float) * 100 - overlap_df["peripheral_interval_size"].astype(float) * 10 - overlap_df["peripheral_internal_overlap"].astype(float) * 1
                overlap_df = overlap_df.sort_values(by=["final_sort_index"], ascending=False).reset_index(drop=True)
                selected_comb = overlap_df.iat[0,0]
                # Make sure all selected intervals are pd.Series objects
                selected_interval = pd.DataFrame([ interval for interval in all_combs[selected_comb] ]).iloc[0, 1:]
                    
        selected_interval.index = ["chr", "start", "end", "strand"]  
        selected_row = gr[gr[['chr', 'start', 'end', 'strand']].apply(lambda row: all(row == selected_interval), axis=1)]
        selected_rows = pd.concat([selected_rows, selected_row], axis=0)
        
    query_nodes = selected_rows.loc[:, ["chr", "start", "end", "strand"]].drop_duplicates()

    # for debugging
    query_nodes_bed_obj = pb.BedTool.from_dataframe(query_nodes)
    query_nodes_bed_obj.intersect(pb.BedTool("selected_regions.bed"), wo=True).saveas("query_nodes.bed")
    logger.info("Writing query nodes to query_nodes.bed")
    logger.info("After picking up merged SD intervals, obtained {} merged SDs from graph with coverage of {} bp.".format(query_nodes.shape[0], query_nodes_bed_obj.sort().total_coverage()))

        
    return selected_rows, query_nodes

