import gc
import sys
import logging
from multiprocessing import Pool
import pandas as pd

import networkx as nx
import graph_tool.all as gt
from intervaltree import Interval, IntervalTree
from pybedtools import BedTool

from src.utils import string_to_tuple

logger = logging.getLogger("SDrecall")

class HashableGraph:
    def __init__(self, graph):
        self.graph = graph

    def __hash__(self):
        node_data = sorted(self.graph.nodes(data=True))
        edge_data = sorted(self.graph.edges(data=True))
        return hash((tuple(node_data), tuple(edge_data)))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __getattr__(self, name):
        return getattr(self.graph, name)
    
    def to_graph(self):
        return self.graph

def graph_to_sorted_bedtool(G):
    # Extract nodes from the graph and convert them into bed format
    bed_format = ["\t".join(map(str, node)) for node in G.nodes]

    # Create a BedTool object
    bedtool = BedTool("\n".join(bed_format), from_string=True).sort()

    return bedtool
    
def compose_PO_graph_per_chr(args):
    chrom, sd_data, graph_file_path = args
    sd_data_chr = sd_data.loc[(sd_data["chr_1"] == chrom) | (sd_data["chr_2"] == chrom), sd_data.columns.tolist()[:8]]

    sd_interval_tree = IntervalTree()
    directed_overlap_graph = nx.DiGraph()
    directed_graph_path = graph_file_path.replace(".graphml", f".directed.overlap.{chrom}.graphml")

    if directed_graph_path == graph_file_path:
        logger.error("The directed graph path is the same as the graph path, please check the input graph path.")
        sys.exit(1)
    
    ## - Properties of PO graph - ##
    ## 1. Potential PO regions are tracked by a growing IntervalTree
    ## 2. PO edges point from large intervals to small intervals
    ## 3. A PO edge is drawn irrespective of extent of overlap
    ## 4. Edge weights are defined as the minimum of interval size divided by overlapping span
    for _, row in sd_data_chr.iterrows():
        directed_overlap_graph.add_node((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]), size = row["end_1"] - row["start_1"])
        directed_overlap_graph.add_node((row["chr_2"], row["start_2"], row["end_2"], row["strand2"]), size = row["end_2"] - row["start_2"])
        interval_1 = Interval(row["start_1"], row["end_1"], (row["chr_1"], row["start_1"], row["end_1"], row["strand1"]))
        interval_2 = Interval(row["start_2"], row["end_2"], (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]))
        
        for interval in [interval_1, interval_2]:
            ## Check if target intervals are located on chrom
            if chrom == interval.data[0]: 
                sd_interval_tree.add(interval)
            ## o_interval represents an interval with overlap with 'interval'
            for o_interval in sd_interval_tree.overlap(interval.begin, interval.end):
                ## Skip if interval and o_interval belong to different chromosomes or represent the same region
                if interval.data[0] != o_interval.data[0] or interval.data == o_interval.data:
                    continue
                overlap_span = min([interval.end, o_interval.end]) - max([interval.begin, o_interval.begin])
                o_size = o_interval.end - o_interval.begin
                interval_size = interval.end - interval.begin
                weight = max(overlap_span/o_size, overlap_span/interval_size)
                nonoverlapping_size = min(o_size, interval_size) - overlap_span
                if interval_size >= o_size:
                    directed_overlap_graph.add_edge(interval.data, o_interval.data, weight=1/weight, overlap="True", nonoverlapping_smaller_interval=nonoverlapping_size)
                else:
                    ## Overlap simply means any intervals overlapping the query interval, not necessarily a subset of all intervals
                    directed_overlap_graph.add_edge(o_interval.data, interval.data, weight=1/weight, overlap="True", nonoverlapping_smaller_interval=nonoverlapping_size)

    # Save the directed graph as a graphml file and visualize the results in cytoscape
    nx.write_graphml(directed_overlap_graph, directed_graph_path)
    gc.collect()
    return (chrom, directed_graph_path)

def create_multiplex_graph(sd_data, graph_filepath=None, threads = 10, target_bed = ""):
    '''
    Inputs:
    sd_data: refined and filtered SD map from previous steps (pd.DataFrame)
    graph_filepath: output graph path (str)

    Returns:
    merged_directed_graph: NetworkX DiGraph with SD edges and PO edges
    '''
    import os

    G = nx.Graph()
    all_chrs = set(sd_data["chr_1"]).union(set(sd_data["chr_2"]))

    ## - Graph construction - ##
    ## 1. Create a graph of nodes with **SD edges**
    ## 2. Create a DiGraph with physical overlap (PO) as edges for **each chromosome**
    ## 3. Merge all PO graphs
    ## 4. In the merged PO graph, draw SD edges between SD nodes or add SD as attribute if PO edge already exists
    for _, row in sd_data.iterrows():
        ## Draw SD edges with mismatch rates as edge weights
        G.add_edge((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]),
                   (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]),
                    type = "segmental_duplication", weight=float(row["mismatch_rate"]))
    logger.debug(f"Constructed the initial graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} SD edges.")
    
    with Pool(threads) as pool:
        # Compose a PO graph for each chromosome
        partitions_by_chr = dict(pool.imap_unordered(compose_PO_graph_per_chr, 
                                                    ((chrom, sd_data, graph_filepath) for chrom in all_chrs)))

    directed_graphs = [read_graphml(dg_path) for dg_path in partitions_by_chr.values()]  # k is node and v is the directed graph path
    
    ## Clean up GraphML files if not in debugging
    if logger.level != logging.DEBUG:
        for dg_path in partitions_by_chr.values():
            os.remove(dg_path)
    
    # Merge the directed PO graphs
    merged_PO_graph = nx.DiGraph()
    for dg in directed_graphs:
        merged_PO_graph = nx.compose(merged_PO_graph, dg)
        
    # Draw SD edges between SD nodes or add SD as attribute if PO edge already exists
    for u, v, d in G.edges(data=True):
        if d.get("type", None) == "segmental_duplication":
            if u in merged_PO_graph.nodes() and v in merged_PO_graph.nodes():
                if not merged_PO_graph.has_edge(u, v):
                    merged_PO_graph.add_edge(u, v, **d)
                else:
                    merged_PO_graph[u][v]["type"] = "segmental_duplication"
            else:
                logger.warning("These two overlapping genomic regions are not found in the overlap graph: \n{}\n{}".format(u, v))
        
    sd_edge_count = [e.get("type", None) for _, _, e in merged_PO_graph.edges(data=True)].count("segmental_duplication")
    logger.info(f"Merged SD + PO graph contains {merged_PO_graph.number_of_nodes()} nodes and {merged_PO_graph.number_of_edges()} edges. ({sd_edge_count} SD edges)")

    if logger.level == logging.DEBUG and graph_filepath is not None:    
        nx.write_graphml(merged_PO_graph, graph_filepath)
        logger.debug("Merged SD + PO graph is saved as {}".format(graph_filepath))
        
    return merged_PO_graph

def read_graphml(graph_path):
    literal_graph = nx.read_graphml(graph_path)

    graph = literal_graph.__class__()
    for node, data in literal_graph.nodes(data=True):
        graph.add_node(string_to_tuple(node), **data)

    for node1, node2, data in literal_graph.edges(data=True):
        graph.add_edge(string_to_tuple(node1), string_to_tuple(node2), **data)
        
    return graph
