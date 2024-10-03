import gc
import logging
from multiprocessing import Pool
import pandas as pd

import networkx as nx
import graph_tool.all
from intervaltree import Interval, IntervalTree
from pybedtools import BedTool

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
    
def compose_multiplex_graph_per_chr(args):
    chr_id, sd_data_chr, overlap_frac_min, avg_frag_size, std_frag_size, resolution, graph_file_path = args
    sd_interval_tree = IntervalTree()
    directed_overlap_graph = nx.DiGraph()
    directed_graph_path = graph_file_path.replace(".graphml", f".directed.overlap.{chr_id}.graphml")
    assert directed_graph_path != graph_file_path, "The directed graph path is the same as the graph path, please check the input graph path"
    
    # Only compose overlap graph per chromosome
    for i in range(0, len(sd_data_chr)):
        row = sd_data_chr.iloc[i, :]
        directed_overlap_graph.add_node((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]), size = row["end_1"] - row["start_1"])
        directed_overlap_graph.add_node((row["chr_2"], row["start_2"], row["end_2"], row["strand2"]), size = row["end_2"] - row["start_2"])
        interval_1 = Interval(row["start_1"], row["end_1"], (row["chr_1"], row["start_1"], row["end_1"], row["strand1"]))
        interval_2 = Interval(row["start_2"], row["end_2"], (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]))
        for interval in [interval_1, interval_2]:
            if chr_id == interval.data[0]: sd_interval_tree.add(interval)
            for overlap in sd_interval_tree.overlap(interval.begin, interval.end):
                overlap_span = min([interval.end, overlap.end]) - max([interval.begin, overlap.begin])
                overlap_size = overlap.end - overlap.begin
                interval_size = interval.end - interval.begin
                weight = max(overlap_span/overlap_size, overlap_span/interval_size)
                uncovered_size = min(overlap_size, interval_size) - overlap_span
                if interval_size >= overlap_size and interval.data[0] == overlap.data[0] and overlap.data != interval.data:
                    directed_overlap_graph.add_edge(interval.data, overlap.data, weight=1/weight, overlap="True", uncovered_smaller_interval=uncovered_size)
                elif interval.data[0] == overlap.data[0] and overlap.data != interval.data:
                    directed_overlap_graph.add_edge(overlap.data, interval.data, weight=1/weight, overlap="True", uncovered_smaller_interval=uncovered_size)

    # Now we would like to save the directed graph as a graphml file and see the results in cytoscape
    nx.write_graphml(directed_overlap_graph, directed_graph_path)
    gc.collect()
    return (chr_id, directed_graph_path)

def create_multiplex_graph(sd_data, graph_filepath=None, 
                            overlap_frac_min = 0.8, 
                            avg_frag_size=350, 
                            std_frag_size = 120,
                            resolution = .1,
                            threads = 10, 
                            target_bed = "",
                            base_dir = ""):
    
    # sd_data is the binary sd map of genomic intervals
    # Input the sd_data to networkx Graph object
    G = nx.Graph()
    all_chrs = set(sd_data["chr_1"]).union(set(sd_data["chr_2"]))
    
    if os.path.exists(target_bed):
        target_chrs = pd.read_table(target_bed, header=None).iloc[:, 0].drop_duplicates().tolist()
        all_chrs = [c for c in all_chrs if c in target_chrs]

    logger.info(f"The input sd_data dataframe looks like: \n{sd_data[:10].to_string(index=False)}\n")
    # adding edges between nodes in the same line
    # The below for loop not just build up a Graph composed of intervals sharing sequence similarity with each other
    # but also builds up an interval tree for each chromosome using the same intervals
    sd_data_by_chr = {chr_id: sd_data.loc[(sd_data["chr_1"]==chr_id) | (sd_data["chr_2"]==chr_id), sd_data.columns.tolist()[:8]] for chr_id in all_chrs}
    for i in range(0, len(sd_data)):
        row = sd_data.iloc[i, :]
        # First adding the binary pairwise SDs to the graph
        G.add_edge((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]),
                   (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]),
                    type = "segmental_duplication", weight=float(row["mismatch_rate"]))
    logger.info(f"How many nodes in the graph: {G.number_of_nodes()} How many edges are SD edges: {G.number_of_edges()} ")
    
    with Pool(threads) as pool:
        # Compose overlap graph within each chromsome
        partitions_by_chr = dict(pool.imap_unordered(compose_multiplex_graph_per_chr, 
                                                    ((chr_id, sd_data_by_chr[chr_id], overlap_frac_min, avg_frag_size, std_frag_size, resolution, graph_filepath) for chr_id in all_chrs)))

    directed_graphs = [read_graphml(dg_path) for chr_id, dg_path in partitions_by_chr.items()]  # k is node and v is the directed graph path
    # Merge the directed graph into one graph
    merged_directed_graph = nx.DiGraph()
    for dg in directed_graphs:
        merged_directed_graph = nx.compose(merged_directed_graph, dg)
        
    # Then we add edges to merged_directed_graph based on the edges in the sd_graph networks which have their type attributes set to "segmental_duplication"
    for u, v, d in G.edges(data=True):
        if d.get("type", None) == "segmental_duplication":
            # Test if u and v are in the directed graph nodes
            if u in merged_directed_graph.nodes() and v in merged_directed_graph.nodes():
                # Test if these two nodes already has an edge between them
                if not merged_directed_graph.has_edge(u, v):
                    merged_directed_graph.add_edge(u, v, **d)
                else:
                    # Just add some attributes to the existing edge instead of force cover the original edge's attributes
                    merged_directed_graph[u][v]["type"] = "segmental_duplication"
            else:
                logger.warning("These two genomic regions sharing sequences seem not contained by the overlap graph: \n{}\n{}\n".format(u, v))
                
    # Now we add some annotations to the directed graph regarding which transcript is overlapped with each node.
    # First convert the graph to bed
    dg_bed = graph_to_sorted_bedtool(merged_directed_graph)
    logger.info(f"The graph converted bed table looks like: \n{dg_bed.to_dataframe()[:5].to_string(index=False)}\n")
    target_bed_obj = BedTool(target_bed)
    dg_exon_df = dg_bed.intersect(target_bed_obj, wo = True).to_dataframe(disable_auto_names = True, 
                                                                          names =["chr_graph", "start_graph", "end_graph", "graph_strand", 
                                                                                  "chr_feature", "start_feature", "end_feature",
                                                                                  "symbol", "tranxID", 
                                                                                  "exon_No", "strand",
                                                                                  "overlap_len"])
    logger.info("Take a look at the overlap dataframe here:\n{}\n".format(dg_exon_df[:5].to_string(index=False)))
    dg_node_overlap_dict = {}
    by_graph_interval = dg_exon_df.groupby(by=["chr_graph", "start_graph", "end_graph", "graph_strand"], as_index=False)
    graph_groups = [group for _, group in by_graph_interval]
    with Pool(threads) as pool:
        results = pool.imap_unordered(summary_exon_overlap, graph_groups)
        summary_df = pd.concat(results, ignore_index=True).drop_duplicates()
    logger.info(f"The summary exon overlap function applied by groupby returned an object of {type(summary_df)}, take a look at the df:\n{summary_df[:5].to_string()}\n")
    for ind, row in summary_df.iterrows():
        node = (row["chr_graph"], int(row["start_graph"]), int(row["end_graph"]), row["graph_strand"])
        dg_node_overlap_dict[node] = row[-1]
    
    for node, overlap_descript in dg_node_overlap_dict.items():
        if node in merged_directed_graph.nodes():
            merged_directed_graph.nodes[node]["FC_overlap"] = overlap_descript
    
    sd_edge_count = [e.get("type", None) for u,v,e in merged_directed_graph.edges(data=True)].count("segmental_duplication")
    logger.info(f"How many nodes are there in the total SD graph with edges implying overlaps and homology? {merged_directed_graph.number_of_nodes()} How many edges in this graph: {merged_directed_graph.number_of_edges()}. How many of them belonged to type segmental_duplication: {sd_edge_count} ")
                    
    if graph_filepath is not None:
        nx.write_graphml(merged_directed_graph, graph_filepath)
        logger.info("The merged directed graph is now saved as {}".format(graph_filepath))
        
    return merged_directed_graph

def read_graphml(graph_path):
    literal_graph = nx.read_graphml(graph_path)
    # We need to convert the string back to tuple
    graph = literal_graph.__class__()
    for node, data in literal_graph.nodes(data=True):
        graph.add_node(string_to_tuple(node), **data)

    for node1, node2, data in literal_graph.edges(data=True):
        graph.add_edge(string_to_tuple(node1), string_to_tuple(node2), **data)
        
    return graph

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