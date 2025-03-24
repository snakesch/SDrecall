import argparse
import logging
import os
from multiprocessing import Pool

import pandas as pd
import pybedtools as pb

from preparation import *

from src.log import logger, configure_logger
from src.suppress_warning import *
from src.utils import is_file_up_to_date, executeCmd

from src.const import SDrecallPaths


def save_connected_qnodes_to_graphml(g, output_path, logger = logger):
    """
    Save connected qnodes graph to GraphML with readable node properties
    
    Args:
        g: graph-tool Graph object with node_data property
        output_path: path to save GraphML file
    """
    # Create separate properties for each component of the node_data tuple
    v_chr = g.new_vertex_property("string")
    v_start = g.new_vertex_property("int")
    v_end = g.new_vertex_property("int")
    v_strand = g.new_vertex_property("string")
    
    # Extract components from node_data property
    node_data_prop = g.vertex_properties["node_data"]
    for v in g.vertices():
        node_tuple = node_data_prop[v]
        if node_tuple:
            chr_val, start_val, end_val, strand_val = node_tuple
            v_chr[v] = str(chr_val)
            v_start[v] = int(start_val)
            v_end[v] = int(end_val)
            v_strand[v] = str(strand_val)
    
    # Add the new properties to the graph
    g.vertex_properties["chr"] = v_chr
    g.vertex_properties["start"] = v_start
    g.vertex_properties["end"] = v_end
    g.vertex_properties["strand"] = v_strand
    
    # Write the graph to GraphML file
    # Include the original node_data property for compatibility
    g.save(output_path, fmt="graphml")
    logger.info(f"Connected qnodes graph saved to {output_path}. It contains {g.num_vertices()} qnodes and {g.num_edges()} edges.")
    return output_path


def prepare_recall_regions( paths: SDrecallPaths,
                            mq_threshold=41,
                            high_quality_depth=10,
                            minimum_depth=5,
                            multialign_frac=0.5,
                            threads=10):
    """
    Main preparation function for SDrecall
    
    Args:
        paths: Initialized SDrecallPaths instance
        mq_threshold: Mapping quality threshold
        high_quality_depth: High quality depth cutoff
        minimum_depth: Minimum depth cutoff
        multialign_frac: Multi-align fraction cutoff
        threads: Number of threads to use
    """
    pb.helpers.set_tempdir(paths.tmp_dir)
    # Access all paths through the paths object
    ref_genome = paths.ref_genome
    input_bam = paths.input_bam
    target_bed = paths.target_bed
    
    # Use path methods for all file access
    multi_align_bed = paths.multi_align_bed_path()
    graph_path = paths.multiplex_graph_path()
    
    # Make sure output directory exists
    outdir = paths.work_dir
    os.makedirs(outdir, exist_ok=True)

    # Step 0: Calculate the distribution of fragment sizes
    avg_frag_size, std_frag_size = paths.median_frag_size, paths.frag_size_std
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")

    
    # Step 1 : Pick the multi_aligned regions within the target regions
    pick_multialign_regions = True
    if os.path.exists(multi_align_bed) and \
       is_file_up_to_date(multi_align_bed, [input_bam, target_bed]) and \
       os.path.getsize(multi_align_bed) > 10:
        logger.info(f"Seems the multialign bed has already been created, skipping the multialign bed creation step.")
        executeCmd(f"ls -lh {multi_align_bed}", logger= logger)
        pick_multialign_regions = False

    if pick_multialign_regions:
        multi_align_bed = pick_multialigned_regions(input_bam, 
                                                    multi_align_bed, 
                                                    MQ_threshold=mq_threshold, 
                                                    high_quality_depth=high_quality_depth, 
                                                    minimum_depth=minimum_depth, 
                                                    multialign_frac=multialign_frac,
                                                    target_region=target_bed, 
                                                    target_tag=paths.target_tag,
                                                    genome_file=paths.ref_genome_fai_path())

    multi_align_bed_obj = pb.BedTool(multi_align_bed).sort()
    logger.info("The multialign bed extracted from {} is {} and it covers {} bp.".format(input_bam,
                                                                                         multi_align_bed,
                                                                                         multi_align_bed_obj.total_coverage()))
    

    # Step 2. Load reference SD map; trim the SD coordinates and filter by size
    # The reference SD map is already expanded, either interval of a pair of SDs are present in the first 3 columns once.
    ref_bed_obj = pb.BedTool(paths.reference_sd_map).sort()
    logger.info("Total coverage of reference SD map: {}bp".format(ref_bed_obj.total_coverage()))

    ## Filter out SD regions smaller than average insert size.
    big_ref_bed_obj = ref_bed_obj.filter(lambda x: len(x) > avg_frag_size).saveas() ## .saveas makes the BedTool object persistent in memory
    logger.info(f"After excluding paired reference SDs with fragment size < {avg_frag_size:.1f}bp, the reference SD map covers {big_ref_bed_obj.total_coverage()}bp. ")

    # Step 3: Compare multialigned BED regions to the reference SD map to identify SD regions with mapping ambiguity
    total_bin_sd_bed = big_ref_bed_obj.intersect(multi_align_bed, wo=True)
    logger.info(f"After overlapping with refined SD reference, {total_bin_sd_bed.sort().total_coverage()}bp SD regions need variant recall.")
    logger.debug("Overlapping SD map: \n{}".format(total_bin_sd_bed.to_dataframe()[:10].to_string(index=False)))
    total_bin_sd_df = total_bin_sd_bed.to_dataframe(disable_auto_names=True, 
                                                    names = ["chr_1", "start_1", "end_1",
                                                             "chr_2", "start_2", "end_2",
                                                             "strand1", "strand2",
                                                             "cigar", "mismatch_rate",
                                                             "chr_bam1", "start_bam1", "end_bam1",
                                                             "overlap_len"]).loc[:,["chr_1", "start_1", "end_1", "strand1",
                                                                                    "chr_2", "start_2", "end_2", "strand2",
                                                                                    "chr_bam1", "start_bam1", "end_bam1", "mismatch_rate", "overlap_len"]].drop_duplicates().dropna()
    total_bin_sd_df.sort_values(by=["chr_bam1", "start_bam1", "end_bam1"], inplace=True)
    # Filter out the SD pairs with either one of them on alternative contigs
    # Need to prepare the regex for main contig names first
    main_contigs_pattern = r'^(chr)?([0-9]+|[XYM]|MT)$'
    # Filter to keep only SD pairs where both contigs are main chromosomes
    main_contig_mask = (
        total_bin_sd_df['chr_1'].str.match(main_contigs_pattern, case=False) & \
        total_bin_sd_df['chr_2'].str.match(main_contigs_pattern, case=False) & \
        total_bin_sd_df['chr_bam1'].str.match(main_contigs_pattern, case=False)
    )
    total_bin_sd_df = total_bin_sd_df[main_contig_mask].reset_index(drop=True)
    logger.info(f"After filtering out alternative contigs, {total_bin_sd_df.shape[0]} SD pairs remain")
    
    # Convert the both negative strand SDs to both positive strand SDs
    both_neg_strand = (total_bin_sd_df.loc[:, "strand1"] == "-") & (total_bin_sd_df.loc[:, "strand2"] == "-")
    reverse_strand_dict = {"-": "+", "+": "-"}
    total_bin_sd_df.loc[both_neg_strand, "strand1"] = total_bin_sd_df.loc[both_neg_strand, "strand1"].map(reverse_strand_dict)
    total_bin_sd_df.loc[both_neg_strand, "strand2"] = total_bin_sd_df.loc[both_neg_strand, "strand2"].map(reverse_strand_dict)
                                                           
    ## One multi-align interval (represented as *_bam1) can overlap multiple SDs, keep only the minimal set of SDs with sufficient overlap
    logger.info(f"Before filtering, the total SD map has {total_bin_sd_df.shape[0]} SDs. Looks like: \n{total_bin_sd_df.head(10).to_string(index=False)}")
    by_bam_region = [ g for _, g in total_bin_sd_df.groupby(["chr_bam1", "start_bam1", "end_bam1"], as_index=False) ]
    with Pool(threads) as pool:
        results = pool.imap_unordered(filter_umbrella_pairs, by_bam_region)
        total_bin_sd_df = pd.concat(results, axis=0, ignore_index=True).drop_duplicates().dropna().sort_values(by=["chr_bam1", "start_bam1", "end_bam1"])
    logger.info("After filtering umbrella SD pairs, the total SD map contains {} SD pairs:\n{}".format(total_bin_sd_df.shape[0], total_bin_sd_df.head(10).to_string(index=False)))

    total_bin_sd_df = total_bin_sd_df.loc[:, ["chr_1", "start_1", "end_1", "strand1",
                                              "chr_2", "start_2", "end_2", "strand2", 
                                              "mismatch_rate"]].drop_duplicates().dropna()
    logger.info(f"Now there are only {total_bin_sd_df.shape[0]} SD regions left.")
    query_nodes = total_bin_sd_df.loc[:, ["chr_1", "start_1", "end_1", "strand1"]].drop_duplicates().rename(columns={"chr_1":"chr", "start_1":"start", "end_1":"end", "strand1":"strand"})
    
    ## Remove the duplicated SD combinations (i.e. same SD pair in different order) for graph construction
    total_bin_sd_df.loc[:, "frozenset_indx"] = total_bin_sd_df.apply(lambda row: frozenset({ row["chr_1"]+":"+str(row["start_1"])+"-"+str(row["end_1"])+":"+row["strand1"], row["chr_2"]+":"+str(row["start_2"])+"-"+str(row["end_2"])+":"+row["strand2"]}), axis=1)
    total_bin_sd_df = total_bin_sd_df.drop_duplicates(subset="frozenset_indx").drop(columns=["frozenset_indx"])
    total_bin_sd_df.to_csv(os.path.join(outdir, "filtered_SD_binary_map.tsv"), sep="\t", index=False)
    logger.info(f"After intersecting with refined SD reference, the total targeted SD region size for sample {input_bam} is {total_bin_sd_bed.sort().total_coverage()}bp.")
    logger.debug(f"Final SD map has shape {total_bin_sd_df.shape} and looks like: \n{total_bin_sd_df.head(5).to_string(index=False)}")
    
    # Step 4: Create a multiplex graph with SD and PO edges
    if os.path.exists(graph_path) and \
       is_file_up_to_date(graph_path, [input_bam, target_bed]) and \
       os.path.getsize(graph_path) > 10:
        logger.info(f"Seems the graph has already been created, skipping the graph creation step.")
        executeCmd(f"ls -lh {graph_path}", logger= logger)
        graph = read_graphml(graph_path)
    else:
        graph = create_multiplex_graph( total_bin_sd_df, graph_path, threads = threads)
    
    if logger.level == logging.DEBUG:
        final_query_bed = pb.BedTool.from_dataframe(query_nodes)
        final_query_bed.intersect(multi_align_bed, wo=True).saveas(paths.target_overlapping_query_sd_bed_path()) # renamed from coding_query_nodes.bed
        logger.debug("Target-overlapping SD regions saved to: {}".format(paths.target_overlapping_query_sd_bed_path()))
    logger.info("Identified {} distinct target-overlapping SDs from graph".format(query_nodes.shape[0]))

    # Step 5. Pool all SD nodes and the corresponding paralogous regions from the SD + PO graph
    final_graph_path = paths.annotated_graph_path()
    sd_paralog_pairs, connected_qnode_graph = extract_SD_paralog_pairs_from_graph( query_nodes, 
                                                                                    graph, 
                                                                                    graph_path = final_graph_path, 
                                                                                    reference_fasta = ref_genome,
                                                                                    tmp_dir = paths.tmp_dir,
                                                                                    avg_frag_size = avg_frag_size, 
                                                                                    std_frag_size = std_frag_size, 
                                                                                    threads=threads )
    # Save the connected qnodes graph (graph-tool graph) to a GRAPHML file
    connected_qnodes_graph_path = paths.qnode_grouping_graph()
    save_connected_qnodes_to_graphml(connected_qnode_graph, connected_qnodes_graph_path, logger = logger)
    logger.info("{} SD-paralog pairs pooled from SD + PO graph. The graph recording similar query nodes is saved to {}".format(len(sd_paralog_pairs), connected_qnodes_graph_path))
    
    # Step 6: Collapse qnodes by identifying qnodes that can be put in the same masked genome
    grouped_qnode_cnodes = query_connected_nodes(sd_paralog_pairs, 
                                                 connected_qnode_graph)

    '''
    (Debug code) Enumerate qnodes and cnodes in the graph. (FC -> qnode; NFC -> cnode)
    This is trying to output a GRAPHML file with robust annotations for debugging purposes.
    final_graph = graph.copy()
    for i, result in enumerate(grouped_qnode_cnodes):
        tag = "PC" + str(i)
        for node in result["SD_qnodes"]:
            if not isinstance(node, tuple):
                raise TypeError(f"Expected {node} (type: {type(node)}) to be tuple.")
            final_graph.nodes[node]["FC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("FC", "").split(",") + [tag])) if len(e) > 0])
        for node in result["SD_counterparts"]:
            if not isinstance(node, HOMOSEQ_REGION):
                raise TypeError(f"Expected {node} (type: {type(node)}) to be HOMOSEQ_REGION. ")
            final_graph.nodes[node.data]["NFC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node.data].get("NFC", "").split(",") + [tag])) if len(e) > 0])
    nx.write_graphml(final_graph, final_graph_path)
    '''
    
    # Step 7: Create beds and masked genomes
    paths = build_beds_and_masked_genomes(grouped_qnode_cnodes = grouped_qnode_cnodes,
                                          sd_paralog_pairs = sd_paralog_pairs,
                                          sdrecall_paths = paths,
                                          nthreads = threads,
                                          avg_frag_size = avg_frag_size,
                                          std_frag_size = std_frag_size)
    logger.info(f"SDrecall has been prepared for {paths.sample_id}, targeting regions specified by {paths.target_bed}.")
    return paths


def main():
    parser = argparse.ArgumentParser(description='Deploy SD_qnodes for SDrecall.')

    parser.add_argument('-r', '--ref_genome', required=True, help='Path to the reference genome. Currently only accept hg19 and hg38')
    parser.add_argument('-o', '--outdir', required=True, help='Base directory for output files.')
    parser.add_argument('-i', '--input_bam', required=True, help='Input BAM file.')
    parser.add_argument('-m', '--reference_sd_map', required=True, help='Reference SD map file.')
    parser.add_argument('-b', '--target_bed', default="", help='Optional target BED file.')
    parser.add_argument('-s', '--sample_id', default=None, help='Sample ID (default: extracted from BAM filename)')
    parser.add_argument('--target_tag', type=str, default=None, help='Optional target tag for filtering.')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads to use.')
    parser.add_argument('--mq_cutoff', type=int, default=41, help='Mapping quality cutoff.')
    parser.add_argument('--high_quality_depth', type=int, default=10, help='High quality depth cutoff.')
    parser.add_argument('--minimum_depth', type=int, default=3, help='Minimum depth cutoff.')
    parser.add_argument('--multialign_frac', type=float, default=0.5, help='Multi-align fraction cutoff.')
    parser.add_argument('-v', '--verbose', type=str, default="INFO", help='Level of verbosity (default = INFO).')

    args = parser.parse_args()

    configure_logger(log_level=args.verbose)
    
    # Initialize the SDrecallPaths singleton instance
    paths = SDrecallPaths.initialize(
        ref_genome=args.ref_genome,
        input_bam=args.input_bam,
        reference_sd_map=args.reference_sd_map,
        target_bed=args.target_bed,
        sample_id=args.sample_id,  # Use user-provided sample_id if available
        target_tag=args.target_tag,  # Use user-provided target_tag if available 
        output_dir=args.outdir
    )
    
    # Log the derived values for confirmation
    logger.info(f"Running with sample ID: {paths.sample_id}")
    logger.info(f"Target tag: {paths.target_tag}")
    logger.info(f"Assembly: {paths.assembly}")
    logger.info(f"Working directory: {paths.work_dir}")

    # Call preparation function with the initialized paths
    prepare_recall_regions(
        paths=paths,  # Pass the paths instance to the preparation function
        mq_threshold=args.mq_cutoff,
        high_quality_depth=args.high_quality_depth,
        minimum_depth=args.minimum_depth,
        multialign_frac=args.multialign_frac,
        threads=args.threads,
    )

if __name__ == "__main__":
    main()
