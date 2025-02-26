import argparse
import logging
import os
from multiprocessing import Pool
import sys

import pandas as pd
import pybedtools as pb
import networkx as nx
import graph_tool.all

from preparation import *

from src.log import logger, configure_logger
from src.suppress_warning import *
from src.utils import is_file_up_to_date, filter_bed_by_interval_size
from src.insert_size import get_insert_size_distribution
from src.const import SDrecallPaths

def preparation(paths: SDrecallPaths,
                mq_threshold=41,
                high_quality_depth=10,
                minimum_depth=3,
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
    _, avg_frag_size,std_frag_size = get_insert_size_distribution(input_bam)
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")

    
    # Step 1 : Pick the multi_aligned regions within the target regions
    if not os.path.exists(multi_align_bed) or not is_file_up_to_date(multi_align_bed, [input_bam, target_bed]):
        multi_align_bed = pick_multialigned_regions(input_bam, 
                                                    multi_align_bed, 
                                                    MQ_threshold=mq_threshold, 
                                                    high_quality_depth=high_quality_depth, 
                                                    minimum_depth=minimum_depth, 
                                                    multialign_frac=multialign_frac,
                                                    target_region=target_bed, 
                                                    genome_file=ref_genome )

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
    logger.debug("Raw SD map contains {} regions. (saved to {})".format(total_bin_sd_df.shape[0], os.path.join(outdir + "raw_SD_binary_map.tsv")))
    
    # Convert the both negative strand SDs to both positive strand SDs
    both_neg_strand = (total_bin_sd_df.loc[:, "strand1"] == "-") & (total_bin_sd_df.loc[:, "strand2"] == "-")
    reverse_strand_dict = {"-": "+", "+": "-"}
    total_bin_sd_df.loc[both_neg_strand, "strand1"] = total_bin_sd_df.loc[both_neg_strand, "strand1"].map(reverse_strand_dict)
    total_bin_sd_df.loc[both_neg_strand, "strand2"] = total_bin_sd_df.loc[both_neg_strand, "strand2"].map(reverse_strand_dict)
                                                           
    ## One multi-align interval (represented as *_bam1) can overlap multiple SDs, keep only the minimal set of SDs with sufficient overlap
    by_bam_region = [ g for _, g in total_bin_sd_df.groupby(["chr_bam1", "start_bam1", "end_bam1"], as_index=False) ]
    with Pool(threads) as pool:
        results = pool.imap_unordered(filter_umbrella_pairs, by_bam_region)
        total_bin_sd_df = pd.concat(results, axis=0, ignore_index=True).loc[:, ["chr_1", "start_1", "end_1", "strand1",
                                                                                "chr_2", "start_2", "end_2", "strand2", 
                                                                                "mismatch_rate"]].drop_duplicates().dropna()
    
    expanded_total_bin_sd_df = total_bin_sd_df.copy()
    logger.debug("Total SD map contains {} regions:\n{}".format(total_bin_sd_df.shape[0], total_bin_sd_df.head(5).to_string(index=False)))
    
    ## Remove the duplicated SD combinations (i.e. same SD pair in different order)
    total_bin_sd_df.loc[:, "frozenset_indx"] = total_bin_sd_df.apply(lambda row: frozenset({ row["chr_1"]+":"+str(row["start_1"])+"-"+str(row["end_1"])+":"+row["strand1"], row["chr_2"]+":"+str(row["start_2"])+"-"+str(row["end_2"])+":"+row["strand2"]}), axis=1)
    total_bin_sd_df = total_bin_sd_df.drop_duplicates(subset="frozenset_indx").drop(columns=["frozenset_indx"])
    total_bin_sd_df.to_csv(os.path.join(outdir, "filtered_SD_binary_map.tsv"), sep="\t", index=False)
    logger.info(f"After intersecting with refined SD reference, the total targeted SD region size for sample {input_bam} is {total_bin_sd_bed.sort().total_coverage()}bp.")
    logger.debug(f"Final SD map has shape {total_bin_sd_df.shape} and looks like: \n{total_bin_sd_df.head(5).to_string(index=False)}")
    
    # Step 4: Create a multiplex graph with SD and PO edges
    if os.path.exists(graph_path) and is_file_up_to_date(graph_path, [input_bam, target_bed]):
        graph = read_graphml(graph_path)
    else:
        graph = create_multiplex_graph( total_bin_sd_df, graph_path, threads = threads)

    # Step 4: Extract SD + PO BED regions with depth issues caused by mapping ambiguity. Mutual overlap exists between intervals.
    ## expanded_total_bin_sd_df is constructed from reference SD map, multi_align_bed_obj is constructed from specified target regions.
    ## Filter out regions with length > 140bp (compared to 150bp reads)
    graph_bed = pb.BedTool.from_dataframe(expanded_total_bin_sd_df)
    graph_bed = filter_bed_by_interval_size(graph_bed, 140)
    
    ## Extract SD + PO BED regions with depth issues caused by mapping ambiguity
    intersect_df = graph_bed.intersect(multi_align_bed.merge(), wo=True).to_dataframe(disable_auto_names=True, names=["chrA", "startA", "endA", "strandA", 
                                                                                                                      "chrB", "startB", "endB", "strandB",
                                                                                                                      "mismatch_rate",
                                                                                                                      "chr_target", "start_target", "end_target", "overlap_len"])
    ## Sanity check
    intersect_df = intersect_df.loc[(intersect_df["endA"] > intersect_df["startA"]) & \
                                    (intersect_df["end_target"] > intersect_df["start_target"]) & \
                                    (intersect_df["endB"] > intersect_df["startB"]), :]

    query_groups = [ g for _, g in intersect_df.groupby(["chr_target", "start_target", "end_target"], as_index=False) ]
    
    logger.debug(f"Before filtering, {intersect_df.shape[0]} intervals overlaps with the target regions: \n{intersect_df[:20].to_string(index=False)}\n")
    with Pool(threads) as pool:
        results = pool.imap_unordered(filter_umbrella_pairs, query_groups)
        intersect_df = pd.concat(results, axis=0, ignore_index=True).drop_duplicates()
    
    logger.debug("Only {} non-umbrella intervals overlaps target region: \n{}\n".format(intersect_df.shape[0], intersect_df[:20].to_string(index=False)))
    
    query_nodes = intersect_df.loc[:, ["chrA", "startA", "endA", "strandA"]].drop_duplicates().rename(columns={"chrA":"chr", "startA":"start", "endA":"end", "strandA":"strand"})
    if logger.level == logging.DEBUG:
        final_query_bed = pb.BedTool.from_dataframe(query_nodes)
        final_query_bed.intersect(multi_align_bed, wo=True).saveas(paths.target_overlapping_query_sd_bed_path()) # renamed from coding_query_nodes.bed
        logger.debug("Target-overlapping SD regions saved to: {}".format(paths.target_overlapping_query_sd_bed_path()))
    logger.info("Identified {} distinct target-overlapping SDs from graph".format(query_nodes.shape[0]))

    # Step 5. Pool all SD nodes and the corresponding paralogous regions from the SD + PO graph
    final_graph_path = paths.annotated_graph_path()
    sd_paralog_pairs, connected_qnode_components = extract_SD_paralog_pairs_from_graph( query_nodes, 
                                                                                        graph, 
                                                                                        graph_path = final_graph_path, 
                                                                                        reference_fasta = ref_genome,
                                                                                        avg_frag_size = avg_frag_size, 
                                                                                        std_frag_size = std_frag_size, 
                                                                                        threads=threads )
    logger.info("{} SD-paralog pairs pooled from SD + PO graph".format(len(sd_paralog_pairs)))
    
    # Step 6: Collapse qnodes by identifying qnodes that can be put in the same masked genome
    grouped_qnode_cnodes = query_connected_nodes(sd_paralog_pairs, connected_qnode_components)

    # # (Debug code) Enumerate qnodes and cnodes in the graph. (FC -> qnode; NFC -> cnode)
    # final_graph = graph.copy()
    # for i, result in enumerate(grouped_qnode_cnodes):
    #     tag = "PC" + str(i)
    #     for node in result["SD_qnodes"]:
    #         if not isinstance(node, tuple):
    #             raise TypeError(f"Expected {node} (type: {type(node)}) to be tuple.")
    #         final_graph.nodes[node]["FC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("FC", "").split(",") + [tag])) if len(e) > 0])
    #     for node in result["SD_counterparts"]:
    #         if not isinstance(node, HOMOSEQ_REGION):
    #             raise TypeError(f"Expected {node} (type: {type(node)}) to be HOMOSEQ_REGION. ")
    #         final_graph.nodes[node.data]["NFC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node.data].get("NFC", "").split(",") + [tag])) if len(e) > 0])
    # nx.write_graphml(final_graph, final_graph_path)
    
    # Step 7: Create beds and masked genomes
    build_beds_and_masked_genomes(grouped_qnode_cnodes = grouped_qnode_cnodes,
                                  sd_paralog_pairs = sd_paralog_pairs,
                                  output_folder = outdir,
                                  ref_genome = ref_genome,
                                  sdrecall_paths = paths,
                                  nthreads = threads)

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
        base_dir=args.outdir
    )
    
    # Log the derived values for confirmation
    logger.info(f"Running with sample ID: {paths.sample_id}")
    logger.info(f"Target tag: {paths.target_tag}")
    logger.info(f"Assembly: {paths.assembly}")
    logger.info(f"Working directory: {paths.work_dir}")

    # Call preparation function with the initialized paths
    preparation(
        paths=paths,  # Pass the paths instance to the preparation function
        mq_threshold=args.mq_cutoff,
        high_quality_depth=args.high_quality_depth,
        minimum_depth=args.minimum_depth,
        multialign_frac=args.multialign_frac,
        threads=args.threads,
    )

if __name__ == "__main__":
    main()
