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

def preparation( ref_genome: str,
                 work_dir: str, 
                 input_bam: str,
                 reference_sd_map: str,
                 target_bed = "",
                 mq_threshold = 41,
                 high_quality_depth = 10,
                 minimum_depth = 3,
                 multialign_frac = 0.5,
                 threads = 10,
                 target_tag = "target"):

    os.makedirs(work_dir, exist_ok=True)

    rg = Genome(path=ref_genome)
    genome_file = rg.fai_index

    # Step 0: Calculate the distribution of fragment sizes
    _, avg_frag_size,std_frag_size = get_insert_size_distribution(input_bam)
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")
    
    ## Define names of intermediate files. Intermediate / final outputs will be written to work_dir.
    basename = os.path.basename(input_bam).replace(".bam", "")
    multi_align_bed = os.path.join(work_dir, basename + f".{target_tag}.multialign.bed")
    
    # Step 1 : Pick the multi_aligned regions within the target regions
    if not os.path.exists(multi_align_bed) or not is_file_up_to_date(multi_align_bed, [input_bam, target_bed]):
        multi_align_bed = pick_region_by_depth(input_bam, 
                                        multi_align_bed, 
                                        MQ_threshold=mq_threshold, 
                                        high_quality_depth=high_quality_depth, 
                                        minimum_depth=minimum_depth, 
                                        multialign_frac=multialign_frac,
                                        target_region=target_bed, 
                                        genome_file=genome_file)

    multi_align_bed_obj = pb.BedTool(multi_align_bed).sort()
    logger.info("The multialign bed extracted from {} is {} and it covers {} bp.".format(input_bam,
                                                                                      multi_align_bed,
                                                                                      multi_align_bed_obj.total_coverage()))
    
    # Step 2: Compare paired SD regions with the reference SD map to identify SD regions associated with mapping ambiguity
    ref_bed_obj = pb.BedTool(reference_sd_map)
    logger.info("Total coverage of reference SD map: {}bp".format(ref_bed_obj.sort().total_coverage()))

    ## Filter out the pairs where paired intervals are smaller than the avg fragment size.
    big_ref_bed_obj = filter_bed_by_interval_size(ref_bed_obj, avg_frag_size)
    logger.info(f"After excluding paired reference SDs with fragment size < {avg_frag_size}bp, the reference SD map covers {big_ref_bed_obj.total_coverage()}bp. ")
    
    ## Intersect with previously derived sample-specific SD map
    total_bin_sd_bed = big_ref_bed_obj.intersect(multi_align_bed_obj, wo=True)
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
    logger.debug("Raw SD map contains {} regions. (saved to {})".format(total_bin_sd_df.shape[0], os.path.join(work_dir + "raw_SD_binary_map.tsv")))
    
    both_neg_strand = (total_bin_sd_df.loc[:, "strand1"] == "-") & (total_bin_sd_df.loc[:, "strand2"] == "-")
    reverse_strand_dict = {"-": "+", "+": "-"}
    total_bin_sd_df.loc[both_neg_strand, "strand1"] = total_bin_sd_df.loc[both_neg_strand, "strand1"].map(reverse_strand_dict)
    total_bin_sd_df.loc[both_neg_strand, "strand2"] = total_bin_sd_df.loc[both_neg_strand, "strand2"].map(reverse_strand_dict)
                                                           
    ## One XA region (represented as *_bam1) can overlap multiple SDs, keep only the minimal set of SDs with sufficient overlap
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
    total_bin_sd_df.to_csv(os.path.join(work_dir, "filtered_SD_binary_map.tsv"), sep="\t", index=False)
    logger.info(f"After intersecting with refined SD reference, the total targeted SD region size for sample {input_bam} is {total_bin_sd_bed.sort().total_coverage()}bp.")
    logger.debug(f"Final SD map has shape {total_bin_sd_df.shape} and looks like: \n{total_bin_sd_df.head(5).to_string(index=False)}")
    
    # Step 3: Create the SD + PO graph
    graph_path = os.path.join(work_dir, basename + "-multiplexed_homologous_sequences.graphml")
    if os.path.exists(graph_path) and is_file_up_to_date(graph_path, [input_bam]):
        graph = read_graphml(graph_path)
    else:
        graph = create_multiplex_graph( total_bin_sd_df, graph_path, threads = threads)

    # Step 4: Extract SD + PO BED regions with depth issues caused by mapping ambiguity. Mutual overlap exists between intervals.
    ## expanded_total_bin_sd_df is constructed from reference SD map, multi_align_bed_obj is constructed from specified target regions.
    ## Filter out regions with length > 140bp (compared to 150bp reads)
    graph_bed = pb.BedTool.from_dataframe(expanded_total_bin_sd_df)
    graph_bed = filter_bed_by_interval_size(graph_bed, 140)
    
    ## Extract SD + PO BED regions with depth issues caused by mapping ambiguity
    intersect_df = graph_bed.intersect(multi_align_bed_obj.merge(), wo=True).to_dataframe(disable_auto_names=True, names=["chrA", "startA", "endA", "strandA", 
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
        final_query_bed.intersect(multi_align_bed_obj, wo=True).saveas(os.path.join(work_dir, basename + "target_overlapping_query_SD.bed")) # renamed from coding_query_nodes.bed
        logger.debug("Target-overlapping SD regions saved to: {}".format(os.path.join(work_dir, basename + "target_overlapping_query_SD.bed")))
    logger.info("Identified {} distinct target-overlapping SDs from graph".format(query_nodes.shape[0]))

    # Step 5. Pool all SD nodes and the corresponding paralogous regions from the SD + PO graph
    final_graph_path = graph_path.replace(".graphml", ".trim.annoPC.graphml")
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
    #     for node in result["PCs"]:
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
                                  output_folder = work_dir,
                                  ref_genome = ref_genome,
                                  nthreads = threads,
                                  avg_frag_size = avg_frag_size,
                                  std_frag_size = std_frag_size)

def main():
    parser = argparse.ArgumentParser(description='Deploy PCs for SDrecall.')

    parser.add_argument('-r', '--ref_genome', required=True, help='Path to the reference genome.')
    parser.add_argument('-d', '--work_dir', required=True, help='Base directory for output files.')
    parser.add_argument('-i', '--input_bam', required=True, help='Input BAM file.')
    parser.add_argument('-m', '--reference_sd_map', required=True, help='Reference SD map file.')
    parser.add_argument('-b', '--target_bed', default="", help='Optional target BED file.')
    # parser.add_argument('-e', '--error_rate', type=float, default=0.05, help='Error rate for determining reads with extreme template lengths.')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads to use.')
    parser.add_argument('--mq_cutoff', type=int, default=41, help='Mapping quality cutoff.')
    parser.add_argument('--high_quality_depth', type=int, default=10, help='High quality depth cutoff.')
    parser.add_argument('--minimum_depth', type=int, default=3, help='Minimum depth cutoff.')
    parser.add_argument('--multialign_frac', type=float, default=0.5, help='Multi-align fraction cutoff.')
    parser.add_argument('--target_tag', type=str, default="target", help='Optional target tag for filtering.')
    parser.add_argument('-v', '--verbose', type=str, default="INFO", help='Level of verbosity (default = INFO).')

    args = parser.parse_args()

    configure_logger(log_level=args.verbose)

    preparation(
        ref_genome=args.ref_genome,
        work_dir=args.work_dir,
        input_bam=args.input_bam,
        reference_sd_map=args.reference_sd_map,
        target_bed=args.target_bed,
        mq_threshold=args.mq_cutoff,
        high_quality_depth=args.high_quality_depth,
        minimum_depth=args.minimum_depth,
        multialign_frac=args.multialign_frac,
        threads=args.threads,
        target_tag=args.target_tag,
    )

if __name__ == "__main__":
    main()
