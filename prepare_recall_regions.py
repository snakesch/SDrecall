import argparse
import cProfile
import logging
import os
from multiprocessing import Pool
import sys

import pandas as pd
import pybedtools as pb
import networkx as nx
import graph_tool.all

from build_regions_from_PC_clusters import establish_beds_per_PC_cluster
from seq import get_bam_frag_size
from extract_SD_pairs_from_bam_json import extract_sd_coordinates_from_json
from pick_multialign_regions import pick_region_by_depth
from genome import Genome
from homoseq_region import HOMOSEQ_REGION
from graph_build import read_graphml, create_multiplex_graph
from graph_query import extract_SD_paralog_pairs_from_graph, query_connected_nodes
from convert_nodes_into_hierachical_beds import convert_nodes_into_hierachical_beds
from sd_pairs import filter_umbrella_pairs

from log import ColoredFormatter
from suppress_warning import *
from utils import is_file_up_to_date, executeCmd, filter_bed_by_interval_size

logger = logging.getLogger('SDrecall')

def deploy_PCs_for_SDrecall_main(ref_genome: str,
                                 work_dir: str, 
                                 input_bam: str,
                                 reference_sd_map: str,
                                 target_bed = "",
                                 err_rate = 0.05, 
                                 threads = 10,
                                 mq_cutoff = 20,
                                 target_tag = "target",
                                 profile_file = None):

    

    # if profile_file:
    #     pr = cProfile.Profile()
    #     pr.enable()
    
    os.makedirs(work_dir, exist_ok=True)

    rg = Genome(path=ref_genome)
    genome_file = rg.fai_index

    # Step 0: Calculate the distribution of fragment sizes
    # avg_frag_size, std_frag_size = get_bam_frag_size(input_bam)
    avg_frag_size, std_frag_size = 570.4, 150.7
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")
    
    ## Define names of intermediate files. Intermediate / final outputs will be written to work_dir.
    basename = os.path.basename(input_bam)
    qname_lst = os.path.join(work_dir, basename).replace(".bam", f".{target_tag}.qname.lst")
    bam_json = os.path.join(work_dir, basename).replace(".bam", f".{target_tag}.json")   
    multi_align_bed = os.path.join(work_dir, basename).replace(".bam", f".{target_tag}.multialign.bed")
    
    # Step 1 : Pick the multi_aligned regions within the target regions
    if not os.path.exists(multi_align_bed) or not is_file_up_to_date(multi_align_bed, [input_bam, target_bed]):
        multi_align_bed = pick_region_by_depth(input_bam, 
                                        multi_align_bed, 
                                        MQ_threshold=41, 
                                        high_quality_depth=10, 
                                        minimum_depth=3, 
                                        multialign_frac=0.5,
                                        target_region=target_bed, 
                                        target_tag="FCRs",
                                        genome_file=genome_file)

    multi_align_bed_obj = pb.BedTool(multi_align_bed).sort()
    logger.info("The multialign bed extracted from {} is {} and it covers {} bp.".format(input_bam,
                                                                                      multi_align_bed,
                                                                                      multi_align_bed_obj.total_coverage()))
    
    # Step 2: Extract from the BAM file multialigned regions identified in step 1.
    ## Then extract the XA read_pairs with low MAPQ and output the results as JSON by sambamba
    cmd = f"sambamba slice -q -L {multi_align_bed} {input_bam} | \
            sambamba view -q -F \"[SA] == null and mapping_quality <= {mq_cutoff + 10} and [XA] != null\" -f sam -t {threads} /dev/stdin | \
            cut -f 1 | sort - | uniq - > {qname_lst} && \
            samtools view -N {qname_lst} -@ {threads} -u {input_bam} | \
            sambamba view -q -f json -o {bam_json} -t {threads} /dev/stdin"
            
    if os.path.exists(bam_json):
        if is_file_up_to_date(bam_json, [input_bam, target_bed]) and os.path.getsize(bam_json) > 28:
            logger.info(f"{bam_json} already exists and updated")
        else:
            executeCmd(cmd)
    else:
        executeCmd(cmd)
    
    # Step 3: Extract paired SD coordinates from JSON alignments
    sd_map_str = extract_sd_coordinates_from_json(bam_json, 
                                             avg_frag_size=avg_frag_size, 
                                             std_frag_size=std_frag_size, 
                                             err_rate=err_rate,
                                             threads=threads)

    ## Exclude non-primary contigs
    regex_pattern = r'(chr)*(X|Y|MT*|[0-9][0-9]*)$'
    sd_map_df = pd.DataFrame([ l.split(",") for l in sd_map_str.split("\n") ], columns = ["chr", "start", "end", "chr_counterparts", "start_counterparts", "end_counterparts"]).drop_duplicates()
    all_involved_chrs = set(sd_map_df.loc[:, "chr"])
    all_involved_chrs.update(set(sd_map_df.loc[:, "chr_counterparts"]))
    logger.debug(f"Prior to excluding non-primary contig SDs, SD map contains {sd_map_df.shape[0]} SDs from {all_involved_chrs}")
    sd_map_df = sd_map_df.loc[sd_map_df.loc[:, "chr"].str.match(regex_pattern) & sd_map_df.loc[:, "chr_counterparts"].str.match(regex_pattern), :]
    all_involved_chrs = set(sd_map_df.loc[:, "chr"])
    all_involved_chrs.update(set(sd_map_df.loc[:, "chr_counterparts"]))
    logger.debug(f"After excluding non-primary contig SDs, SD map contains {sd_map_df.shape[0]} SDs from {all_involved_chrs}")

    if len(all_involved_chrs) > 25:
        logger.error("Unable to filter alternative contigs.")
        sys.exit(1)

    ## Next identify the target SDs of interest from the primary SDs (first three columns of sd_map_df) by intersecting multi_aiign-bed_obj
    sd_map_bed = pb.BedTool.from_dataframe(sd_map_df).intersect(multi_align_bed_obj,
                                                             wa=True)
    sd_map_df = sd_map_bed.to_dataframe(disable_auto_names=True, names = ["chr", "start", "end", 
                                                                          "chr_counterparts", "start_counterparts", "end_counterparts"]).drop_duplicates()
    
    ## Expand SD map to include paralogous regions in subsequent variant recall
    reverse_sd_df = sd_map_df.loc[:, ["chr_counterparts", "start_counterparts", "end_counterparts", "chr", "start", "end"]].rename(columns={"chr_counterparts":"chr", "start_counterparts":"start", "end_counterparts":"end", "chr":"chr_counterparts", "start":"start_counterparts", "end":"end_counterparts"})
    merged_sd_df = pd.concat([sd_map_df, reverse_sd_df], axis=0, ignore_index=True).loc[:,["chr", "start", "end"]].drop_duplicates()
    sd_map_bed = pb.BedTool.from_dataframe(merged_sd_df)
    bam_sd_bed = os.path.join(work_dir, basename).replace(".bam", f".{target_tag}.bamsd.bed")
    sd_map_bed.saveas(bam_sd_bed)
    logger.info(f"SD regions from primary contigs are written to {bam_sd_bed} (covering {sd_map_bed.sort().total_coverage()}bp). ")
    logger.debug(f"Primary contig SD map: {sd_map_df.head(10).to_string(index=False)}")
    
    # Step 4: Compare paired SD regions with the reference SD map to identify SD regions associated with mapping ambiguity
    ref_bed_obj = pb.BedTool(reference_sd_map)
    logger.info("Total coverage of reference SD map: {}bp".format(ref_bed_obj.sort().total_coverage()))

    ## Filter out the pairs where paired intervals are smaller than the avg fragment size.
    big_ref_bed_obj = filter_bed_by_interval_size(ref_bed_obj, avg_frag_size)
    logger.info(f"After excluding paired reference SDs with fragment size < {avg_frag_size}bp, the reference SD map covers {big_ref_bed_obj.total_coverage()}bp. ")
    
    ## Intersect with previously derived sample-specific SD map
    total_bin_sd_bed = big_ref_bed_obj.intersect(sd_map_bed, wo=True)
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
    
    # Step 5: Create the SD + PO graph
    graph_path = os.path.join(work_dir, basename + "-multiplexed_homologous_sequences.graphml")
    if os.path.exists(graph_path) and is_file_up_to_date(graph_path, [input_bam, target_bed, bam_json]):
        graph = read_graphml(graph_path)
    else:
        graph = create_multiplex_graph( total_bin_sd_df, graph_path, threads = threads, target_bed = target_bed)

    # Step 6: Extract SD + PO BED regions with depth issues caused by mapping ambiguity. Mutual overlap exists between intervals.
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
        final_query_bed.intersect(multi_align_bed_obj, wo=True).saveas(os.path.join(work_dir, basename).replace(".bam", "target_overlapping_query_SD.bed")) # renamed from coding_query_nodes.bed
        logger.debug("Target-overlapping SD regions saved to: {}".format(os.path.join(work_dir, basename).replace(".bam", "target_overlapping_query_SD.bed")))
    logger.info("Identified {} distinct target-overlapping SDs from graph".format(query_nodes.shape[0]))

    # Step 7. Pool all SD nodes and the corresponding paralogous regions from the SD + PO graph
    final_graph_path = graph_path.replace(".graphml", ".trim.annoPC.graphml")
    sd_paralog_pairs, connected_qnode_components = extract_SD_paralog_pairs_from_graph( query_nodes, 
                                                                                graph, 
                                                                                graph_path = final_graph_path, 
                                                                                avg_frag_size = avg_frag_size, 
                                                                                std_frag_size = std_frag_size, 
                                                                                threads=threads )
    logger.info("{} SD-paralog pairs pooled from SD + PO graph".format(len(sd_paralog_pairs)))
    
    # Step 8: Collapse qnodes by identifying qnodes that can be put in the same masked genome
    grouped_qnode_cnodes = query_connected_nodes(sd_paralog_pairs, connected_qnode_components)

    # Step 8.5: Enumerate qnodes and cnodes in the graph. (FC -> qnode; NFC -> cnode)
    final_graph = graph.copy()
    for i, result in enumerate(grouped_qnode_cnodes):
        tag = "PC" + str(i)
        for node in result["PCs"]:
            if not isinstance(node, tuple):
                raise TypeError(f"Expected {node} (type: {type(node)}) to be tuple.")
            final_graph.nodes[node]["FC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("FC", "").split(",") + [tag])) if len(e) > 0])
        for node in result["SD_counterparts"]:
            if not isinstance(node, HOMOSEQ_REGION):
                raise TypeError(f"Expected {node} (type: {type(node)}) to be HOMOSEQ_REGION. ")
            final_graph.nodes[node.data]["NFC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node.data].get("NFC", "").split(",") + [tag])) if len(e) > 0])
    nx.write_graphml(final_graph, final_graph_path)

    # Step 9: Create beds and masked genomes
    convert_nodes_into_hierachical_beds(grouped_qnode_cnodes = grouped_qnode_cnodes,
                                        sd_paralog_pairs = sd_paralog_pairs,
                                        output_folder = work_dir,
                                        ref_genome = ref_genome,
                                        target_region_bed=target_bed,
                                        nthreads = threads,
                                        avg_frag_size = avg_frag_size,
                                        std_frag_size = std_frag_size)

    # if profile_file:
    #     pr.disable()
    #     pr.dump_stats(profile_file)
    
    
def main():
    parser = argparse.ArgumentParser(description='Deploy PCs for SDrecall.')

    parser.add_argument('-r', '--ref_genome', required=True, help='Path to the reference genome.')
    parser.add_argument('-d', '--work_dir', required=True, help='Base directory for output files.')
    parser.add_argument('-i', '--input_bam', required=True, help='Input BAM file.')
    parser.add_argument('-m', '--reference_sd_map', required=True, help='Reference structural variant map file.')
    parser.add_argument('-b', '--target_bed', default="", help='Optional target BED file.')
    parser.add_argument('-e', '--error_rate', type=float, default=0.05, help='Error rate for determining reads with extreme template lengths.')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads to use.')
    parser.add_argument('--mq_cutoff', type=int, default=20, help='Mapping quality cutoff.')
    parser.add_argument('--target_tag', type=str, default="target", help='Optional target tag for filtering.')
    parser.add_argument('-p', '--profile_file', type=str, default=None, help='Optional profile file for tuning.')
    parser.add_argument('-v', '--verbose', type=str, default="INFO", help='Level of verbosity (default = INFO).')

    args = parser.parse_args()

    logger.setLevel(getattr(logging, args.verbose))
    console_handler=logging.StreamHandler()
    console_handler.setLevel(getattr(logging, args.verbose))
    formatter = ColoredFormatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Call the main function with parsed arguments
    deploy_PCs_for_SDrecall_main(
        ref_genome=args.ref_genome,
        work_dir=args.work_dir,
        input_bam=args.input_bam,
        reference_sd_map=args.reference_sd_map,
        target_bed=args.target_bed,
        err_rate=args.error_rate,
        threads=args.threads,
        mq_cutoff=args.mq_cutoff,
        target_tag=args.target_tag,
        profile_file=args.profile_file
    )

if __name__ == "__main__":
    main()
