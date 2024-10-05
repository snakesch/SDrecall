import argparse
import cProfile
import logging
import os
# import re
from multiprocessing import Pool
import sys

import pandas as pd
import pybedtools as pb
import networkx as nx
import graph_tool.all

from graph_traversal import traverse_network_to_get_homology_counterparts
from build_regions_from_PC_clusters import establish_beds_per_PC_cluster
from seq import get_bam_frag_size
from extract_SD_pairs_from_bam_json import extract_sd_coordinates_from_json
from pick_multialign_regions import pick_region_by_depth
from genome import Genome
from homoseq_region import HOMOSEQ_REGION
from graph import read_graphml, create_multiplex_graph
from graph_query import extract_FC_NFC_pairs_from_graph, query_connected_nodes
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
                                 fraction_cutoff = 0.7,
                                 err_rate = 0.05, 
                                 threads = 10,
                                 aggregation_resolution = .5,
                                 mq_cutoff = 20,
                                 target_tag = "target",
                                 profile_file = None):

    def imap_establish(tup_args):
        return establish_beds_per_PC_cluster(*tup_args)

    def imap_traverse(tup_args):
        return traverse_network_to_get_homology_counterparts(*tup_args)

    # if profile_file:
    #     pr = cProfile.Profile()
    #     pr.enable()
    
    os.makedirs(work_dir, exist_ok=True)

    rg = Genome(path=ref_genome)
    genome_file = rg.fai_index

    # Step 0: Calculate the distribution of fragment sizes
    # avg_frag_size, std_frag_size = get_bam_frag_size(input_bam)
    avg_frag_size, std_frag_size = 571.1, 120
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}")
    
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

    multi_align_bed_obj = pb.BedTool(multi_align_bed)
    logger.info("The multialign bed extracted from {} is {} and it covers {} bp.".format(input_bam,
                                                                                      multi_align_bed,
                                                                                      multi_align_bed_obj.sort().total_coverage()))
    
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
    
    # Remove the duplicated SD combinations (i.e. same SD pair in different order)
    total_bin_sd_df.loc[:, "frozenset_indx"] = total_bin_sd_df.apply(lambda row: frozenset({ row["chr_1"]+":"+str(row["start_1"])+"-"+str(row["end_1"])+":"+row["strand1"], row["chr_2"]+":"+str(row["start_2"])+"-"+str(row["end_2"])+":"+row["strand2"]}), axis=1)
    total_bin_sd_df = total_bin_sd_df.drop_duplicates(subset="frozenset_indx").drop(columns=["frozenset_indx"])
    total_bin_sd_df.to_csv(os.path.join(work_dir, "filtered_SD_binary_map.tsv"), sep="\t", index=False)
    logger.info(f"After intersecting with refined SD reference, the total targeted SD region size for sample {input_bam} is {total_bin_sd_bed.sort().total_coverage()}bp.")
    logger.debug(f"Final SD map has shape {total_bin_sd_df.shape} and looks like: \n{total_bin_sd_df.head(5).to_string(index=False)}")
    
    # Step 5: Create a graph on the total_bin_sd_df
    graph_path = os.path.join(work_dir, basename + "-multiplexed_homologous_sequences.graphml")
    if os.path.exists(graph_path) and is_file_up_to_date(graph_path, [input_bam, target_bed, bam_json]):
        graph = read_graphml(graph_path)
    else:
        graph = create_multiplex_graph( total_bin_sd_df, graph_path, 
                                        overlap_frac_min = fraction_cutoff, 
                                        avg_frag_size = avg_frag_size, 
                                        std_frag_size = std_frag_size, 
                                        resolution = aggregation_resolution,
                                        threads = threads, 
                                        target_bed = target_bed)
    
    # Step 6: Pick out the nodes that overlaps with target region, these intervals come from filtered reference SD map. Intervals might tangle(overlap) with each other
    graph_bed = pb.BedTool.from_dataframe(expanded_total_bin_sd_df.drop_duplicates())
    # Filter out the graph bed where the interval length should be larger than 150 bp (common read length)
    graph_bed = filter_bed_by_interval_size(graph_bed, 140)
    
    # This is to identify which graph nodes overlap with the target regions to pool the reads to
    intersect_df = graph_bed.intersect(multi_align_bed_obj.sort().merge(), wo=True).to_dataframe(disable_auto_names=True, names=["chrA", "startA", "endA", "strandA", 
                                                                                                                                 "chrB", "startB", "endB", "strandB",
                                                                                                                                 "mismatch_rate",
                                                                                                                                 "chr_query", "start_query", "end_query", "overlap_len"])
    intersect_df = intersect_df.loc[(intersect_df["endA"] > intersect_df["startA"]) & \
                                    (intersect_df["end_query"] > intersect_df["start_query"]) & \
                                    (intersect_df["endB"] > intersect_df["startB"]), :]

    by_query_interval = intersect_df.groupby(["chr_query", "start_query", "end_query"], as_index=False)
    query_groups = [group for _, group in by_query_interval]
    
    # We just need to select the smallest FC_enwrapped HS intervals 
    # Update: This filtering might filter out SD regions that map to other genomic regions so temporarily decaprecate it and observe the computation performance downstream
    logger.info(f"Before the filtering, there are {intersect_df.shape[0]} intervals overlaps with the target regions: \n{intersect_df[:20].to_string(index=False)}\n")
    with Pool(threads) as pool:
        results = pool.imap_unordered(filter_umbrella_pairs, query_groups)
        intersect_df = pd.concat(results, axis=0, ignore_index=True).drop_duplicates()
    
    logger.info("Only {} necessary intervals overlaps with the target region intervals: \n{}\n".format(intersect_df.shape[0], intersect_df[:20].to_string(index=False)))
    
    query_nodes = intersect_df.loc[:, ["chrA", "startA", "endA", "strandA"]].drop_duplicates().rename(columns={"chrA":"chr", "startA":"start", "endA":"end", "strandA":"strand"})
    final_query_bed = pb.BedTool.from_dataframe(query_nodes)
    final_query_bed.intersect(pb.BedTool(multi_align_bed), wo=True).saveas(os.path.join(work_dir, "coding_query_nodes.bed"))
    logger.info("Save the coding overlap query nodes to this file: {}\n".format(os.path.join(work_dir, "coding_query_nodes.bed")))
    
    logger.info("Finally after picking up the merged the SD intervals, there are {} merged SDs from graph overlapped with exon, and it covers a region of {}bp.".format(query_nodes.shape[0],
                                                                                                                                                                        final_query_bed.sort().total_coverage()))
    
    # Step 7. Now use a function to generate FC and NFC pairs for all the query nodes
    final_graph_path = graph_path.replace(".graphml", ".trim.annoPC.graphml")
    fc_nfc_pairs, connected_qnode_components = extract_FC_NFC_pairs_from_graph( query_nodes, 
                                                                                graph, 
                                                                                graph_path = final_graph_path, 
                                                                                avg_frag_size = avg_frag_size, 
                                                                                std_frag_size = std_frag_size, 
                                                                                threads=threads )
    logger.info("After extracting the FC-NFC pairs, we got {} FC-NFC pairs".format(len(fc_nfc_pairs)))
    
    # Now we need to test which FC-NFC pairs can be merged together (meaning FCs cannot appear in NFCs covered region)
    # Step 8: Pass the query nodes to function and get the cluster of SD intervals (The query nodes here is a pandas dataframe object)
    pair_graph_path = graph_path.replace(".graphml", ".pair.graphml")
    logger.info("Now we about to use the extract FC-NFC pairs to find out which of them can be put in the same masked genome")
    all_results, fc_nfc_dict = query_connected_nodes(fc_nfc_pairs, 
                                                     connected_qnode_components, 
                                                     multiplex_graph = graph,
                                                     threads = threads,
                                                     avg_frag_size = avg_frag_size,
                                                     std_frag_size = std_frag_size,
                                                     graph_path = pair_graph_path)

    # Step 8.5: We need to tag the nodes with FC-NFC pair ID
    final_graph = graph.copy()
    for i, result in enumerate(all_results):
        tag = "PC" + str(i)
        for node in result["PCs"]:
            # Add an attribute "FC" to the node
            assert type(node) == tuple, f"The counterpart node {node} type is not tuple, the FC nodes are {result['PCs']}"
            final_graph.nodes[node]["FC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("FC", "").split(",") + [tag])) if len(e) > 0])
        for node in result["SD_counterparts"]:
            # Add an attribute "NFC" to the node
            assert type(node) == HOMOSEQ_REGION, f"The counterpart node {node} type is not HOMOSEQ_REGION, the FC nodes are {result['PCs']}"
            final_graph.nodes[node.data]["NFC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node.data].get("NFC", "").split(",") + [tag])) if len(e) > 0])
    nx.write_graphml(final_graph, final_graph_path)

    # Step 9: Create beds and masked genomes
    convert_nodes_into_hierachical_beds(all_results = all_results,
                                        fc_nfc_dict = fc_nfc_dict,
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
    parser.add_argument('--fraction_cutoff', type=float, default=0.7, help='Fraction cutoff for analysis.')
    parser.add_argument('-e', '--error_rate', type=float, default=0.05, help='Error rate for determining reads with extreme template lengths.')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads to use.')
    parser.add_argument('--aggregation_resolution', type=float, default=0.5, help='Resolution for aggregation.')
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
        fraction_cutoff=args.fraction_cutoff,
        err_rate=args.error_rate,
        threads=args.threads,
        aggregation_resolution=args.aggregation_resolution,
        mq_cutoff=args.mq_cutoff,
        target_tag=args.target_tag,
        profile_file=args.profile_file
    )

if __name__ == "__main__":
    main()
