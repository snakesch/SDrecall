## Realign --> realignment + variant calling
import os
import io
import sys
import multiprocessing as mp
ctx = mp.get_context("spawn")

import gc
import pysam
import pandas as pd
import pybedtools as pb

from typing import List, Set, Tuple, Dict
from itertools import repeat
from multiprocessing import Pool

from src.utils import executeCmd, prepare_tmp_file
from src.const import shell_utils, SDrecallPaths
from src.log import logger

from realign_recall.stat_realign_group_regions import stat_all_RG_region_size
from realign_recall.realign_per_RG import imap_process_masked_bam
from realign_recall.cal_edge_NM_values import calculate_NM_distribution_poisson


def pool_init():
    global imap_prepare_masked_align_region_per_RG
    from realign_recall.prepare_masked_align_region import imap_prepare_masked_align_region_per_RG


def SDrecall_per_sample(sdrecall_paths: SDrecallPaths,
                        threads = 12,
                        numba_threads = 4,
                        mq_cutoff = 20,
                        conf_level = 0.01,
                        varno_cutoff = 3,
                        histo_fig = None):
    # First calculate input bam fragment size distribution
    input_bam = sdrecall_paths.input_bam
    avg_frag_size, std_frag_size = sdrecall_paths.avg_frag_size, sdrecall_paths.frag_size_std
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")

    prepared_arguments_df = stat_all_RG_region_size(sdrecall_paths, threads)

    fc_size_stat_tab = os.path.join(sdrecall_paths.tmp_dir, "FC_size_stat.tsv")
    logger.info(f"The prepared arguments for the PC subgroups will be saved to {fc_size_stat_tab} looks like:\n{prepared_arguments_df[:20].to_string(index=False)}\n\n")
    prepared_arguments_df.to_csv(fc_size_stat_tab, sep="\t", index=False)

    uniq_rg_labels = prepared_arguments_df.loc[:, "rg_label"].drop_duplicates().tolist()
    rg_subids_tup_list = [tuple(prepared_arguments_df.loc[prepared_arguments_df["rg_label"] == rg, "subgroup_id"].drop_duplicates().tolist()) for rg in uniq_rg_labels]
    all_region_beds = [sdrecall_paths.all_homo_regions_bed_path(rg) for rg in uniq_rg_labels]
    
    ref_genome = sdrecall_paths.ref_genome
    target_region = sdrecall_paths.target_bed
    pool = ctx.Pool(threads, initializer=pool_init)
    prepared_beds = pool.imap_unordered(imap_prepare_masked_align_region_per_RG, zip(uniq_rg_labels,
                                                                                    rg_subids_tup_list,
                                                                                    repeat(sdrecall_paths.tmp_dir),
                                                                                    repeat(target_region),
                                                                                    all_region_beds,
                                                                                    repeat(ref_genome)))
    
    result_beds = []
    i=0
    for success, result, log_contents in prepared_beds:
        i+=1
        print(f"\n************************************{i}_subprocess_start_for_prepare_masked_align_beds************************************\n", file=sys.stderr)
        if success:
            for record in result:
                result_beds.append(record)
            print(f"Successfully prepared the masked bam file. The log info are:\n{log_contents}\n", file=sys.stderr)
        else:
            error_mes, tb_str = result
            result_beds.append(tb_str)
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")
        print(f"\n************************************{i}_subprocess_end_for_prepare_masked_align_beds************************************\n", file=sys.stderr)

    pool.close()
    gc.collect()

    prepared_beds_df = pd.read_csv(io.StringIO("\n".join(result_beds)), header=None, names=["rg_label", "rg_subgroup_id", "fc_bed", "nfc_bed", "fc_bed_size", "nfc_bed_size"], na_values="NaN")
    logger.info(f"The prepared beds for the masked alignment looks like:\n{prepared_beds_df[:20].to_string(index=False)}\n\n")
    
    # Prepare the input arguments for parallel execution of realignment and recall
    uniq_rgs = prepared_beds_df.loc[:, "rg_label"].drop_duplicates().tolist()
    uniq_rgs = sorted(uniq_rgs, key = lambda x: prepared_beds_df.loc[prepared_beds_df["rg_label"] == x, "nfc_bed_size"].sum(), reverse=True)  # Load balancing
    rg_query_beds = [ sdrecall_paths.rg_query_bed_path(rg) for rg in uniq_rgs ]
    rg_fc_beds = [ prepared_beds_df.loc[prepared_beds_df["rg_label"] == rg, "fc_bed"].tolist() for rg in uniq_rgs ]
    rg_nfc_beds = [ prepared_beds_df.loc[prepared_beds_df["rg_label"] == rg, "nfc_bed"].tolist() for rg in uniq_rgs]
    rg_raw_masked_bams = [ sdrecall_paths.rg_raw_masked_bam_path(rg) for rg in uniq_rgs ]
    rg_masked_genomes = [ sdrecall_paths.masked_genome_path(rg) for rg in uniq_rgs ]
    fastq_tuple_list = [sdrecall_paths.rg_realign_fastqs_path(rg_label) for rg_label in uniq_rgs]
    sd_freads, sd_rreads = zip(*fastq_tuple_list)

    # Perform the realignment and recall
    pool = ctx.Pool(threads)
    results = pool.imap_unordered(imap_process_masked_bam, zip( uniq_rgs,
                                                                rg_query_beds,
                                                                rg_fc_beds,
                                                                rg_nfc_beds,
                                                                repeat(input_bam),
                                                                rg_raw_masked_bams,
                                                                rg_masked_genomes,
                                                                repeat(sdrecall_paths.sample_id),
                                                                sd_freads,
                                                                sd_rreads,
                                                                repeat(1),
                                                                repeat(mq_cutoff),
                                                                repeat(ref_genome) ))

    result_records = []
    i=0
    for success, result, log_contents in results:
        i+=1
        print(f"\n************************************{i}_subprocess_start_for_process_masked_align************************************\n", file=sys.stderr)
        if success:
            result_records.append(result)
            print(f"Successfully filtered the masked bam file and call variants on it: {result}. The log info are:\n{log_contents}\n", file=sys.stderr)
        else:
            error_mes, tb_str = result
            result_records.append(tb_str)
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")
        print(f"\n************************************{i}_subprocess_end_for_process_masked_align************************************\n", file=sys.stderr)

    pool.close()
    gc.collect()

    realign_meta_table = sdrecall_paths.realign_meta_table_path()
    post_result_df = pd.read_csv(io.StringIO("\n".join(result_records)), header=None, names=["rg_label", "rg_subgroup_id", "raw_masked_bam", "raw_masked_vcf"], na_values="NaN")
    post_result_df.to_csv(realign_meta_table, sep="\t", index=False)
    logger.info(f"The processed masked bam and vcf files stored in {realign_meta_table} and it looks like:\n{post_result_df.to_string(index=False)}\n\n")

