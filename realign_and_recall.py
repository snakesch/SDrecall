## Realign --> realignment + variant calling
import os
import io
import sys
import multiprocessing as mp
ctx = mp.get_context("spawn")

import gc
import pandas as pd

from itertools import repeat

from src.utils import executeCmd, merge_bams, configure_parallelism
from src.const import SDrecallPaths, shell_utils
from src.log import logger
from src.merge_variants_with_priority import merge_with_priority

from realign_recall.stat_realign_group_regions import stat_all_RG_region_size
from realign_recall.realign_per_RG import imap_process_masked_bam
from realign_recall.annotate_HP_tag_to_vars import annotate_vcf as annotate_vcf_HP_tag
from realign_recall.prepare_masked_align_region import imap_prepare_masked_align_region_per_RG
from misalignment_elimination import eliminate_misalignments


def pool_init():
    global imap_prepare_masked_align_region_per_RG
    from realign_recall.prepare_masked_align_region import imap_prepare_masked_align_region_per_RG


def SDrecall_per_sample(sdrecall_paths: SDrecallPaths,
                        threads = 12,
                        numba_threads = 4,
                        mq_cutoff = 20,
                        conf_level = 0.01):
    # First calculate input bam fragment size distribution
    print("\n"*2, "*"*100, file=sys.stderr)
    logger.info(f"START PERFORMING REALIGNMENT AND RECALL FOR {sdrecall_paths.sample_id} based on the REALIGNMENT GROUP REGIONS in {sdrecall_paths.work_dir}")
    print("*"*100, "\n"*2, file=sys.stderr)
    input_bam = sdrecall_paths.input_bam
    avg_frag_size, std_frag_size = sdrecall_paths.avg_frag_size, sdrecall_paths.frag_size_std
    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}bp (std: {std_frag_size}bp)")

    prepared_arguments_df = stat_all_RG_region_size(sdrecall_paths)

    fc_size_stat_tab = os.path.join(sdrecall_paths.tmp_dir, "FC_size_stat.tsv")
    logger.info(f"The prepared arguments for the PC subgroups will be saved to {fc_size_stat_tab} looks like:\n{prepared_arguments_df[:10].to_string(index=False)}\n\n")
    prepared_arguments_df.to_csv(fc_size_stat_tab, sep="\t", index=False)

    uniq_rg_labels = prepared_arguments_df.loc[:, "rg_label"].drop_duplicates().tolist()
    rg_subids_tup_list = [tuple(prepared_arguments_df.loc[prepared_arguments_df["rg_label"] == rg, "subgroup_id"].drop_duplicates().tolist()) for rg in uniq_rg_labels]
    all_region_beds = [sdrecall_paths.all_homo_regions_bed_path(rg) for rg in uniq_rg_labels]
    
    ref_genome = sdrecall_paths.ref_genome
    target_region = sdrecall_paths.target_bed

    # Prepare the regions to recruit reads for downstream realignment
    num_jobs, _ = configure_parallelism(threads, 1)
    pool = ctx.Pool(num_jobs, initializer=pool_init)
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
    num_jobs, threads_per_job = configure_parallelism(threads, 4)
    pool = ctx.Pool(num_jobs, initializer=pool_init)
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
                                                                repeat(threads_per_job),
                                                                repeat(mq_cutoff),
                                                                repeat(ref_genome) ))

    result_records = []
    i=0
    for success, result, log_contents in results:
        i+=1
        print(f"\n************************************{i}_subprocess_start_for_process_masked_align************************************\n", file=sys.stderr)
        if success:
            result_records.append(result)
            print(f"Successfully filtered the masked bam file: {result}. The log info are:\n{log_contents}\n", file=sys.stderr)
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
    logger.info(f"The processed masked bam files stored in {realign_meta_table} and it looks like:\n{post_result_df.to_string(index=False)}\n\n")

    # Merge the pooled raw bam files
    pooled_raw_bam = sdrecall_paths.pooled_raw_bam_path()
    execute_merging = merge_bams(bam_list = post_result_df["raw_masked_bam"].dropna().unique().tolist(), 
                                 merged_bam = pooled_raw_bam, 
                                 ref_fasta = ref_genome,
                                 tmp_dir = sdrecall_paths.tmp_dir,
                                 threads = threads,
                                 logger = logger)
    deduped_raw_bam = pooled_raw_bam.replace(".bam", ".deduped.bam")
    cmd = f"samtools collate -@ {threads} -O -u {pooled_raw_bam} | \
            samtools fixmate -@ {threads} -m -u - - | \
            samtools sort -T {sdrecall_paths.tmp_dir} -@ {threads} -u - | \
            samtools markdup -@ {threads} -r - {deduped_raw_bam} && \
            ls -lh {deduped_raw_bam}"
    if execute_merging:
        executeCmd(cmd, logger=logger)

    # Merge variants to have the raw vcf file directly from unfiltered realignments
    pooled_raw_vcf = sdrecall_paths.recall_raw_vcf_path()
    vcf_list = post_result_df["raw_masked_vcf"].dropna().unique().tolist()

    if len(vcf_list) > 0:
        # Create a temporary file with the list of VCFs to merge
        vcf_list_file = os.path.join(sdrecall_paths.tmp_dir, "vcfs_to_merge.txt")
        with open(vcf_list_file, 'w') as f:
            for vcf in vcf_list:
                f.write(f"{vcf}\n")
        
        logger.info(f"Merging {len(vcf_list)} VCF files into {pooled_raw_vcf}")
        cmd = f"bash {shell_utils} bcftools_concatvcfs \
                -v {vcf_list_file} \
                -o {pooled_raw_vcf} \
                -c {threads} \
                -t {sdrecall_paths.tmp_dir} \
                -s {sdrecall_paths.sample_id}"
        executeCmd(cmd, logger=logger)
    else:
        logger.warning("No valid VCF files found to merge.")

    # - Now we start filtering the pooled raw bam file by phasing and mislalignment elimination - #
    total_intrinsic_bam = sdrecall_paths.total_intrinsic_bam_path()
    pooled_filtered_bam = sdrecall_paths.pooled_filtered_bam_path()
    print("\n"*2, "*"*100, file=sys.stderr)
    logger.info(f"START PERFORMING MISALIGNMENT ELIMINATION FOR {sdrecall_paths.sample_id} on file {deduped_raw_bam}")
    print("*"*100, "\n"*2, file=sys.stderr)
    deduped_raw_bam, pooled_filtered_bam = eliminate_misalignments(deduped_raw_bam,
                                                                  pooled_filtered_bam,
                                                                  total_intrinsic_bam,
                                                                  ref_genome,
                                                                  avg_frag_size=avg_frag_size,
                                                                  threads=threads,
                                                                  numba_threads=numba_threads,
                                                                  conf_level=conf_level,
                                                                  mapq_cutoff=mq_cutoff,
                                                                  cache_dir=sdrecall_paths.tmp_dir,
                                                                  logger=logger)

    # Call variants on the filtered BAM file after misalignment elimination
    pooled_filtered_vcf = sdrecall_paths.recall_filtered_vcf_path()
    logger.info(f"Calling variants on the filtered BAM file {pooled_filtered_bam}")
    cmd = f"bash {shell_utils} bcftools_call_per_RG \
            -m {ref_genome} \
            -b {pooled_filtered_bam} \
            -o {pooled_filtered_vcf} \
            -c {threads} \
            -p {sdrecall_paths.sample_id}"
    executeCmd(cmd, logger=logger)

    # Now we need to merge the pooled filtered vcf and the pooled raw vcf to identify which variants might be derived from misalignments
    pooled_filtered_vcf = annotate_vcf_HP_tag(pooled_filtered_vcf, 
                                              pooled_filtered_vcf, 
                                              pooled_filtered_bam, 
                                              "HP", 
                                              logger = logger)
    pooled_raw_vcf = annotate_vcf_HP_tag(pooled_raw_vcf, 
                                         pooled_raw_vcf, 
                                         deduped_raw_bam, 
                                         "HP", 
                                         logger = logger)
    
    final_recall_vcf = sdrecall_paths.final_recall_vcf_path()
    merge_with_priority(query_vcf = pooled_raw_vcf, 
                        reference_vcf = pooled_filtered_vcf, 
                        output_vcf = final_recall_vcf, 
                        added_filter = "MISALIGNED", 
                        qv_tag = "RAW", 
                        rv_tag = "CLEAN", 
                        ref_genome = ref_genome, 
                        threads = threads)

    print("\n"*2, "*"*100, file=sys.stderr)
    logger.info(f"FINISHED PERFORMING REALIGNMENT AND RECALL FOR {sdrecall_paths.sample_id} based on the REALIGNMENT GROUP REGIONS in {sdrecall_paths.work_dir}, the final recall vcf is {final_recall_vcf}")
    print("*"*100, "\n"*2, file=sys.stderr)

    return final_recall_vcf

