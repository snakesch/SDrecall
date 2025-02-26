#!/usr/bin/env python
import os
os.environ["OMP_NUM_THREADS"] = "72"
# os.environ["MKL_NUM_THREADS"] = "20"
# os.environ["OPENBLAS_NUM_THREADS"] = "20"
os.environ["NUMBA_NUM_THREADS"] = "4"
# os.environ["VECLIB_MAXIMUM_THREADS"] = "20"
# os.environ["NUMEXPR_NUM_THREADS"] = "20"
os.environ["TBB_NUM_THREADS"] = "48"


import multiprocessing as mp
ctx = mp.get_context("spawn")

import gc
import io
import subprocess
import pandas as pd
import pybedtools as pb
pb.set_tempdir("/tmp")
import numpy as np
import logging
import re
import uuid
import glob
from scipy import stats
import pysam
from itertools import repeat
import argparse as ap
from python_utils import convert_input_value
import inspect
import time
from datetime import datetime
from pathlib import Path
from io import StringIO
import traceback
import sys
import numba
from numba import set_num_threads, get_num_threads

from annotate_bam_tag_to_vcf import main_annotate as annotate_HP_tag_to_vcf
from misalignment_elimination import main_function as gt_read_clustering_filter
from identify_misaligned_reads_SDrecall import main_process as identify_misaligned_reads
from python_utils import split_bed_by_size


script_path = Path( __file__ ).absolute()
bash_utils_hub = "/paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh"



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# This is a new python script trying to run SDrecall's variants calling workflow for a given sample:
# 1. BAM Reads selection  (Select read pairs from the entire SDs related regions)
# 2. FASTQ Reads selection
# 3. Reads pooling
# 4. Multiploidy short variants calling with HaplotypeCaller
# 5. Compare with intrinsic vcf
# 6. Compare with in-house database


def gt_filter_init(threads):
    os.environ["NUMBA_NUM_THREADS"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["TBB_NUM_THREADS"] = str(threads)
    set_num_threads(4)
    numba.config.THREADING_LAYER = 'omp'
    print(f"Init the subprocess with {get_num_threads()} threads for numba", file=sys.stderr)
    

def executeCmd(cmd, stdout_only = False, shell="/home/yangyxt/miniforge3/envs/ngs_pipeline/bin/bash", logger=logger) -> None:
    from subprocess import PIPE
    logger.info("About to run this command in shell invoked within python: \n{}\n".format(cmd))
    
    if stdout_only:
        result = subprocess.run(cmd, shell=True, executable=shell, stdout=PIPE, stderr=PIPE)
    else:
        result = subprocess.run(cmd, shell=True, executable=shell, stderr=subprocess.STDOUT, stdout=PIPE)
        
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        logger.error("Error in {}\nAnd the output goes like:\n{}\n".format(" ".join(cmd_lst), result.stdout.decode()))
        if cmd_lst[1][0] != "-":
            raise RuntimeError
        else:
            raise RuntimeError
        
    logger.info(f"Ran the following shell command inside python:\n{cmd}\nAnd it receives a return code of {code}, the output goes like this:\n{result.stdout.decode()}\n**********************END_OF_LOG**********************\n\n")
    return result.stdout.decode()


def prepare_tmp_file(tmp_dir="/paedyl01/disk1/yangyxt/test_tmp", **kwargs):
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass
    
    import tempfile
    return tempfile.NamedTemporaryFile(dir = "/paedyl01/disk1/yangyxt/test_tmp", delete = False, **kwargs)


def stat_file_rows(file):
    with open(file, "r") as f:
        lines = f.readlines()
    return len(lines)


def bam_reads_selection_by_region(input_bam: str, 
                                  region_bed: str, 
                                  output_bam = None,
                                  output_freads = None,
                                  output_rreads = None,
                                  input_freads = None,
                                  input_rreads = None,
                                  multi_aligned = False,
                                  threads = 1,
                                  mq_cutoff = 20, 
                                  gatk = f"bash {bash_utils_hub} gatk_wrapper",
                                  output_qnames = "",
                                  logger=logger):
    if os.path.exists(os.path.dirname(output_qnames)):
        tmp_file = output_qnames
    else:
        tmp_file = prepare_tmp_file().name
    
    if multi_aligned:
        cmd = f"""sambamba view -q -t {threads} -L {region_bed} {input_bam} | \
                  mawk -F '\\t' '($0 !~ /SA:Z:/) && (($5 < 60) || ($0 ~ /XA:Z:/)) {{print $1;}}' | tail -n +1 | sort - | uniq - > {tmp_file}"""
    else:
        cmd = f"""sambamba view -q -t {threads} -L {region_bed} {input_bam} | cut -f 1 | tail -n +1 | sort - | uniq - > {tmp_file}"""
    executeCmd(cmd, logger=logger)

    '''
    # The below code block run into some performance issue. Some of the bam files take forever to be sliced.
    if multi_aligned:
        cmd = f"""sambamba slice -q -L {region_bed} {input_bam} | \
                  sambamba view -q -F \"[XA] != null and [SA] == null and mapping_quality <= {mq_cutoff + 10}\" -f sam -t {threads} /dev/stdin | \
                  cut -f 1 | sort - | uniq - > {tmp_file}"""
    else:
        cmd = f"""sambamba slice -q -L {region_bed} {input_bam} | \
                  sambamba view -t {threads} -f sam /dev/stdin | cut -f 1 | sort - | uniq - > {tmp_file}"""
    executeCmd(cmd, logger=logger)
    '''
    
    if stat_file_rows(tmp_file) < 1:
        logger.error(f"The region {region_bed} does not have any poor aligned reads in BAM ${input_bam}")
        if output_qnames:
            return None
        elif output_bam:
            return None
        else:
            return None, None
    
    if output_qnames:
        return tmp_file

    if output_bam:
        cmd = f"""{gatk} FilterSamReads -I ${input_bam} -O ${output_bam} --FILTER includeReadList -RLF ${tmp_file} -SO coordinate"""
        executeCmd(cmd, logger=logger)
        os.remove(tmp_file)
        return output_bam

    if output_freads and output_rreads and input_freads and input_rreads:
        cmdf = f"seqtk subseq {input_freads} {tmp_file} > {output_freads}"
        cmdr = f"seqtk subseq {input_rreads} {tmp_file} > {output_rreads}"
        executeCmd(cmdf, logger=logger)
        executeCmd(cmdr, logger=logger)
        # Make sure the output fastq files are consistent regarding the read qnames
        cmd = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
        executeCmd(cmd, logger=logger)
        # os.remove(tmp_file)
        return output_freads, output_rreads


def extract_reads_by_qnames(qname_lst_f, qname_lst_r,
                            input_freads, input_rreads, 
                            output_freads, output_rreads, 
                            logger = logger):

    cmdf = f"seqtk subseq {input_freads} {qname_lst_f} > {output_freads}"
    cmdr = f"seqtk subseq {input_rreads} {qname_lst_r} > {output_rreads}"

    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)
    # Make sure the output fastq files are consistent regarding the read qnames
    cmds = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
    executeCmd(cmds, logger=logger)

    return output_freads, output_rreads



def imap_process_masked_bam(tup_args):
    return process_masked_bam(*tup_args)

def process_masked_bam( rg_tag,
                        rg_fc_beds,
                        rg_nfc_beds,
                        original_bam,
                        rg_bed,
                        sample_ID,
                        total_freads,
                        total_rreads,
                        max_varno=5,
                        threads=1,
                        mq_cutoff=20,
                        ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                        logger = logger ):

    # First merge the FC and NFC regions and masked genomes
    mq_cutoff = int(mq_cutoff)
    rg_dir = os.path.dirname(rg_bed)
    rg_fc_beds = list(dict.fromkeys(rg_fc_beds))
    rg_nfc_beds = list(dict.fromkeys(rg_nfc_beds))
    genome_name = os.path.basename(ref_genome).replace(".fasta", f".sub.{rg_tag}.masked.fasta")
    masked_genome = os.path.join(rg_dir, genome_name)
    rg_fc_bed = rg_bed.replace(".bed", ".targeted.bed")
    rg_nfc_bed = rg_bed.replace(".bed", ".counterparts_regions.targeted.bed")
    cmd = f"cat {' '.join(rg_fc_beds)} | bedtools sort -i - | bedtools merge -i - > {rg_fc_bed}"
    executeCmd(cmd, logger=logger)
    cmd = f"cat {' '.join(rg_nfc_beds)} | bedtools sort -i - | bedtools merge -i - > {rg_nfc_bed}"
    executeCmd(cmd, logger=logger)

    # Now select the reads from the total fastq files
    per_sample_rg_dir = os.path.join(os.path.dirname(original_bam), f"{sample_ID}_ref_SD_recall")
    os.makedirs(per_sample_rg_dir, exist_ok = True)

    sd_freads = os.path.join(per_sample_rg_dir, os.path.basename(total_freads.replace(".fastq", f".{rg_tag}_sdrecall.fastq")))
    sd_rreads = os.path.join(per_sample_rg_dir, os.path.basename(total_rreads.replace(".fastq", f".{rg_tag}_sdrecall.fastq")))

    assert sd_freads != total_freads, f"The output fastq file {sd_freads} should not be the same as the input fastq file {total_freads}"
    assert sd_rreads != total_rreads, f"The output fastq file {sd_rreads} should not be the same as the input fastq file {total_rreads}"

    masked_bam = os.path.join(per_sample_rg_dir, os.path.basename(original_bam.replace(".bam", f".only_{rg_tag}.bam")))
    raw_masked_bam = masked_bam.replace(".bam", ".raw.bam")
    raw_masked_vcf = masked_bam.replace(".bam", ".raw.vcf.gz")

    # Extract the qnames
    rg_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                 rg_fc_bed, 
                                                 output_qnames = "Yes",
                                                 threads=threads,
                                                 multi_aligned=False,
                                                 logger=logger )

    counterpart_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                          rg_nfc_bed, 
                                                          output_qnames = "Yes",
                                                          threads=threads,
                                                          multi_aligned=True,
                                                          mq_cutoff=mq_cutoff,
                                                          logger=logger )


    # merged_qname_lst = prepare_tmp_file().name
    qname_lst_f = prepare_tmp_file().name
    qname_lst_r = prepare_tmp_file().name
    cmdf = f"cat {rg_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_f}"
    cmdr = f"cat {rg_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_r}"
    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)

    sd_freads, sd_rreads = extract_reads_by_qnames( qname_lst_f, qname_lst_r,
                                                    total_freads, total_rreads, 
                                                    sd_freads, sd_rreads, 
                                                    logger=logger )

    # Now perform the mapping
    cmd = f"bash {bash_utils_hub} independent_minimap2_masked \
            -b {original_bam} \
            -a {masked_genome} \
            -s {sample_ID} \
            -f {sd_freads} \
            -r {sd_rreads} \
            -g {ref_genome} \
            -o {raw_masked_bam} \
            -t {threads} \
            -i {rg_tag} \
            -c {max_varno} && ls -lh {raw_masked_bam}"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        return f"{rg_tag},{rg_bed},NaN,NaN"
    else:
        # Now perform the variants calling
        # First determine the names of the output VCFs
        logger.info(f"Now we have prepared the masked BAM file {raw_masked_bam} for sample {sample_ID} for {rg_bed}, we need to generate the VCF file.")

        sub_running_log = raw_masked_vcf.replace(".vcf.gz", ".log")
        cmd = f"bash {bash_utils_hub} call_polyploidy_per_PC \
                -a {original_bam} \
                -m {ref_genome} \
                -c {threads} \
                -o {raw_masked_vcf} \
                -p {rg_tag} \
                -b {raw_masked_bam} > {sub_running_log} 2>&1"
            # We need to allow this function to be failed sometimes
        try:
            executeCmd(cmd, logger=logger)
        except RuntimeError:
            logger.error(f"Failed to generate {raw_masked_vcf}, running log in {sub_running_log}")
            cmd = f"bash {bash_utils_hub} check_vcf_validity {raw_masked_vcf} 1 && [[ {raw_masked_vcf} -nt {script_path} ]]"
            try:
                executeCmd(cmd, logger=logger)
            except RuntimeError:
                return f"{rg_tag},{rg_bed},{raw_masked_bam},NaN"
            else:
                logger.warning(f"Though the process reported error. The VCF file {raw_masked_vcf} is still valid and updated")
                return f"{rg_tag},{rg_bed},{raw_masked_bam},{raw_masked_vcf}"
        else:
            logger.info(f"Succesfully generate the vcf file {raw_masked_vcf} for region ${rg_bed}")
            return f"{rg_tag},{rg_bed},{raw_masked_bam},{raw_masked_vcf}"



def merge_bams(bam_list: list, 
               merged_bam: str, 
               ref_fasta = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
               threads = 2,
               logger = logger):
    # For data visualization and debugging
    merged_bam_header = merged_bam.replace(".bam", ".header")
    cmd = f"bash {bash_utils_hub} modify_bam_sq_lines {bam_list[0]} {ref_fasta} {merged_bam_header}"
    executeCmd(cmd, logger=logger)
    
    test_cmd = f"bash {bash_utils_hub} quick_check_bam_validity {merged_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        logger.info(f"The merged pooled BAM file {merged_bam} is not valid")
        execute = True
    else:
        if all([os.path.getmtime(merged_bam) > os.path.getmtime(vb) for vb in bam_list]):
            execute = False
        else:
            logger.info(f"The merged pooled BAM file {merged_bam} is valid but not updated")
            execute = True

    merged_bam_list = merged_bam.replace(".bam", ".bams.list.txt")
    with open(merged_bam_list, "w") as f:
        f.write("\n".join(bam_list))
    
    cmd = f"samtools merge -c -@ {threads} -h {merged_bam_header} -b {merged_bam_list} -o - | \
            samtools sort -O bam -o {merged_bam} -@ {threads} && \
            samtools index {merged_bam} && \
            ls -lht {merged_bam} || \
            echo Failed to concatenate all the pooled BAM files. It wont be a fatal error but brings troubles to debugging and variant tracing."
    if execute:
        executeCmd(cmd, logger = logger)

    return execute
        
        
def merge_vcfs(vcf_list, merged_vcf, threads=4, check_error = True, chunk_size=300, logger = logger):
    # Then we need to use bcftools concat to merge the vcfs together
    tmp_list = prepare_tmp_file().name
    if check_error:
        ce_arg = ""
    else:
        ce_arg = "-e yes"

    if len(vcf_list) < chunk_size:
        with open(tmp_list, "w") as f:
            for v in vcf_list:
                f.write(v + "\n")

        cmd = f"bash {bash_utils_hub} bcftools_concatvcfs \
                -v {tmp_list} \
                -o {merged_vcf} \
                -t {threads} \
                {ce_arg} && \
                bash {bash_utils_hub} display_table {merged_vcf} 20"
        executeCmd(cmd, logger = logger)
        return
    
    i = 0
    middle_vcfs = []
    while i < len(vcf_list):
        use_list = vcf_list[i:min(i+chunk_size, len(vcf_list))]
        with open(tmp_list, "w") as f:
            for v in use_list: 
                f.write(v + "\n")
        
        middle_vcf = merged_vcf.replace(".vcf.gz", f".{i}.vcf.gz")
        # The input vcf validity is already done checking in the below bash function
        cmd = f"bash {bash_utils_hub} bcftools_concatvcfs \
                -v {tmp_list} \
                -o {middle_vcf} \
                -t {threads} \
                {ce_arg} && \
                bash {bash_utils_hub} display_table {middle_vcf} 20"

        executeCmd(cmd, logger = logger)
        middle_vcfs.append(middle_vcf)
        i += chunk_size

    # Now merge the middle vcfs
    with open(tmp_list, "w") as f:
        for v in middle_vcfs:
            f.write(v + "\n")

    cmd = f"bash {bash_utils_hub} bcftools_concatvcfs \
            -v {tmp_list} \
            -o {merged_vcf} \
            -t {threads} \
            {ce_arg} && \
            bash {bash_utils_hub} display_table {merged_vcf} 20"
    executeCmd(cmd)

def imap_slice_bam_per_bed(tup_args):
    return slice_bam_per_bed(*tup_args)

def slice_bam_per_bed(bed, bam, ref_genome, threads = 4, logger = logger):

    chunk_id = os.path.basename(bed).split(".")[-2]
    cov_bam = bam.replace(".bam", f".{chunk_id}.bam")
    cov_bam_header = cov_bam.replace(".bam", ".header")
    if not os.path.exists(bam + ".bai") or os.path.getmtime(bam) > os.path.getmtime(bam + ".bai"):
        tmp_index = prepare_tmp_file(suffix=".bai").name
        cmd = f"samtools index -b -o {tmp_index} {bam} && \
                mv {tmp_index} {bam + '.bai'}"
        executeCmd(cmd, logger = logger)
        # We need to make sure the index update is atomic

    cmd = f"sambamba slice -q -L {bed} {bam} | \
            sambamba sort -q -t {threads} -o {cov_bam} /dev/stdin && \
            bash {bash_utils_hub} modify_bam_sq_lines {cov_bam} {ref_genome} {cov_bam_header} && \
            samtools reheader {cov_bam_header} {cov_bam} > {cov_bam}.tmp && \
            mv {cov_bam}.tmp {cov_bam} && \
            samtools index {cov_bam}"
    try:
        executeCmd(cmd, logger = logger)
    except RuntimeError:
        logger.warning(f"Failed to slice the bam file {bam} by bed {bed} and generate a {cov_bam}")
        return "NaN"
    else:
        return cov_bam
    


def split_bam_by_cov(bam, 
                     chunksize=10000, 
                     delimiter_size=1000,
                     base_qual = 0,
                     map_qual = 0,
                     beds = [],
                     threads = 4,
                     ref_genome = "",
                     logger = logger):
    
    if len(beds) == 0:
        cov_bed = bam.replace(".bam", f".bed")
        cmd = f"bash {bash_utils_hub} samtools_bam_coverage \
                -i {bam} \
                -d {delimiter_size} \
                -o {cov_bed} \
                -b {base_qual} \
                -m {map_qual}"
        executeCmd(cmd, logger = logger)

        beds = split_bed_by_size(cov_bed, 
                                chunksize = chunksize, 
                                logger = logger)

        logger.info(f"Getting the coverage bed file {cov_bed} and split it into {len(beds)} parts")
    else:
        logger.info(f"Using the provided coverage bed files({len(beds)} beds) to split the BAM file {bam}")

    with ctx.Pool(threads - 1) as pool:
        result_records = pool.imap(imap_slice_bam_per_bed, zip( beds, 
                                                                repeat(bam),
                                                                repeat(ref_genome), 
                                                                repeat(1)))

        splitted_bams = []
        i=0
        for success, result, log_contents in result_records:
            i+=1
            print(f"\n************************************{i}_subprocess_start_for_slice_bam************************************\n", file=sys.stderr)
            if success:
                splitted_bams.append(result)
                print(f"Successfully slice the bam file {bam}. The log info are:\n{log_contents}\n", file=sys.stderr)
            else:
                error_mes, tb_str = result
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")

            print(f"\n************************************{i}_subprocess_end_for_slice_bam************************************\n", file=sys.stderr)
     
    logger.info(f"Successfully splitted the BAM file {bam} into {len(splitted_bams)} parts")
    return splitted_bams, beds



def post_process_bam(raw_bam: str,
                     raw_vcf: str,
                     output_vcf: str,
                     original_bam: str,
                     intrinsic_bam: str,
                     threads = 12,
                     numba_threads = 4,
                     conf_level = 0.01,
                     stat_sample_size = 1000000,
                     histo_fig = None,
                     varno_cutoff = 3,
                     max_varno = 5,
                     ref_genome = "",
                     avg_frag_size = 500,
                     logger = logger):
    '''
    Input raw bam is the sample-level merged raw bam
    Now we need to split it by chromosomes and parallel run filter_out_noisy_reads to save time and increase computation efficiency 
    '''
    
    remapped_bam = original_bam.replace(".bam", ".pooled.bam")
    remapped_raw_bam = original_bam.replace(".bam", ".pooled.raw.bam")
    remapped_clean_bam = original_bam.replace(".bam", ".pooled.clean.bam")
    remapped_noise_bam = original_bam.replace(".bam", ".pooled.noise.bam")
    clean_vcf = output_vcf.replace(".vcf.gz", ".clean.vcf.gz")
    
    splitted_bams, splitted_beds = split_bam_by_cov(raw_bam, 
                                                    delimiter_size=int(np.ceil(2*avg_frag_size)), 
                                                    logger = logger, 
                                                    threads = threads,
                                                    ref_genome = ref_genome)
    splitted_intrin_bams, splitted_intrin_beds = split_bam_by_cov(intrinsic_bam, 
                                                                  beds = splitted_beds, 
                                                                  logger = logger, 
                                                                  threads = threads,
                                                                  ref_genome = ref_genome)

    # There is a possiblity that the splitted_intrin_bams might be empty, so we need to check it
    # Sort the list for load balancing
    tup_bams = list(zip(splitted_bams, splitted_intrin_bams, splitted_beds))
    assert len(splitted_intrin_bams) == len(splitted_bams), f"The splitted intrin BAMs number is not equal to splitted BAMs number"
    tb_removed_indices = set()
    for i in range(len(splitted_bams)):
        if splitted_bams[i] == "NaN":
            logger.warning(f"The BAM file {splitted_bams[i]} for {splitted_beds[i]} is not valid, so we need to remove it from the list")
            tb_removed_indices.add(i)
        elif splitted_intrin_bams[i] == "NaN":
            logger.warning(f"The intrinsic BAM file {splitted_intrin_bams[i]} for {splitted_beds[i]} is not valid, so we need to remove it from the list")
            tb_removed_indices.add(i)

    tup_bams = [tup_bams[i] for i in range(len(tup_bams)) if i not in tb_removed_indices]
    tup_bams = sorted(tup_bams, key=lambda x: os.path.getsize(x[0]), reverse=True)
    raw_bams = [cb[0] for cb in tup_bams]
    intrinsic_bams = [cb[1] for cb in tup_bams]
    raw_bam_regions = [cb[2] for cb in tup_bams]
    masked_bams = [ re.sub(r"\.raw\.", ".", cb) for cb in raw_bams]
    noisy_bams = [ re.sub(r"\.raw\.", ".noise.", cb) for cb in raw_bams]

    # Filter out the misaligned_reads per chromosome
    numba.set_num_threads(numba_threads)
    logger.info(f"Now start the parallel filtering on the raw realigned BAM files. There are in total {len(raw_bams)} BAM files to be processed. Each subprocess has up to {get_num_threads()} threads to use.")
    parallelism = int(np.ceil(threads-numba_threads))

    with ctx.Pool(parallelism, initializer=gt_filter_init, initargs=(numba_threads,)) as pool:
        result_records = pool.imap_unordered(imap_filter_out, zip(raw_bams,
                                                                  repeat(raw_vcf),
                                                                  raw_bam_regions,
                                                                  repeat(original_bam),
                                                                  masked_bams,
                                                                  noisy_bams,
                                                                  intrinsic_bams,
                                                                  repeat(conf_level),
                                                                  repeat(stat_sample_size),
                                                                  repeat(histo_fig),
                                                                  repeat(1),
                                                                  repeat(numba_threads),
                                                                  repeat(varno_cutoff),
                                                                  repeat(max_varno),
                                                                  repeat(ref_genome)))

        result_record_strs = []
        i=0
        for result in result_records:
            i+=1
            print(f"\n************************************{i}_subprocess_start_for_filtering************************************\n", file=sys.stderr)
            result_record_strs.append(result)
            assert len(result.split(",")) == 5, f"The returned result {result} should only have 5 fields separated by comma"     
            raw_bam = result.split(",")[0]       
            print(f"{datetime.now()}: ************************************{i}_subprocess_end_for_filtering_{raw_bam}************************************", file=sys.stderr)

    gc.collect()
    logger.info("Start to compose the parallel returned results into a dataframe")
    post_result_tab = os.path.join(os.path.dirname(raw_bam), "post_result.tsv")
    post_result_df = pd.read_csv(io.StringIO("\n".join(result_record_strs)), header=None, names=["raw_masked_bam", 
                                                                                                "masked_bam", 
                                                                                                "clean_masked_bam", 
                                                                                                "noise_masked_bam", 
                                                                                                "masked_vcf"], na_values="NaN")
    post_result_df.to_csv(post_result_tab, sep="\t", index=False)
    logger.info(f"The post processed masked bam and vcf files (stored in {post_result_tab}) across different chromosomes looks like:\n{post_result_df.to_string(index=False)}\n\n")

    # Now we merge the masked bams and vcfs across different chromosomes

    merge_bams(post_result_df.loc[:, "masked_bam"].dropna().drop_duplicates().tolist(), remapped_bam, ref_fasta = ref_genome, threads=threads - 1, logger=logger)
    merge_bams(post_result_df.loc[:, "clean_masked_bam"].dropna().drop_duplicates().tolist(), remapped_clean_bam, ref_fasta = ref_genome, threads=threads - 1, logger=logger)
    merge_bams(post_result_df.loc[:, "noise_masked_bam"].dropna().drop_duplicates().tolist(), remapped_noise_bam, ref_fasta = ref_genome, threads=threads - 1, logger=logger)
    merge_bams(post_result_df.loc[:, "raw_masked_bam"].dropna().drop_duplicates().tolist(), remapped_raw_bam, ref_fasta = ref_genome, threads=threads - 1, logger=logger)

    # Merge the vcfs
    merge_vcfs( post_result_df.loc[:, "masked_vcf"].dropna().drop_duplicates().tolist(), 
                clean_vcf, 
                threads = threads - 1, 
                logger = logger,
                check_error=False )

    
def filter_out_reads_per_region(raw_bam: str,
                                raw_vcf: str,
                                bam_region: str,
                                original_bam: str,
                                masked_bam: str,
                                noisy_bam: str,
                                intrinsic_bam: str,
                                conf_level = 0.01,
                                sample_size = 1000000,
                                histo_fig = None,
                                threads = 4,
                                numba_threads = 4,
                                varno_cutoff = 3,
                                max_varno = 5,
                                ref_genome = "",
                                logger = logger):

    logger.info(f"Input max_varno is {max_varno} (take the larger one from a pair with 5), varno_cutoff is {varno_cutoff}. ")
    max_varno = max(5, max_varno)
    
    # Extract the file path from the logger
    filter_logger_file = re.sub(r"\.bam", ".filter.log", raw_bam)
    executeCmd(f": > {filter_logger_file}", logger = logger)
    logger = file_logger(filter_logger_file)
    sl_stderr = StreamToLogger(logger, logging.WARNING)
    sys.stderr = sl_stderr
    set_num_threads(numba_threads)
    print(f"{datetime.now()}: Start the filtering on raw_bam {raw_bam}, the numba threads are set to {get_num_threads()}. log file path: {filter_logger_file}", file = sys.stderr)
    
    ref_gen_tag = os.path.basename(ref_genome).split(".")[1]
    chrom = os.path.basename(raw_bam).split(".")[-2]
    assert masked_bam != original_bam, f"Masked bam file {masked_bam} should not be the same as the original bam file {original_bam}"
    # Here raw_bam is .raw.bam files

    phasing_graph = gt_read_clustering_filter(raw_bam,
                                            output_bam=masked_bam,
                                            filter_out_bam=noisy_bam,
                                            intrinsic_bam=intrinsic_bam,
                                            raw_vcf=raw_vcf,
                                            bam_region_bed=bam_region,
                                            varno_cutoff=varno_cutoff,
                                            max_varno=max_varno,
                                            raw_intrinsic_bam=f"/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/hg19/test_BISER/assembly_intrinsic_align.WGAC.{ref_gen_tag}.reformat.bam",
                                            logger = logger)

    if phasing_graph == None:
        logger.warning(f"No variants detected in {raw_bam} for region {bam_region}. Skip the filtering process.")
        close_logger(logger)
        return f"{raw_bam},NaN,NaN,NaN,NaN"

    cmd = f"ls -lht {masked_bam} {noisy_bam}"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.warning(f"Failed to generate filtered BAM file {masked_bam} for region {chrom}")
        close_logger(logger)
        return f"{raw_bam},NaN,NaN,NaN,NaN"
    
    noisy_vars = noisy_bam.replace(".bam", ".vcf.gz")
    clean_bam = masked_bam.replace(".bam", ".clean.bam")
    
    # Call variants from the noisy bam and these variants will be need to generate a new clean bam file
    cmd = f'''bcftools mpileup --threads {threads} --indels-2.0 -a FORMAT/AD,FORMAT/DP -q 10 -Q 15 -f {ref_genome} {noisy_bam} | \
    bcftools call --threads {threads} -mv -P 0 -Ou -f GQ | \
    bcftools norm --threads {threads} -m -both -f {ref_genome} --multi-overlaps 0 -a -Ou - | \
    bcftools norm --threads {threads} -d exact - | \
    bcftools view --threads {threads} -i 'ALT!="*"' -Ov | \
    bcftools filter -s "LIKELY_INTRINSIC" -e 'ALT != "*"' -Ov - | \
    mawk 'BEGIN{{FS=OFS="\\t";}} {{printf "%s\\t%s\\t%s\\t%s\\t%s", $1, $2, $3, toupper($4), toupper($5); \
                            for(i=6;i<=NF;i++) {{printf "\\t%s", $i; }} \
                            printf "\\n"; }}' - | \
    bcftools sort -Oz -o {noisy_vars} && \
    tabix -f -p vcf {noisy_vars} && \
    ls -lht {noisy_vars}'''
    
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.warning(f"Failed to generate noisy VCF file {noisy_vars} for region {chrom}")
        close_logger(logger)
        return f"{raw_bam},{masked_bam},NaN,{noisy_bam},NaN"
    
    # Use the noisy variants bam to identify the misaligned reads and filter out them to generate an even more clean bam
    cleaned_bam = identify_misaligned_reads(noisy_vars, 
                                            masked_bam, 
                                            output_bam = clean_bam, 
                                            occurence_cutoff = max_varno,
                                            threads = threads,
                                            logger = logger)
    assert clean_bam == cleaned_bam, f"Returned file name {cleaned_bam} is different with the preset clean bam name {clean_bam}"
    cleaned_bam = clean_bam.replace(".bam", ".tmp.bam")
    cmd = f"samtools sort -O bam -o {cleaned_bam} -@ {threads} {clean_bam} && \
            mv {cleaned_bam} {clean_bam} && \
            samtools index {clean_bam}"
            
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.warning(f"Failed to generate cleaned BAM file {clean_bam} for region {chrom}")
        close_logger(logger)
        return f"{raw_bam},{masked_bam},NaN,{noisy_bam},NaN"
    
    # Now perform the variant calling on two bam files
    # First determine the names of the output VCFs
    logger.info(f"Now we have prepared the masked BAM file {masked_bam} for {chrom}, we need to generate the VCF file.")
    masked_vcf = masked_bam.replace(".bam", ".vcf.gz")
    
    # Only .raw.vcf and .vcf two versions of vcf. 
    cmd = f"bash {bash_utils_hub} call_polyploidy_per_PC \
            -a {original_bam} \
            -m {ref_genome} \
            -c {threads} \
            -o {masked_vcf} \
            -b {masked_bam} \
            -p \".\" && ls -lht {masked_vcf}"
    
    try:
        executeCmd(cmd, logger=logger)
        # Annotate the masked_vcf with HP tag from the masked_bam file
        tmp_vcf = masked_vcf.replace(".vcf.gz", ".tmp.vcf.gz")
        tmp_vcf = annotate_HP_tag_to_vcf(masked_vcf, masked_bam, tmp_vcf, "HP", logger = logger)
        cmd = f"mv {tmp_vcf} {masked_vcf} && bcftools index -f -t {masked_vcf}"
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.warning(f"Failed to generate {masked_vcf} for region {chrom}")
        close_logger(logger)
        return f"{raw_bam},{masked_bam},{clean_bam},{noisy_bam},NaN"
    else:
        logger.info(f"Succesfully generate the vcf file {masked_vcf} for region {chrom}")      
        close_logger(logger)
        return f"{raw_bam},{masked_bam},{clean_bam},{noisy_bam},{masked_vcf}"
    
    

def imap_filter_out(args):
    try:
        return filter_out_reads_per_region(*args)
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}\nTraceback: {traceback.format_exc()}")
        logger.error(f"The input arguments are {args}")
        raise e
    

