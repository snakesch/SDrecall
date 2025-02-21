#!/usr/bin/env python
import os
os.environ["OMP_NUM_THREADS"] = "4"
# os.environ["MKL_NUM_THREADS"] = "20"
# os.environ["OPENBLAS_NUM_THREADS"] = "20"
os.environ["NUMBA_NUM_THREADS"] = "4"
# os.environ["VECLIB_MAXIMUM_THREADS"] = "20"
# os.environ["NUMEXPR_NUM_THREADS"] = "20"
os.environ["TBB_NUM_THREADS"] = "4"


import multiprocessing as mp
ctx = mp.get_context("spawn")

import gc
import io
import subprocess
import pandas as pd
import pybedtools as pb
pb.set_tempdir("/paedyl01/disk1/yangyxt/test_tmp")
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
import cProfile
import numba
from subprocess import PIPE
from numba import set_num_threads, get_num_threads


from compare_vcf_recs_pysam import main as compare_vcf_recs
from annotate_bam_tag_to_vcf import main_annotate as annotate_HP_tag_to_vcf
from gt_read_clustering_filter import main_function as gt_read_clustering_filter
from identify_misaligned_reads_SDrecall import main_process as identify_misaligned_reads
from python_utils import split_bed_by_size
from prepare_masked_align_region import imap_prepare_masked_align_region_per_PC


script_path = Path( __file__ ).absolute()
bash_utils_hub = "/paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh"


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# This workflow is supposed to be run at the CPOS server using the main conda environment
# This is a new python script trying to run SDrecall's variants calling workflow for a given sample:
# 1. BAM Reads selection  (Select read pairs from the entire SDs related regions)
# 2. FASTQ Reads selection
# 3. Reads pooling
# 4. Multiploidy short variants calling with HaplotypeCaller
# 5. Compare with intrinsic vcf
# 6. Compare with in-house database


def init_logger(handler = logging.StreamHandler(), tag = ""):
    logger = logging.getLogger(f"Process-{tag}")
    handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def close_logger(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)


def file_logger(log_file):
    logger = logging.getLogger(f"Process-{uuid.uuid4()}")
    handler = logging.FileHandler(log_file)
    handler.setFormatter(logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass


def log_command(func):
    def wrapper(*args, **kwargs):
        arg_names = list(inspect.signature(func).parameters.keys())
        arg_values = list(args)
        arg_dict = dict(zip(arg_names, arg_values))
        arg_dict.update(kwargs)

        tmp_tag = str(uuid.uuid4())
        log_stream = StringIO()
        ch = logging.StreamHandler(log_stream)
        logger = init_logger(handler = ch, tag = tmp_tag)

        # Get the default values of the function's keyword arguments
        sig = inspect.signature(func)
        defaults = {k: v.default for k, v in sig.parameters.items() if v.default is not inspect.Parameter.empty}

        # Fill in missing keyword arguments with their default values
        for arg_name, default_value in defaults.items():
            if arg_name not in arg_dict:
                arg_dict[arg_name] = default_value

        # Convert arguments and their values to a string
        args_str = ', '.join([f"{k}={arg_dict[k]}" for k in arg_names])
        defaults_str = ', '.join([f"{k}={arg_dict[k]}" for k in defaults])

        logger.info(f"Executing: {func.__name__}({args_str}, {defaults_str})")
        start = time.time()
        try:
            result = func(*args, logger=logger, **kwargs)
            end = time.time()
            logger.info(f"Finished: {func.__name__}({args_str}, {defaults_str}) in {end - start} seconds")
            log_contents = log_stream.getvalue()
            log_stream.close()
        except Exception as e:
            end = time.time()
            logger.info(f"Failed: {func.__name__}({args_str}, {defaults_str}) in {end - start} seconds")
            tb_str = traceback.format_exc()
            log_contents = log_stream.getvalue()
            log_stream.close()
            return (False, (e, tb_str), log_contents)
        else:
            return (True, result, log_contents)

    return wrapper


def pool_init():
    global imap_prepare_masked_align_region_per_PC
    from prepare_masked_align_region import imap_prepare_masked_align_region_per_PC
    
    
def gt_filter_init(threads):
    os.environ["NUMBA_NUM_THREADS"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["TBB_NUM_THREADS"] = str(threads)
    set_num_threads(threads)
    numba.config.THREADING_LAYER = 'omp'
    print(f"Init the subprocess with {get_num_threads()} threads for numba", file=sys.stderr)
    


def executeCmd(cmd, stdout_only = False, logger = logger) -> None:
    logger.info(f"##### About to run the following shell command inside python:\n{cmd}\n")
    
    if stdout_only:
        result = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    else:
        result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=PIPE)
        
    logger.info(f"##### Just ran the following shell command inside python:\n{cmd}\nAnd the output goes like this:\n{result.stdout.decode()}\n\n")
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        if cmd_lst[1][0] != "-":
            raise RuntimeError(f"Error in {cmd}:\n{result.stdout.decode()}") if not stdout_only else RuntimeError(f"Error in {cmd}:\n{result.stderr.decode()}")
        else:
            raise RuntimeError(f"Error in {cmd}:\n{result.stdout.decode()}") if not stdout_only else RuntimeError(f"Error in {cmd}:\n{result.stderr.decode()}")
    
    return result.stdout.decode()


# def executeCmd(cmd, stdout_only = False, shell="/home/yangyxt/miniforge3/envs/ngs_pipeline/bin/bash", logger=logger) -> None:
#     from subprocess import PIPE
#     logger.info("About to run this command in shell invoked within python: \n{}\n".format(cmd))
    
#     if stdout_only:
#         result = subprocess.run(cmd, shell=True, executable=shell, stdout=PIPE, stderr=PIPE)
#     else:
#         result = subprocess.run(cmd, shell=True, executable=shell, stderr=subprocess.STDOUT, stdout=PIPE)
        
#     code = result.returncode
#     cmd_lst = cmd.split(" ")
#     if code != 0:
#         logger.error("Error in {}\nAnd the output goes like:\n{}\n".format(" ".join(cmd_lst), result.stdout.decode()))
#         if cmd_lst[1][0] != "-":
#             raise RuntimeError
#         else:
#             raise RuntimeError
        
#     logger.info(f"Ran the following shell command inside python:\n{cmd}\nAnd it receives a return code of {code}, the output goes like this:\n{result.stdout.decode()}\n**********************END_OF_LOG**********************\n\n")
#     return result.stdout.decode()


def update_plain_file_on_md5(old_file, new_file):
    import hashlib
    import uuid
    
    if os.path.exists(old_file):
        old_md5 = hashlib.md5(open(old_file,'r').read().encode()).hexdigest()
    else:
        old_md5 = str(uuid.uuid4())

    if not os.path.exists(new_file):
        raise FileNotFoundError("The new file {} does not exist so no updates should be carried out.".format(new_file))
    
    if os.stat(new_file).st_size == 0:
        raise FileExistsError("The new file {} input is completely empty. Quit using it to update the original file {}".format(new_file, old_file))

    new_md5 = hashlib.md5(open(new_file,'r').read().encode()).hexdigest()
    
    if new_md5 == old_md5:
        logger.warning("The new file {} shares the identical content with the old one {} so no updates should be carried out. And the new file {} should be deleted".format(new_file, 
                                                                                                                                                                                old_file,
                                                                                                                                                                                new_file))
        os.remove(new_file)
        return False
    else:
        import shutil
        logger.info("The new file {} is different from the old one {} so the old one will be replaced with the new one.".format(new_file, old_file))
        shutil.move(new_file, old_file)
        executeCmd(f"ls -lht {old_file}")
        return True

    
def update_gzip_file_on_md5(old_file, new_file):
    # use gzip python package to open the file content
    import gzip
    import hashlib
    import uuid
    
    if os.path.exists(old_file):
        old_md5 = hashlib.md5(gzip.open(old_file,'rb').read()).hexdigest()
    else:
        old_md5 = str(uuid.uuid4())

    if not os.path.exists(new_file):
        raise FileNotFoundError("The new file {} does not exist so no updates should be carried out.".format(new_file))
    
    if os.stat(new_file).st_size == 0:
        raise FileExistsError("The new file {} input is completely empty. Quit using it to update the original file {}".format(new_file, old_file))

    new_md5 = hashlib.md5(gzip.open(new_file,'rb').read()).hexdigest()
    
    if new_md5 == old_md5:
        logger.warning("The new file {} shares the identical content with the old one {} so no updates should be carried out. And the new file {} should be deleted".format(new_file, 
                                                                                                                                                                                old_file,
                                                                                                                                                                                new_file))
        os.remove(new_file)
        return False
    else:
        import shutil
        logger.info("The new file {} is different from the old one {} so the old one will be replaced with the new one.".format(new_file, old_file))
        shutil.move(new_file, old_file)
        executeCmd(f"ls -lht {old_file}")
        return True


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
                                  output_freads = None,
                                  output_rreads = None,
                                  input_freads = None,
                                  input_rreads = None,
                                  multi_aligned = False,
                                  threads = 1,
                                  mq_cutoff = 20, 
                                  output_qnames = "",
                                  logger=logger):
    if os.path.exists(os.path.dirname(output_qnames)):
        tmp_file = output_qnames
    else:
        tmp_file = prepare_tmp_file().name
    
    if multi_aligned:
        cmd = f"""sambamba view -q -t {threads} -L {region_bed} {input_bam} | \
                  mawk -F '\\t' '($0 !~ /SA:Z:/) && (($5 < 40) || ($0 ~ /XA:Z:/)) {{print $1;}}' | tail -n +1 | sort - | uniq - > {tmp_file}"""
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
        else:
            return None, None
    
    if output_qnames:
        return tmp_file

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
    # cmda = f"awk '{{printf \"%s/2\n\", $0;}}' {qname_lst} > temp.txt"
    cmdf = f"seqtk subseq {input_freads} {qname_lst_f} > {output_freads}"
    cmdr = f"seqtk subseq {input_rreads} {qname_lst_r} > {output_rreads}"
    # add "/1" and "/2" in the end
    # cmdf = f"seqtk subseq {input_freads} {qname_lst_f} | sed -e 's+/1++g' > {output_freads}"
    # cmdr = f"seqtk subseq {input_rreads} {qname_lst_r} | sed -e 's+/2++g' > {output_rreads}"
    # executeCmd(cmda, logger=logger)
    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)
    # Make sure the output fastq files are consistent regarding the read qnames
    cmds = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
    executeCmd(cmds, logger=logger)
    # os.remove(tmp_file)
    return output_freads, output_rreads





def split_PC_bed(pc_bed: str):
    pc_bed_dir = os.path.dirname(pc_bed)
    pc_bed_name = os.path.basename(pc_bed)
    pc_bed_tag = str(pc_bed_name.replace(".bed", ""))
    
    sub_beds_dir = os.path.join(pc_bed_dir, pc_bed_name.replace(".bed", ".sub_beds"))
    try:
        os.mkdir(sub_beds_dir)
    except FileExistsError:
        pass
    
    result_str = executeCmd(f"bedtools sort -i {pc_bed} | bedtools merge -i -", stdout_only=True)
    lines = [l for l in result_str.split("\n") if len(l) > 3]
    chromosomes = set(line.split('\t')[0] for line in lines)
    # Then split the PC bed one line at a time
    sub_bed_paths = []
    for chromosome in chromosomes:
        sub_bed_path = os.path.join(sub_beds_dir, f"{pc_bed_tag}_{chromosome}.bed")
        with open(sub_bed_path, "w") as sub_bed_file:
            for line in lines:
                if line.split('\t')[0] == chromosome:
                    sub_bed_file.write(line + "\n")
        sub_bed_paths.append(sub_bed_path)
        
    # Remove the historical beds
    existed_beds = next(os.walk(sub_beds_dir))[-1]
    historical_beds = [ os.path.join(sub_beds_dir, b) for b in existed_beds if str(os.path.join(sub_beds_dir, b)) not in sub_bed_paths ]
    for hb in historical_beds: os.remove(hb)
        
    return sub_bed_paths


def imap_process_masked_bam(tup_args):
    return process_masked_bam(*tup_args)


@log_command
def process_masked_bam( pc_tag,
                        pc_fc_beds,
                        pc_nfc_beds,
                        original_bam,
                        pc_bed,
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
    pc_dir = os.path.dirname(pc_bed)
    pc_fc_beds = list(dict.fromkeys(pc_fc_beds))
    pc_nfc_beds = list(dict.fromkeys(pc_nfc_beds))
    genome_name = os.path.basename(ref_genome).replace(".fasta", f".sub.{pc_tag}.masked.fasta")
    masked_genome = os.path.join(pc_dir, genome_name)
    pc_fc_bed = pc_bed.replace(".bed", ".targeted.bed")
    pc_nfc_bed = pc_bed.replace(".bed", ".counterparts_regions.targeted.bed")
    cmd = f"cat {' '.join(pc_fc_beds)} | bedtools sort -i - | bedtools merge -i - > {pc_fc_bed}"
    executeCmd(cmd, logger=logger)
    cmd = f"cat {' '.join(pc_nfc_beds)} | bedtools sort -i - | bedtools merge -i - > {pc_nfc_bed}"
    executeCmd(cmd, logger=logger)

    # Now select the reads from the total fastq files
    per_sample_pc_dir = os.path.join(os.path.dirname(original_bam), f"{sample_ID}_ref_SD_recall")
    os.makedirs(per_sample_pc_dir, exist_ok = True)

    sd_freads = os.path.join(per_sample_pc_dir, os.path.basename(total_freads.replace(".fastq", f".{pc_tag}_sdrecall.fastq")))
    sd_rreads = os.path.join(per_sample_pc_dir, os.path.basename(total_rreads.replace(".fastq", f".{pc_tag}_sdrecall.fastq")))

    assert sd_freads != total_freads, f"The output fastq file {sd_freads} should not be the same as the input fastq file {total_freads}"
    assert sd_rreads != total_rreads, f"The output fastq file {sd_rreads} should not be the same as the input fastq file {total_rreads}"

    masked_bam = os.path.join(per_sample_pc_dir, os.path.basename(original_bam.replace(".bam", f".only_{pc_tag}.bam")))
    raw_masked_bam = masked_bam.replace(".bam", ".raw.bam")
    raw_masked_vcf = masked_bam.replace(".bam", ".raw.vcf.gz")

    cmd = f"bash {bash_utils_hub} quick_check_bam_validity {raw_masked_bam} && \
            [[ {raw_masked_bam} -nt {sd_freads} ]] && \
            [[ {raw_masked_bam} -nt {sd_rreads} ]] && \
            bash {bash_utils_hub} check_vcf_validity {raw_masked_vcf} 1 && \
            [[ {raw_masked_vcf} -nt {raw_masked_bam} ]] && \
            [[ {sd_freads} -nt {pc_bed} ]] && \
            [[ {sd_freads} -nt {original_bam} ]]"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.info(f"The masked bam {raw_masked_bam} is not valid, so we need to generate it.")
        execute = True
    else:
        logger.info(f"The masked bam {raw_masked_bam} is already up-to-date. No need to map again.")
        execute = False

    if execute:
        # Extract the qnames
        pc_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                     pc_fc_bed, 
                                                     output_qnames = "Yes",
                                                     threads=threads,
                                                     multi_aligned=False,
                                                     logger=logger )

        counterpart_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                              pc_nfc_bed, 
                                                              output_qnames = "Yes",
                                                              threads=threads,
                                                              multi_aligned=True,
                                                              mq_cutoff=mq_cutoff,
                                                              logger=logger )

        
        # merged_qname_lst = prepare_tmp_file().name
        qname_lst_f = prepare_tmp_file().name
        qname_lst_r = prepare_tmp_file().name
        cmdf = f"cat {pc_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_f}"
        cmdr = f"cat {pc_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_r}"
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
                -i {pc_tag} \
                -c {max_varno} && ls -lh {raw_masked_bam}"
        try:
            executeCmd(cmd, logger=logger)
        except RuntimeError:
            return f"{pc_tag},{pc_bed},NaN,NaN"
        else:
            # Now perform the variants calling
            # First determine the names of the output VCFs
            logger.info(f"Now we have prepared the masked BAM file {raw_masked_bam} for sample {sample_ID} for {pc_bed}, we need to generate the VCF file.")
            
            sub_running_log = raw_masked_vcf.replace(".vcf.gz", ".log")
            cmd = f"bash {bash_utils_hub} bcftools_call_per_PC \
                    -a {original_bam} \
                    -m {ref_genome} \
                    -c {threads} \
                    -o {raw_masked_vcf} \
                    -p {pc_tag} \
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
                    return f"{pc_tag},{pc_bed},{raw_masked_bam},NaN"
                else:
                    logger.warning(f"Though the process reported error. The VCF file {raw_masked_vcf} is still valid and updated")
                    return f"{pc_tag},{pc_bed},{raw_masked_bam},{raw_masked_vcf}"
            else:
                logger.info(f"Succesfully generate the vcf file {raw_masked_vcf} for region {pc_bed}")
                return f"{pc_tag},{pc_bed},{raw_masked_bam},{raw_masked_vcf}"
    else:
        logger.info(f"The masked bam {raw_masked_bam} is already up-to-date. No need to map again.")
        return f"{pc_tag},{pc_bed},{raw_masked_bam},{raw_masked_vcf}"





def extract_PC_beds_and_genome(total_folder: str):
    import re
    PC_beds = []
    total_region_beds = []
    masked_genomes = []

    for root, dirs, files in os.walk(total_folder):
        if re.search(r"/PC[0-9]+$", root):
            # Extract the PC beds and masked genome
            PC_tag = re.search(r"/(PC[0-9]+)$", root).group(1)
            for file in files:
                if file.endswith(f'{PC_tag}.bed'):
                    PC_beds.append(os.path.join(root, file))
                elif file.endswith(f'{PC_tag}.masked.fasta'):
                    masked_genomes.append(os.path.join(root,file))
        elif re.search(r"/PC[0-9]+_related_homo_regions$", root):
            # Extract the total_regions beds
            whole_region_tag = re.search(r"/(PC[0-9]+_related_homo_regions)$", root).group(1)
            for file in files:
                if file.endswith(f'{whole_region_tag}.bed'):
                    total_region_beds.append(os.path.join(root,file))
    return PC_beds, total_region_beds, masked_genomes



def annotate_intrinsic_vars(query_vcf:str,
                            output_vcf:str,
                            sd_working_dir="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{samp_ID}_refSD_priority_component_pairs",
                            intrinsic_vcf_name="all_pc_region_intrinsic_variants.vcf.gz",
                            conf_level = 0.001):
    
    tmp_output = prepare_tmp_file(suffix=".vcf.gz").name
    intrinsic_vcf = os.path.join(sd_working_dir, intrinsic_vcf_name)
    from compare_with_intrinsic_vcf import main_process
    main_process(query_vcf_path=query_vcf, 
                 intrin_vcf_path=intrinsic_vcf,
                 conf_level=conf_level,
                 output_vcf=tmp_output)
    
    update_gzip_file_on_md5(output_vcf, tmp_output)
    executeCmd(f"tabix -f -p vcf {output_vcf}")
    
    
def annotate_inhouse_common_vars(query_vcf: str,
                                 output_vcf: str,
                                 ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                 inhouse_sdrecall_vcf="/paedyl01/disk1/yangyxt/ngs.SDrecall.control.{ref_gen_tag}.vcf.gz"):
    tmp_output = prepare_tmp_file(suffix=".vcf.gz").name
    ref_gen_tag = os.path.basename(ref_genome).split(".")[1]

    inhouse_sdrecall_vcf = inhouse_sdrecall_vcf.format(ref_gen_tag=ref_gen_tag)
    
    from filter_var_by_inhouse_AF import identify_inhouse_common
    identify_inhouse_common(query_vcf,
                            output_vcf = tmp_output,
                            inhouse_common_cutoff = 0.01,
                            internal_cohort_shortv = inhouse_sdrecall_vcf)

    update_gzip_file_on_md5(output_vcf, tmp_output)
    executeCmd(f"tabix -f -p vcf {output_vcf}")



def stat_fc_nfc_size_per_bed(bed: str):
    pc_tag = os.path.basename(bed).split("_")[0]
    bedf = pd.read_table(bed, header=None)
    fc_bedf = bedf.loc[bedf.iloc[:, -1].str.contains(r"^FC:"), :].drop_duplicates().sort_values(by=bedf.columns.tolist()[-1])
    nfc_bedf = bedf.loc[bedf.iloc[:, -1].str.contains(r"^NFC:"), :].drop_duplicates()
    subgroup_fc_bedfs = [fc_bedf.iloc[i, :] for i in range(0, len(fc_bedf))]
    # display(subgroup_fc_bedfs[0].to_frame().T)
    records = []
    for subgroup_fc_bedf in subgroup_fc_bedfs:
        subgroup_id = subgroup_fc_bedf.iloc[-1].split(":")[-1].split("_")[-1]
        nfc_bool = nfc_bedf.iloc[:, -1].str.contains(fr"^NFC:{pc_tag}_{subgroup_id}$")
        subgroup_nfc_bedf = nfc_bedf.loc[nfc_bool, :]
        fc_subgroup_size = pb.BedTool.from_dataframe(subgroup_fc_bedf.to_frame().T).sort().total_coverage()
        nfc_subgroup_size = pb.BedTool.from_dataframe(subgroup_nfc_bedf).sort().total_coverage()
        
        records.append((pc_tag, subgroup_id, fc_subgroup_size, nfc_subgroup_size))
    return records


def stat_all_PC_subgroups(wkd, threads=10):
    # Test this function, only takes 3 min to finish
    # use glob module to list out all the bed files under wkd named as "PC*_related_homo_regions.raw.bed"
    beds = glob.glob(f"{wkd}/PC*_related_homo_regions/PC*_related_homo_regions.raw.bed")
    logger.info(f"The beds are {beds}")
    
    beds = sorted(beds, key = lambda f: os.path.getsize(f), reverse=True)
    with ctx.Pool(threads) as pool:
        results = pool.imap_unordered(stat_fc_nfc_size_per_bed, beds)
        records = [t for l in results for t in l]
        
    # Convert the list of tuples into a dataframe
    # You need to sort the df by the fc_size (proportional to the downstream function run time) in desending order to ensure the load balancing for subsequent parallel computation
    stat_df = pd.DataFrame(records, columns=["pc_name", "subgroup_id", "fc_size", "nfc_size"])
    stat_df["sorting_index"] = stat_df.loc[:, "fc_size"] * stat_df.loc[:, "nfc_size"]
    stat_df = stat_df.sort_values(by="sorting_index", ascending=False)
    return stat_df



def stat_all_PC_region_size(wkd, threads=10):
    # Test this function, only takes 3 min to finish
    # use glob module to list out all the bed files under wkd named as "PC*_related_homo_regions.raw.bed"
    beds = glob.glob(f"{wkd}/PC*_related_homo_regions/PC*_related_homo_regions.raw.bed")
    logger.info(f"The beds are {beds}")
    
    beds = sorted(beds, key = lambda f: pb.BedTool(f).sort().total_coverage(), reverse=True)

    records = []
    for bed in beds:
        bedf = pd.read_table(bed, header=None)
        fc_bedf = bedf.loc[bedf.iloc[:, -1].str.contains(r"^FC:"), :].drop_duplicates().sort_values(by=bedf.columns.tolist()[-1])
        sub_ids = [i for i in range(fc_bedf.shape[0])]
        pc_name = os.path.basename(bed).split("_")[0]
        records.append((pc_name,sub_ids))

    df = pd.DataFrame(records, columns=["pc_name", "subgroup_id"])
    df_exploded = df.explode("subgroup_id")
    return df_exploded




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


def calculate_NM_distribution_poisson(bam, conf_level=0.01, sample_size=1000000, histo_fig=None, logger=logger):
    '''
    The input vcf_rec is simply a tuple containing(chr, start, ref, alt), remember the coordinates from VCF are 1-indexed
    '''
    nm_array = np.empty(sample_size)
    bam_handle = pysam.AlignmentFile(bam, "rb")
    n = 0
    for read in bam_handle.fetch():
        if read.is_proper_pair and \
           not read.is_secondary and \
           not read.is_supplementary and \
           not read.is_duplicate and \
           not read.is_qcfail and \
           len([t for t in read.get_tags() if t[0] in ["XA"]]) == 0 and \
           read.mapping_quality == 60:
            gap_sizes = [t[1] for t in read.cigartuples if t[0] in [1,2] and t[1] > 1]
            max_gap_size = max(gap_sizes) if len(gap_sizes) > 0 else 0
            edit_dist = read.get_tag("NM")
            scatter_edit_dist = edit_dist - max_gap_size
            nm_array[n] = scatter_edit_dist
            n += 1
            if n >= sample_size:
                break

    bam_handle.close()

    # Deal with the edge case where the input BAM does not have enough reads
    if n < sample_size:
        logger.warning(f"Input BAM {bam} does not have enough reads, only {n} reads are used for NM distribution calculation")
        nm_array = nm_array[:n]
    
    nm_mean = np.mean(nm_array)
    cutoff = 0
    while stats.poisson.cdf(cutoff, nm_mean) < 1 - conf_level:
        cutoff += 1
    logger.info(f"NM distribution for {bam} has mean {nm_mean}, the value for {1-conf_level} percentile is {cutoff}")

    if histo_fig:
        logger.info(f"Plotting NM distribution histogram for {bam} to Figure {histo_fig}")
        import matplotlib.pyplot as plt
        num_bins = 100
        n, bins, patches = plt.hist(nm_array, num_bins, density=True, alpha=0.6, color='g', edgecolor='black')
        # Add a best fit Poisson distribution line
        xmin, xmax = plt.xlim()
        x = np.arange(xmin, xmax)
        p = stats.poisson.pmf(x, nm_mean)
        plt.plot(x, p, 'k', linewidth=2)

        title = "NM distribution of %s" % bam
        plt.title(title)
        plt.savefig(histo_fig, dpi=300)

    return cutoff, nm_mean


    

def calculate_NM_distribution(bam, conf_level=0.01, sample_size = 1000000, histo_fig = None, logger = logger):
    '''
    The input vcf_rec is simply a tuple containing(chr, start, ref, alt), remember the coordinates from VCF are 1-indexed
    '''
    nm_array = np.empty(sample_size)
    bam_handle = pysam.AlignmentFile(bam, "rb")
    n = 0
    for read in bam_handle.fetch():
        if read.is_proper_pair and \
           not read.is_secondary and \
           not read.is_supplementary and \
           not read.is_duplicate and \
           not read.is_qcfail and \
           len([t for t in read.get_tags() if t[0] in ["XA"]]) == 0 and \
           read.mapping_quality == 60:
            gap_sizes = [t[1] for t in read.cigartuples if t[0] in [1,2] and t[1] > 1]
            max_gap_size = max(gap_sizes) if len(gap_sizes) > 0 else 0
            edit_dist = read.get_tag("NM")
            scatter_edit_dist = edit_dist - max_gap_size
            nm_array[n] = scatter_edit_dist
            n += 1
            if n >= sample_size:
                break

    bam_handle.close()

    # Deal with the edge case where the input BAM does not have enough reads
    if n < sample_size:
        logger.warning(f"Input BAM {bam} does not have enough reads, only {n} reads are used for NM distribution calculation")
        nm_array = nm_array[:n]
    
    nm_mean = np.mean(nm_array)
    nm_std = np.std(nm_array)
    cutoff = stats.norm.ppf(1-conf_level, loc=nm_mean, scale=nm_std)
    logger.info(f"NM distribution for {bam} has mean {nm_mean} and std {nm_std}, the value for {1-conf_level} percentile is {cutoff}")

    if histo_fig:
        logger.info(f"Plotting NM distribution histogram for {bam} to Figure {histo_fig}")
        import matplotlib.pyplot as plt
        num_bins = 100
        n, bins, patches = plt.hist(nm_array, num_bins, density=True, alpha=0.6, color='g', edgecolor='black')
        # Add a best fit normal distribution line
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = stats.norm.pdf(x, nm_mean, nm_std)
        plt.plot(x, p, 'k', linewidth=2)

        title = "NM distribution of %s" % bam
        plt.title(title)
        plt.savefig(histo_fig, dpi= 300)

    return cutoff
    

'''
I have no idea that why I always run into this error: 

OpenBLAS blas_thread_init: pthread_create failed for thread 3 of 80: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 3089471 max
OpenBLAS blas_thread_init: pthread_create failed for thread 4 of 80: Resource temporarily unavailable
OpenBLAS blas_thread_init: pthread_create failed for thread 3 of 80: Resource temporarily unavailable
OpenBLAS blas_thread_init: pthread_create failed for thread 4 of 80: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 3089471 max
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 3089471 max
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 3089471 max
OpenBLAS blas_thread_init: pthread_create failed for thread 12 of 80: Resource temporarily unavailable
OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 3089471 max

I run current script in 6 concurrency, and each process will use 12 threads, so in total 72 threads are used.
But I have in total 80 CPUs at disposal

I have tried to set the environment variable OMP_NUM_THREADS=1 and many other env variables at the beg, but it does not work
Now I try to migrate the post_process_bam function to a bash function and migrate the parallelism part to GNU parallel instead of python multiprocessing
As I recall, I used to run GNU parallel for filter_out_noisy_reads function and it does not have this problem

Update:
It turns out bcftools will initiate an openblas thread for every input vcf file. So when I try to merge over 6000 files to 6 file. The openblas threads hit limit.
'''


def imap_slice_bam_per_bed(tup_args):
    return slice_bam_per_bed(*tup_args)



@log_command
def slice_bam_per_bed(bed, bam, ref_genome, threads = 4, logger = logger):
    # assert os.path.basename(bed).split(".")[-2] != chunk_id, f"The bed file {bed} should not have the same chunk_id as the previous bed file"
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
    
    cmd = f"bash {bash_utils_hub} quick_check_bam_validity {remapped_bam} && \
            bash {bash_utils_hub} quick_check_bam_validity {remapped_clean_bam} && \
            bash {bash_utils_hub} quick_check_bam_validity {remapped_noise_bam} && \
            bash {bash_utils_hub} check_vcf_validity {clean_vcf} && \
            [[ {clean_vcf} -nt {remapped_clean_bam} ]] && \
            [[ {remapped_clean_bam} -nt {remapped_bam} ]] && \
            [[ {remapped_bam} -nt {raw_bam} ]] && \
            [[ {remapped_clean_bam} -nt /paedyl01/disk1/yangyxt/ngs_scripts/gt_read_clustering_filter.py ]] && \
            [[ {clean_vcf} -nt /paedyl01/disk1/yangyxt/ngs_scripts/gt_read_clustering_filter.py ]] && false"

    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.info(f"Start to post process the bam file {raw_bam}")
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
        parallelism = 4 if parallelism < 4 else parallelism
        
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
    else:
        logger.info(f"The post processed masked bam and vcf files are already up-to-date. No need to run the post processing again.")


    
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
    '''
    Deprecated in 13/05/2024
    vis_dict = filter_out_noisy_reads(raw_bam,
                                        output_bam = masked_bam, 
                                        filter_out_bam = noisy_bam,
                                        intrinsic_bam = intrinsic_bam,
                                        conf_level = conf_level,
                                        sample_size = sample_size,
                                        histo_fig = histo_fig,
                                        threads = threads,
                                        varno_cutoff = varno_cutoff,
                                        max_varno = max_varno,
                                        logger = logger)
    '''

    try:
        cmd = f"[[ -f {masked_bam} ]] && \
                [[ {masked_bam} -nt {raw_bam} ]] && \
                [[ -f {noisy_bam} ]] && \
                [[ {noisy_bam} -nt {raw_bam} ]]"
        executeCmd(cmd, logger = logger)
    except RuntimeError:
        execute = True
    else:
        execute = False

    if execute:
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
    cmd = f"bash {bash_utils_hub} bcftools_call_per_PC \
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
    


def get_bam_frag_size(input_bam):
    cmd = f"samtools stats {input_bam}"
    output = executeCmd(cmd, stdout_only=True)
    for line in output.split('\n'):
        if line.startswith('SN'):
            fields = line.split('\t')
            if fields[1].startswith('insert size average:'):
                avg_insert_size = float(fields[2])
            elif fields[1].startswith('insert size standard deviation:'):
                std_insert_size = float(fields[2])
    return avg_insert_size, std_insert_size



def SDrecall_per_sample(input_bam: str, 
                        input_freads: str, 
                        input_rreads: str, 
                        output_vcf: str,
                        threads = 12,
                        numba_threads = 4,
                        mq_cutoff = 20,
                        conf_level = 0.01,
                        varno_cutoff = 3,
                        stat_sample_size = 3000000,
                        histo_fig = None,
                        target_region = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.coding.func.pad20.bed",
                        all_PC_folder = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{samp_ID}_refSD_priority_component_pairs",
                        ref_fasta = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                        profile_file = None):
    set_num_threads(numba_threads)
    logger.info(f"The NUMBA threads is set to {get_num_threads()}")

    if profile_file:
        pr = cProfile.Profile()
        pr.enable()

    logger.info(f"The threads deployed to this function is {threads}")
    avg_frag_size, std_frag_size = get_bam_frag_size(input_bam)
    # First we need to extract the PCs and the total regions into two lists.
    # Then we need to perform the reads pooling per PC-LPC pairs, and we can do it in mutliple subprocesses
    # The bcftools_call_per_PC function contains several steps:
    # 1. Extract the reads from input BAM file at input target regions, into a pair of FASTQ files
    # 2. map the extracted FASTQ files against the masked genome
    # 3. Call GVCF on the masked genome BAM file
    # 4. Genotype GVCF on the GVCF file
    # depth_tab = calculate_inferred_coverage(input_bam, min_mapq = 0)
    output_dir = os.path.dirname(output_vcf)
    samp_ID = os.path.basename(input_bam).split(".")[0]
    all_PC_folder = all_PC_folder.format(samp_ID=samp_ID)
    intrinsic_bam = os.path.join(all_PC_folder, "total_intrinsic_alignments.bam")
    prepared_arguments_df = stat_all_PC_region_size(all_PC_folder, threads=threads)
    gc.collect()

    fc_size_stat_tab = os.path.join(all_PC_folder, "FC_size_stat.tsv")
    logger.info(f"The prepared arguments for the PC subgroups will be saved to {fc_size_stat_tab} looks like:\n{prepared_arguments_df[:20].to_string(index=False)}\n\n")
    prepared_arguments_df.to_csv(fc_size_stat_tab, sep="\t", index=False)
    
    uniq_pc_names = prepared_arguments_df.loc[:, "pc_name"].drop_duplicates().tolist()
    pc_subids_tup_list = [tuple(prepared_arguments_df.loc[prepared_arguments_df["pc_name"] == pc, "subgroup_id"].drop_duplicates().tolist()) for pc in uniq_pc_names]
    
    pool = ctx.Pool(threads, initializer=pool_init)
    prepared_beds = pool.imap_unordered(imap_prepare_masked_align_region_per_PC, zip(uniq_pc_names,
                                                                                    pc_subids_tup_list,
                                                                                    repeat(all_PC_folder),
                                                                                    repeat(target_region),
                                                                                    repeat(ref_fasta)))
    
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

    prepared_beds_df = pd.read_csv(io.StringIO("\n".join(result_beds)), header=None, names=["pc_index", "pc_subgroup_id", "fc_bed", "nfc_bed", "fc_bed_size", "nfc_bed_size"], na_values="NaN")
    logger.info(f"The prepared beds for the masked alignment looks like:\n{prepared_beds_df[:20].to_string(index=False)}\n\n")
    uniq_pcs = prepared_beds_df.loc[:, "pc_index"].drop_duplicates().tolist()
    uniq_pcs = sorted(uniq_pcs, key = lambda x: prepared_beds_df.loc[prepared_beds_df["pc_index"] == x, "nfc_bed_size"].sum(), reverse=True)
    
    # Merge bam files, note that bam file list will be too long to be passed as a command line argument
    # This is for data visualization and debug purpose
    nm_cutoff, nm_mean = calculate_NM_distribution_poisson(input_bam, conf_level, stat_sample_size, histo_fig, logger)
    pc_specific_regions = [ os.path.join(all_PC_folder, f"{pc}_related_homo_regions", f"{pc}", f"{pc}.bed") for pc in uniq_pcs ]
    pc_specific_fc_beds = [ prepared_beds_df.loc[prepared_beds_df["pc_index"] == pc, "fc_bed"].tolist() for pc in uniq_pcs ]
    pc_specific_nfc_beds = [ prepared_beds_df.loc[prepared_beds_df["pc_index"] == pc, "nfc_bed"].tolist() for pc in uniq_pcs]
    # pc_specific_masked_genomes = [ glob.glob(f"{all_PC_folder}/{pc}_related_homo_regions/{pc}/{pc}_*.fc.fasta") for pc in uniq_pcs ]
    
    # unzip the gzipped input_freads will cause 10 fold of extra time in running so we need to create a decompressed version of the fastq files first for downstream analysis
    if input_freads.endswith(".gz"):
        decomp_if = input_freads.replace(".gz", "")
        cmd = f"gunzip -c {input_freads} > {decomp_if} && ls -lht {decomp_if}"
        executeCmd(cmd, logger=logger)
        input_freads = decomp_if
    
    if input_rreads.endswith(".gz"):
        decomp_ir = input_rreads.replace(".gz", "")
        cmd = f"gunzip -c {input_rreads} > {decomp_ir} && ls -lht {decomp_ir}"
        executeCmd(cmd, logger=logger)
        input_rreads = decomp_ir

    pool = ctx.Pool(threads)
    results = pool.imap_unordered(imap_process_masked_bam, zip( uniq_pcs,
                                                                pc_specific_fc_beds,
                                                                pc_specific_nfc_beds,
                                                                repeat(input_bam),
                                                                pc_specific_regions,
                                                                repeat(samp_ID),
                                                                repeat(input_freads),
                                                                repeat(input_rreads),
                                                                repeat(nm_cutoff),
                                                                repeat(1),
                                                                repeat(mq_cutoff),
                                                                repeat(ref_fasta) ))

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
    
    process_meta_table = os.path.join(output_dir, f"{samp_ID}_process_raw_bam_results.tsv")
    post_result_df = pd.read_csv(io.StringIO("\n".join(result_records)), header=None, names=["pc_index", "pc_bed", "raw_masked_bam", "raw_masked_vcf"], na_values="NaN")
    post_result_df.to_csv(process_meta_table, sep="\t", index=False)
    logger.info(f"The processed masked bam and vcf files stored in {process_meta_table} and it looks like:\n{post_result_df.to_string(index=False)}\n\n")
    

    # Merge the BAM files
    # bam_regex = r"^(?P<name>.+)\.bam$"
    # valid_raw_bams = post_result_df.loc[post_result_df["masked_bam"].notna(), "masked_bam"].str.replace(bam_regex, lambda m:m.group("name") + ".raw.bam", regex=True).tolist()
    # valid_clean_bams = post_result_df.loc[post_result_df["masked_bam"].notna(), "masked_bam"].str.replace(bam_regex, lambda m:m.group("name") + ".clean.bam", regex=True).tolist()
    remapped_raw_bam = input_bam.replace(".bam", ".pooled.raw.bam")
    # remapped_bam = input_bam.replace(".bam", ".pooled.bam")
    # remapped_clean_bam = input_bam.replace(".bam", ".pooled.clean.bam")
    execute_merging = merge_bams(post_result_df.loc[:, "raw_masked_bam"].dropna().drop_duplicates().tolist(), remapped_raw_bam, ref_fasta=ref_fasta, threads=threads)
    deduped_raw_bam = remapped_raw_bam.replace(".bam", ".deduped.bam")
    cmd = f"samtools collate -@ {threads} -O -u {remapped_raw_bam} | \
            samtools fixmate -@ {threads} -m -u - - | \
            samtools sort -@ {threads} -u - | \
            samtools markdup -@ {threads} -r - {deduped_raw_bam} && \
            ls -lh {deduped_raw_bam}"

    if execute_merging:
        executeCmd(cmd, logger=logger)
    
    # Drop dup and redundant sequence in total_intrin_bam file
    # deduped_intrinsic_bam = intrinsic_bam.replace(".bam", ".deduped.bam")
    # drop_dup_read_intrin_bam(intrinsic_bam, deduped_intrinsic_bam, threads=threads, logger=logger)
    
    # merge_bams(post_result_df.loc[:, "masked_bam"].dropna().tolist(), remapped_bam, ref_fasta=ref_fasta, threads=threads)
    # merge_bams(valid_clean_bams, remapped_clean_bam, ref_fasta=ref_fasta, threads=threads)
    
    # Merge VCF files
    # vcf_regex = r"^(?P<name>.+)\.vcf\.gz$"
    # raw_vcfs = post_result_df.loc[post_result_df["masked_vcf"].notna(), "masked_vcf"].str.replace(vcf_regex, lambda m:m.group("name") + ".raw.vcf.gz", regex=True).tolist()
    # unfiltered_vcfs = post_result_df.loc[post_result_df["masked_vcf"].notna(), "masked_vcf"].str.replace(vcf_regex, lambda m:m.group("name") + ".unfiltered.vcf.gz", regex=True).tolist()
    all_pc_intrin_vcfs = [os.path.join(all_PC_folder, f"{pc}_related_homo_regions", f"{pc}", f"{pc}.raw.vcf.gz") for pc in uniq_pcs]
    total_intrin_vcf = os.path.join(all_PC_folder, "total_intrinsic_variants.vcf.gz")
    total_raw_vcf = output_vcf.replace(".vcf.gz", ".raw.vcf.gz")
    
    merge_vcfs(post_result_df.loc[:, "raw_masked_vcf"].dropna().drop_duplicates().tolist(), total_raw_vcf, threads = threads, check_error=False )
    merge_vcfs(all_pc_intrin_vcfs, total_intrin_vcf, threads = threads, check_error=False)
    gc.collect()
    # merge_vcfs(unfiltered_vcfs, output_vcf.replace(".vcf.gz", ".unfiltered.vcf.gz"))
    
    '''
    After Merging the raw BAMs and raw VCFs, we can perform variant filtering by chromosome on the merged BAM 
    '''
    post_process_bam(deduped_raw_bam,
                     total_raw_vcf,
                     output_vcf,
                     input_bam,
                     intrinsic_bam,
                     threads = threads,
                     numba_threads = numba_threads, 
                     conf_level = conf_level,
                     stat_sample_size = stat_sample_size,
                     histo_fig = histo_fig,
                     varno_cutoff = varno_cutoff,
                     max_varno = nm_cutoff,
                     ref_genome = ref_fasta,
                     avg_frag_size = avg_frag_size,
                     logger = logger)
    gc.collect()
    
    # cmd = f"bash {hs_bash_utils} \
    #         post_process_bam \
    #         -r {remapped_raw_bam} \
    #         -o {output_vcf} \
    #         -a {input_bam} \
    #         -i {intrinsic_bam} \
    #         -t {threads} \
    #         -c {conf_level} \
    #         -s {stat_sample_size} \
    #         -f {histo_fig} \
    #         -v {varno_cutoff} \
    #         -m {nm_cutoff} \
    #         -g {ref_fasta}"

    # executeCmd(cmd, logger=logger)
    
    '''
    Below there is a section moved from each bcftools_call_per_PC subprocess to here for the improved computation efficiency
    1. Compare with intrinsic vcf to add LIKELY_INTRINSIC FILTER tag to total VCF file
    2. Compare with raw vcf to add MISALIGNED FILTER tag to total VCF file
    '''
    
    # compare_with_intrinsic_vcf(query_vcf = output_vcf.replace(".vcf.gz", ".raw.vcf.gz"),
    #                            intrinsic_vcf = total_intrin_vcf,
    #                            output_vcf = output_vcf.replace(".vcf.gz", ".raw.intrinsic.vcf.gz"),
    #                            conf_level = conf_level, 
    #                            threads = threads)
    
    compare_vcf_recs(query_vcf = output_vcf.replace(".vcf.gz", ".raw.vcf.gz"),
                     reference_vcf = output_vcf.replace(".vcf.gz", ".clean.vcf.gz"),
                     added_filter = "MISALIGNED",
                     qv_tag = "RAW", 
                     rv_tag = "CLEAN",
                     output_vcf = output_vcf, 
                     threads = threads)
    
    # Then we need to annotate intrinsic variants in the output vcfs
    tmp_vcf = output_vcf.replace(".vcf.gz", ".tmp.vcf.gz")
    assert tmp_vcf != output_vcf
    executeCmd(f"mv {output_vcf} {tmp_vcf} && tabix -f -p vcf {tmp_vcf}")
    # annotate_intrinsic_vars(output_vcf, filtered_vcf, sd_working_dir=all_PC_folder)
    # logger.info("The merged VCF file with annotations on intrinsic variants for sample {} is {}".format(input_bam, filtered_vcf))
    
    # Then we need to tag out the variants that might be common according to in-house SD variants database
    annotate_inhouse_common_vars(tmp_vcf, output_vcf, ref_genome=ref_fasta)
    executeCmd(f"bash {bash_utils_hub} display_table {output_vcf} 20")
    logger.info(f"The merged VCF annotated with potential in-house common variants is named after {output_vcf}")
    
    if profile_file:
        pr.disable()
        pr.dump_stats(profile_file)
    
    

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--function", type=str, help="The function name", required=True)
    parser.add_argument("-a", "--arguments", type=str, help="The function's input arguments, delimited by semi-colon ;", required=False, default=None)
    parser.add_argument("-k", "--key_arguments", type=str, help="Keyword arguments for the function, delimited by semi-colon ;", required=False, default=None)
    
    args = parser.parse_args()
    try:
        fargs = [ convert_input_value(a) for a in args.arguments.split(";") ] if type(args.arguments) == str else []
        fkwargs = { t.split("=")[0]: convert_input_value(t.split("=")[1]) for t in args.key_arguments.split(";") } if type(args.key_arguments) == str else {}
        logger.info("Running function: {}, input args are {}, input kwargs are {}".format(args.function, fargs, fkwargs))
    except Exception as e:
        logger.error("Input argument does not meet the expected format, encounter Parsing error {}, Let's check the input:\n-f {}, -a {}, -k {}".format(
            e,
            args.function,
            args.arguments,
            args.key_arguments
        ))
        raise e
    globals()[args.function](*fargs, **fkwargs)
    
    
    
    
    
    
    
    


