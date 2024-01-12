#!/usr/bin/env python
import os
import io
import subprocess
import pandas as pd
import numpy as np
import logging
import re
import uuid
import multiprocessing as mp
from itertools import repeat
import argparse as ap
from python_utils import convert_input_value, calculate_inferred_coverage, retrieve_infer_depth_bam
import inspect
import time
from pathlib import Path
from io import StringIO
import traceback
import sys
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



def executeCmd(cmd, stdout_only = False, shell="/home/yangyxt/anaconda3/bin/bash", logger=logger) -> None:
    from subprocess import PIPE
    logger.info("About to run this command in shell invoked within python: \n{}\n".format(cmd))
    
    if stdout_only:
        result = subprocess.run(cmd, shell=True, executable=shell, stdout=PIPE, stderr=PIPE)
    else:
        result = subprocess.run(cmd, shell=True, executable=shell, stderr=subprocess.STDOUT, stdout=PIPE)
        
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        logger.error("Error in \n{}\nAnd the output goes like:\n{}\n".format(" ".join(cmd_lst), result.stdout.decode()))
        if cmd_lst[1][0] != "-":
            raise RuntimeError
        else:
            raise RuntimeError
        
    logger.info(f"Running the following shell command inside python:\n{cmd}\nAnd it receives a return code of {code}, the output goes like this:\n{result.stdout.decode()}\n\n")
    return result.stdout.decode()


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
        cmd = f"""samtools view -@ {threads} -L {region_bed} {input_bam} | mawk -F '\t' '($0 !~ /SA:Z:/) && (($5 <= {mq_cutoff + 10}) || ($0 ~ /XA:Z:/)) {{print $1;}}' | tail -n +1 | sort - | uniq - > {tmp_file}"""
    else:
        cmd = f"""samtools view -@ {threads} -L {region_bed} {input_bam} | cut -f 1 | tail -n +1 | sort - | uniq - > {tmp_file}"""
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


def extract_reads_by_qnames(qname_lst, 
                            input_freads, input_rreads, 
                            output_freads, output_rreads, 
                            logger = logger):
    cmdf = f"seqtk subseq {input_freads} {qname_lst} > {output_freads}"
    cmdr = f"seqtk subseq {input_rreads} {qname_lst} > {output_rreads}"
    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)
    # Make sure the output fastq files are consistent regarding the read qnames
    cmd = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
    executeCmd(cmd, logger=logger)
    # os.remove(tmp_file)
    return output_freads, output_rreads



def estimate_avg_depth(input_bam: str, region_bed: str, MQ_filter=False) -> float:
    if MQ_filter:
        cmd = f"bash {bash_utils_hub} extract_avg_depth_with_highqual {input_bam} {region_bed}"
    else:
        cmd = f"bash {bash_utils_hub} extract_avg_depth_without_MQ {input_bam} {region_bed}"
    return executeCmd(cmd, stdout_only=True)



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



def imap_prepare_masked_align(tup_args):
    return prepare_masked_align_per_PC(*tup_args)


@log_command
def prepare_masked_align_per_PC(pc_bed: str,
                                whole_region_bed: str,
                                masked_genome: str,
                                input_bam: str,
                                freads: str,
                                rreads: str,
                                mq_cutoff = 20,
                                threads=2,
                                logger = logger):
    pc_tag = os.path.basename(pc_bed).split(".")[0]
    sample_ID = str(os.path.basename(input_bam)).split(".")[0]
    per_sample_pc_dir = os.path.join(os.path.dirname(input_bam), sample_ID + "_ref_SD_recall")
    
    try:
        os.mkdir(per_sample_pc_dir)
    except FileExistsError:
        pass
    
    # First extract the reads from the input freads and rreads
    sd_freads = os.path.join(per_sample_pc_dir, os.path.basename(freads.replace(".fastq.gz", f".{pc_tag}_sdrecall.fastq")))
    sd_rreads = os.path.join(per_sample_pc_dir, os.path.basename(rreads.replace(".fastq.gz", f".{pc_tag}_sdrecall.fastq")))
    
    assert sd_rreads != rreads
    assert sd_freads != freads
    
    test_cmd = f"bash {bash_utils_hub} check_fastq_pair_consistency {sd_freads} {sd_rreads}"
    try:
        executeCmd(test_cmd, logger=logger)
    except RuntimeError:
        logger.info(f"The extracted fastqs {sd_freads} and {sd_rreads} for {pc_bed} are not consistent regarding qnames between them. We need to extract again.")
        execute = True
    else:
        logger.info(f"The extracted fastqs {sd_freads} and {sd_rreads} for {pc_bed} are a pair of fastqs having consistent query names")
        if  os.path.getmtime(sd_freads) > os.path.getmtime(freads) and \
            os.path.getmtime(sd_rreads) > os.path.getmtime(rreads) and \
            os.path.getmtime(sd_freads) > os.path.getmtime(whole_region_bed) and \
            os.path.getmtime(sd_freads) > os.path.getmtime(input_bam) and \
            os.path.getmtime(sd_rreads) > os.path.getmtime(input_bam):
                logger.info(f"The extracted fastqs {sd_freads} and {sd_rreads} for {pc_bed} are already up-to-date. No need to extract again.")
                execute = False
        else:
            execute = True
        
    if execute:
        pc_qname_lst = bam_reads_selection_by_region( input_bam, 
                                                      pc_bed, 
                                                      output_qnames = "Yes",
                                                      threads=threads,
                                                      multi_aligned=False,
                                                      logger=logger )

        counterpart_qname_lst = bam_reads_selection_by_region( input_bam, 
                                                                whole_region_bed, 
                                                                output_qnames = "Yes",
                                                                threads=threads,
                                                                multi_aligned=True,
                                                                mq_cutoff=mq_cutoff + 10,
                                                                logger=logger )

        merged_qname_lst = prepare_tmp_file().name
        cmd = f"cat {pc_qname_lst} {counterpart_qname_lst} | sort - | uniq > {merged_qname_lst}"
        executeCmd(cmd, logger=logger)

        sd_freads, sd_rreads = extract_reads_by_qnames( merged_qname_lst, 
                                                        freads, rreads, 
                                                        sd_freads, sd_rreads, 
                                                        logger=logger )
    
    # Now we use the extracted fastqs to map to masked genome.
    masked_bam = os.path.join(per_sample_pc_dir, os.path.basename(input_bam.replace(".bam", f".only_{pc_tag}.bam")))
    clean_masked_bam = os.path.join(per_sample_pc_dir, os.path.basename(input_bam.replace(".bam", f".only_{pc_tag}.clean.bam")))
    logger.info(f"The masked bam file (realigned bam file) is {masked_bam}, the reads used to generate this masked bam file are {sd_freads} and {sd_rreads}")
    
    test_cmd = f"bash {bash_utils_hub} check_bam_validity {masked_bam} && bash {bash_utils_hub} check_bam_mapping {masked_bam}"
    try:
        executeCmd(test_cmd, logger=logger)
    except RuntimeError:
        logger.info(f"The masked bam {masked_bam} is not valid, so we need to generate it.")
    else:
        if os.path.getmtime(masked_bam) > os.path.getmtime(sd_freads) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(sd_rreads) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(pc_bed) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(masked_genome) and \
           os.path.getmtime(masked_genome) > os.path.getmtime(whole_region_bed) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(input_bam) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(bash_utils_hub) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(script_path) and \
           os.path.getmtime(masked_bam.replace(".bam", ".clean.bam")) > os.path.getmtime(masked_bam) and \
           os.path.getmtime(masked_bam) > os.path.getmtime(masked_bam.replace(".bam", ".raw.bam")):
            logger.info(f"The masked bam {masked_bam} is already up-to-date. No need to map again.")
            return f"{pc_bed},{sd_freads},{sd_rreads},{clean_masked_bam},{input_bam},{whole_region_bed},{masked_genome}"

    cmd = f"bash {bash_utils_hub} independent_minimap2_masked -b {input_bam} -a {masked_genome} -s {sample_ID} -f {sd_freads} -r {sd_rreads} -o {masked_bam} -t {threads} -i {pc_tag} && ls -lh {masked_bam}"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        return f"{pc_bed} failed"
    else:
        test_cmd = f"bash {bash_utils_hub} check_bam_mapping {masked_bam}"
        try:
            executeCmd(test_cmd, logger=logger)
        except RuntimeError:
            raise ValueError(f"{masked_bam} has more unmapped reads than mapped reads. Try run command python3 {script_path} -f prepare_masked_align_per_PC -k \"pc_bed={pc_bed};whole_region_bed={whole_region_bed};masked_genome={masked_genome};input_bam={input_bam};freads={freads};rreads={rreads};threads={threads}\"")
        else:
            return f"{pc_bed},{sd_freads},{sd_rreads},{clean_masked_bam},{input_bam},{whole_region_bed},{masked_genome}"
    


def imap_call_polyploidy(tup_args):
    return call_polyploidy_per_PC(*tup_args)


@log_command
def call_polyploidy_per_PC (freads: str,
                            rreads: str,
                            masked_bam: str,
                            raw_bam: str,
                            target_region_bed: str,
                            masked_genome: str,
                            pc_bed: str,
                            output_vcf = None,
                            threads = 3,
                            ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                            logger = logger):
    if not output_vcf:
        output_vcf = masked_bam.replace(".bam", ".vcf.gz")
    
    sample_name = str(os.path.basename(output_vcf.split(".")[0]))
    output_dir = os.path.dirname(output_vcf)
    
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
        
    running_log = re.sub(r"\.vcf.*", ".log", output_vcf)
    pc_name = str(os.path.basename(pc_bed)).replace(".bed", "")
    pc_vcf_name = os.path.basename(output_vcf.replace(".vcf", f".{pc_name}.vcf"))
    pc_vcf = os.path.join(output_dir, pc_vcf_name)
    sub_running_log = running_log.replace(".log", f".{pc_name}.log")
    reassembly_out_bam = sub_running_log.replace(".log", ".realigned.bam")
    reassembly_filtered_bam = reassembly_out_bam.replace(".bam", ".filtered.bam")
    masked_filtered_bam = masked_bam.replace(".bam", ".filtered.bam")
    pc_dir = os.path.dirname(pc_bed)
    pc_tag = os.path.basename(pc_bed.replace(".bed", ""))
    pc_intrin_vcf = None
    for root, dirs, files in os.walk(pc_dir):
        for file in files:
            if re.search(r"{}\.[0-9]+\.vcf\.gz$".format(pc_tag), file):
                if not pc_intrin_vcf:
                    pc_intrin_vcf = os.path.join(pc_dir, file)
                else:
                    if os.path.getmtime(pc_intrin_vcf) < os.path.getmtime(os.path.join(pc_dir, file)):
                        pc_intrin_vcf = os.path.join(pc_dir, file)
    
    test_cmd = f"bash {bash_utils_hub} check_vcf_validity {pc_vcf} && [[ $(bash {bash_utils_hub} check_vcf_ploidy {pc_vcf}) -eq 2 ]] && bash {bash_utils_hub} check_vcf_multiallelics {pc_vcf}"
    

    try:
        executeCmd(test_cmd, logger=logger)
    except RuntimeError:
        logger.info(f"{pc_vcf} is not ready")
    else:
        if os.path.exists(reassembly_filtered_bam) and os.path.exists(pc_vcf) and os.path.exists(pc_bed):
            if os.path.getmtime(pc_vcf) > os.path.getmtime(masked_bam) and \
               os.path.getmtime(pc_vcf) > os.path.getmtime(target_region_bed) and \
               os.path.getmtime(pc_vcf) > os.path.getmtime(freads) and \
               os.path.getmtime(reassembly_filtered_bam) > os.path.getmtime(masked_bam) and \
               os.path.getmtime(pc_bed) > os.path.getmtime(target_region_bed) and \
               os.path.getmtime(pc_vcf) > os.path.getmtime(script_path) and \
               os.path.getmtime(pc_vcf) > os.path.getmtime(bash_utils_hub):
                logger.info(f"The recall VCF {pc_vcf} for region {pc_bed} is valid and updated. No need to recall the VCF file again.")
                return f"{pc_vcf},{pc_bed},{masked_filtered_bam},{reassembly_filtered_bam}"

    cmd = f"bash {bash_utils_hub} call_polyploidy_per_PC \
                    -s {pc_bed} \
                    -f {freads} \
                    -r {rreads} \
                    -a {raw_bam} \
                    -t {target_region_bed} \
                    -m {ref_genome} \
                    -c {threads} \
                    -o {pc_vcf} \
                    -g {reassembly_out_bam} \
                    -i {pc_intrin_vcf} \
                    -b {masked_bam} > {sub_running_log} 2>&1"
        # We need to allow this function to be failed sometimes
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError as rexception:
        logger.warning(f"Failed to generate {pc_vcf}, running log in {sub_running_log}")
        cmd = f"bash {bash_utils_hub} check_vcf_validity {pc_vcf} 1 && [[ {pc_vcf} -nt {script_path} ]]"
        try:
            executeCmd(cmd, logger=logger)
        except RuntimeError:
            return f"{pc_bed} failed"
        else:
            logger.warning(f"Though the process reported error. The VCF file {pc_vcf} is still valid and updated")
            return f"{pc_vcf},{pc_bed},{masked_filtered_bam},{reassembly_filtered_bam}"
    else:
        logger.info(f"Succesfully generate the vcf file {pc_vcf} for region ${pc_bed}")
        return f"{pc_vcf},{pc_bed},{masked_filtered_bam},{reassembly_filtered_bam}"

    
    
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
                            conf_level = 0.01):
    
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
                                 inhouse_sdrecall_vcf="/paedyl01/disk1/yangyxt/ngs.SDrecall.control.vcf.gz"):
    tmp_output = prepare_tmp_file(suffix=".vcf.gz").name
    
    from filter_var_by_inhouse_AF import identify_inhouse_common
    identify_inhouse_common(query_vcf,
                            output_vcf = tmp_output,
                            inhouse_common_cutoff = 0.01,
                            internal_cohort_shortv = inhouse_sdrecall_vcf)

    update_gzip_file_on_md5(output_vcf, tmp_output)
    executeCmd(f"tabix -f -p vcf {output_vcf}")


    
def SDrecall_per_sample(input_bam: str, input_freads: str, input_rreads: str, output_vcf: str,
                        threads = 12,
                        mq_cutoff = 20,
                        target_region = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.coding.func.pad20.bed",
                        all_PC_folder = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{samp_ID}_refSD_priority_component_pairs",
                        ref_fasta = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"):
    all_PC_region_bed = os.path.join(all_PC_folder, "all_PC_related_homo_regions.bed")
    
    logger.info(f"The threads deployed to this function is {threads}")
    # First we need to extract the PCs and the total regions into two lists.
    PC_beds, total_region_beds, masked_genomes = extract_PC_beds_and_genome(all_PC_folder)
    
    # Then we need to perform the reads pooling per PC-LPC pairs, and we can do it in mutliple subprocesses
    # The call_poly_ploidy_per_PC function contains several steps:
    # 1. Extract the reads from input BAM file at input target regions, into a pair of FASTQ files
    # 2. map the extracted FASTQ files against the masked genome
    # 3. Call GVCF on the masked genome BAM file
    # 4. Genotype GVCF on the GVCF file
    depth_tab = calculate_inferred_coverage(input_bam, min_mapq = 0)
    
    pool = mp.Pool(int(np.ceil(threads/2)))
    prepared_files = pool.imap_unordered(imap_prepare_masked_align, zip(PC_beds,
                                                                        total_region_beds,
                                                                        masked_genomes,
                                                                        repeat(input_bam),
                                                                        repeat(input_freads),
                                                                        repeat(input_rreads),
                                                                        repeat(mq_cutoff)))
    
    prepared_bams = []
    i=0
    for success, result, log_contents in prepared_files:
        i+=1
        print(f"\n************************************{i}_subprocess_start************************************\n", file=sys.stderr)
        if success:
            prepared_bams.append(result)
            print(f"Successfully prepared the masked bam file {result}. The log info are:\n{log_contents}\n", file=sys.stderr)
        else:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")
        print(f"\n************************************{i}_subprocess_end************************************\n", file=sys.stderr)

    prepared_files = prepared_bams
    pool.close()
    
    valid_files = [ f for f in prepared_files if not re.search(r".* failed$", f) ]
    valid_file_df = pd.read_csv(io.StringIO("\n".join(valid_files)), header=None, names=["PC_bed", "forward_read", "reverse_read", "masked_bam", "raw_bam", "whole_region_bed", "masked_genome"])
    valid_file_df["target_region_bed"] = target_region
    assert len(valid_file_df) > 0, "No valid files are generated from the input BAM file"

    logger.info("The metatable here is :\n{}\n".format(valid_file_df[:5].to_string(index=False)))
    
    # Now we merge the masked bam files into one to have a per-sample BAM file
    valid_remapped_bams = " ".join(valid_file_df["masked_bam"].drop_duplicates().dropna().tolist())
    remapped_bam = input_bam.replace(".bam", ".pooled.bam")
    remapped_bam_header = remapped_bam.replace(".bam", ".header")
    cmd = f"bash {bash_utils_hub} modify_bam_sq_lines {valid_remapped_bams.split()[0]} {ref_fasta} {remapped_bam_header}"
    executeCmd(cmd, logger=logger)
    
    test_cmd = f"bash {bash_utils_hub} check_bam_validity {remapped_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        logger.info(f"The merged pooled BAM file {remapped_bam} is not ready")
        execute = True
    else:
        if all([os.path.getmtime(remapped_bam) > os.path.getmtime(vb) for vb in valid_remapped_bams.split(" ")]):
            execute = False
        else:
            execute = True

    remapped_bam_list = remapped_bam.replace(".bam", ".bams.list.txt")
    with open(remapped_bam_list, "w") as f:
        f.write("\n".join(valid_remapped_bams.split(" ")))
    
    cmd = f"samtools merge -@ {threads} -h {remapped_bam_header} -b {remapped_bam_list} -o - | \
            samtools sort -O bam -o {remapped_bam} -@ {threads} && \
            samtools index {remapped_bam} && \
            ls -lht {remapped_bam} || \
            echo Failed to concatenate all the pooled BAM files. It wont be a fatal error but brings troubles to debugging and variant tracing."
    if execute:
        executeCmd(cmd)
    
    pool = mp.Pool(int(np.ceil(threads/3)))
    output_results = pool.imap_unordered(imap_call_polyploidy, zip( valid_file_df["forward_read"].to_list(), 
                                                                    valid_file_df["reverse_read"].to_list(), 
                                                                    valid_file_df["masked_bam"].to_list(),
                                                                    valid_file_df["raw_bam"].to_list(),
                                                                    valid_file_df["target_region_bed"].to_list(),
                                                                    valid_file_df["masked_genome"].to_list(),
                                                                    valid_file_df["PC_bed"].to_list() ))

    output_vcfs = []
    i=0
    for success, result, log_contents in output_results:
        i+=1
        print(f"\n************************************{i}_subprocess_start************************************\n", file=sys.stderr)
        if success:
            output_vcfs.append(result)
            print(f"Successfully call vcf from masked bam file: {result}. The log info are:\n{log_contents}\n", file=sys.stderr)
        else:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe full error log is :\n{log_contents}\n")
        print(f"\n************************************{i}_subprocess_end************************************\n", file=sys.stderr)
    pool.close()

    valid_vcfs = [v.split(",")[0] for v in output_vcfs if not re.search(r" failed$", v)]  # Filter out the failed samples
    logger.info("These VCFs are valid VCFs called from all the sub PC beds:\n{}\n".format("\n".join(valid_vcfs)))
    
    failed_regions = [r.split(" ")[0] for r in output_vcfs if re.search(r" failed$", r)]
    logger.info(f"These regions failed to generate a valid SDrecall VCF file for sample {input_bam}: \n{failed_regions}\n\n")


    execute = False
    valid_filtered_bams = " ".join(list(dict.fromkeys([v.split(",")[-2] for v in output_vcfs if not re.search(r" failed$", v)])))
    remapped_filtered_bam = remapped_bam.replace(".bam", ".filtered.bam")
    test_cmd = f"bash {bash_utils_hub} check_bam_validity {remapped_filtered_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        logger.info(f"The merged remapped filtered BAM file {remapped_filtered_bam} is not ready")
        execute = True
    else:
        if all([os.path.getmtime(remapped_filtered_bam) > os.path.getmtime(vb) for vb in valid_filtered_bams.split(" ")]):
            execute = False
        else:
            execute = True

    remapped_filtered_bam_list = remapped_filtered_bam.replace(".bam", ".bams.list.txt")
    with open(remapped_filtered_bam_list, "w") as f:
        f.write("\n".join(valid_filtered_bams.split(" ")))

    cmd = f"samtools merge -@ {threads} -h {remapped_bam_header} -b {remapped_filtered_bam_list} -o - | \
            samtools sort -O bam -o {remapped_filtered_bam} && \
            samtools index {remapped_filtered_bam} && \
            ls -lht {remapped_filtered_bam} || \
            echo Failed to concatenate all the filtered realigned BAM files. It wont be a fatal error but brings troubles to debugging and variant tracing."
    if execute:
        executeCmd(cmd)

    
    assert len(valid_vcfs) > 0
    # Then we need to use bcftools concat to merge the vcfs together
    tmp_list = prepare_tmp_file().name
    with open(tmp_list, "w") as f:
        for v in valid_vcfs: f.write(v + "\n")
        
    cmd = f"bash {bash_utils_hub} bcftools_concatvcfs -v {tmp_list} -o {output_vcf} && bash {bash_utils_hub} display_table {output_vcf} 20"
    executeCmd(cmd)


    # Below the code block is just for benchmarking and debugging purpose
    total_str = ""
    header = ""
    for v in valid_vcfs:
        meta_tab = v.replace(".vcf.gz", ".benchmark.tsv")
        if os.path.exists(meta_tab):
            with open(meta_tab, "r") as f:
                lines = f.readlines()
                header = lines[0]
                contents = "\n".join(lines[1:])
                total_str = total_str + contents + "\n"
    # Convert the total_str into a dataframe
    if len(total_str) > 0:
        total_str = total_str.replace("\n\n", "\n").replace("\n\n", "\n")
        total_df = pd.read_csv(io.StringIO(total_str), sep=",", header=None, names=header.split(","))
        total_df.to_csv(os.path.join(all_PC_folder, "sample_level_benchmark.tsv"), sep="\t", index=False)

    
    valid_vcf_str = " ".join(valid_vcfs)
    cmd = f"if bash {bash_utils_hub} check_vcf_validity {output_vcf}; then rm -f {valid_vcf_str}; fi"
    # executeCmd(cmd)
    
    # Then we need to annotate intrinsic variants in the output vcfs
    tmp_vcf = output_vcf.replace(".vcf.gz", ".tmp.vcf.gz")
    assert tmp_vcf != output_vcf
    executeCmd(f"mv {output_vcf} {tmp_vcf} && tabix -f -p vcf {tmp_vcf}")
    # annotate_intrinsic_vars(output_vcf, filtered_vcf, sd_working_dir=all_PC_folder)
    # logger.info("The merged VCF file with annotations on intrinsic variants for sample {} is {}".format(input_bam, filtered_vcf))
    
    # Then we need to tag out the variants that might be common according to in-house SD variants database
    annotate_inhouse_common_vars(tmp_vcf, output_vcf)
    executeCmd(f"bash {bash_utils_hub} display_table {output_vcf} 20")
    logger.info(f"The merged VCF annotated with potential in-house common variants is named after {output_vcf}")
    
    
    

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
    
    
    
    
    
    
    
    


