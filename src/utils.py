import os
import subprocess
import logging
import tempfile
import re
import numpy as np
from typing import List, Tuple

from pybedtools import BedTool
from src.log import logger
from src.const import shell_utils


# - Miscellaneous helper functions - #
def executeCmd(cmd, logger = logger) -> None:

    proc = subprocess.run(cmd, shell=True, capture_output=True)
        
    logger.debug(f"Command run:{cmd}")
    logger.debug(f"Return code: {proc.returncode}")

    err_msg = proc.stderr.decode()
    if err_msg != "":
        logger.debug(proc.stderr.decode())
        
    if proc.returncode != 0:
        raise RuntimeError(err_msg)

    return proc.stdout.decode()


def prepare_tmp_file(tmp_dir="/tmp", **kwargs):
    os.makedirs(tmp_dir, exist_ok=True)
    return tempfile.NamedTemporaryFile(dir = tmp_dir, delete = False, **kwargs)


def is_file_up_to_date(file_to_check, list_of_dependency_files):
    file_to_check_time = os.path.getmtime(file_to_check)
    return all(file_to_check_time > os.path.getmtime(dep_file) for dep_file in list_of_dependency_files)


def update_plain_file_on_md5(old_file, new_file, logger=logger):
    import hashlib

    old_md5 = hashlib.md5(open(old_file,'r').read().encode()).hexdigest() if os.path.exists(old_file) else "non_existent"

    if not os.path.exists(new_file):
        raise FileNotFoundError("File {} does not exist.".format(new_file))

    if os.stat(new_file).st_size == 0:
        raise FileExistsError("File {} is empty. Not using it to update the original file {}".format(new_file, old_file))

    new_md5 = hashlib.md5(open(new_file,'r').read().encode()).hexdigest()

    if new_md5 == old_md5:
        logger.warning("The new file {} is identical to the old one {} so no updates will be carried out. Deleting new file {}".format(new_file, old_file, new_file))
        os.remove(new_file)
        return False
    else:
        import shutil
        logger.debug("The new file {} is different from the old one {} so the old one will be replaced with the new one.".format(new_file, old_file))
        shutil.move(new_file, old_file)
        executeCmd(f"ls -lht {old_file}", logger=logger)
        return True


# - BED file manipulation using pybedtools - #
def sortBed_and_merge(bed_file, output=None):
    
    # Load your BED file
    bed = BedTool(bed_file)
    num_cols = bed.field_count()
    c_arg = ','.join(str(i) for i in range(4, num_cols + 1))  # '4,5,6,...'
    if num_cols <= 6:
        o_arg = ','.join('first' for _ in range(4, num_cols + 1))  # 'first,first,first,...'
    elif num_cols == 7:
        o_arg = 'distinct,distinct,distinct,distinct'
   
    merged_bed = bed.sort().merge(s=True, c=c_arg, o=o_arg)
    
    # Inplace file sort and merge
    if output is None:
        output = bed_file

    # You can save the results to a new file
    merged_bed.saveas(output)
    
    return None
    
def merge_bed_files(bed_files: List[str]) -> BedTool:
    """Merges multiple BED files using pybedtools. Returns BedTool."""
    bed_files = list(dict.fromkeys(bed_files))  # Remove duplicates
    if not bed_files:
        return BedTool("", from_string=True)  # Return empty BedTool

    bedtools_list = [BedTool(f) for f in bed_files]
    merged_bedtool = bedtools_list[0]
    for bt in bedtools_list[1:]:
        merged_bedtool = merged_bedtool.cat(bt, postmerge=False)
    return merged_bedtool.sort()


def filter_bed_by_interval_size(bed_obj, interval_size_cutoff):
    return bed_obj.filter(lambda x: len(x) > interval_size_cutoff).sort()


# - VCF/BCF file manipulation using bcftools -# 
def combine_vcfs(*vcfs, output=None, threads=4):
    ## Given a list of VCFs, runs bcftools concat + tabix + sort, returns the processed VCF path
    import subprocess
    
    if output is None:
        raise ValueError("No output specified. ")
    elif not output.endswith(".gz"):
        output += ".gz"
        
    cmd = f"bcftools concat -a --threads {threads} -d exact "
    
    for vcf in vcfs:
        cmd += f"{vcf} "
        
    cmd += f"2>/dev/null | bcftools sort /dev/stdin -o {output} -Oz 2>/dev/null && tabix -f -p vcf {output} 2>/dev/null"
    
    proc = subprocess.run(cmd, shell=True, capture_output=False)
    
    if proc.returncode != 0: raise RuntimeError("Error combining VCFs")


def configure_parallelism(total_threads: int,
                              threads_per_job: int = 4) -> Tuple[int, int]:
    """
    Configures parallelisation for a given number of total threads and threads per job.
    
    Args:
        total_threads (int): The total number of threads available.
    """
    num_jobs = np.ceil(total_threads / threads_per_job)
    return num_jobs, threads_per_job


def merge_bams(bam_list: list, 
               merged_bam: str, 
               ref_fasta = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
               threads = 2,
               logger = logger):
    # For data visualization and debugging
    merged_bam_header = merged_bam.replace(".bam", ".header")
    cmd = f"bash {shell_utils} modify_bam_sq_lines {bam_list[0]} {ref_fasta} {merged_bam_header}"
    executeCmd(cmd, logger=logger)
    
    test_cmd = f"bash {shell_utils} quick_check_bam_validity {merged_bam}"
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






