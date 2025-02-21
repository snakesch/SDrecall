import os
import subprocess
import logging
import tempfile
import re
from typing import List

from pybedtools import BedTool

logger = logging.getLogger("SDrecall")

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

def update_plain_file_on_md5(old_file, new_file, logger=logging.getLogger('SDrecall')):
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

# - Construct file tree for homologous region BEDs - #

def construct_folder_struc(base_folder,
                           label="",
                           logger = logging.getLogger('SDrecall')):
    parent_folder_name = label+ "_related_homo_regions"
    parent_folder_full_path = os.path.join(base_folder, parent_folder_name)

    os.makedirs(parent_folder_full_path, exist_ok=True)

    total_bed_name = label + "_related_homo_regions.bed"
    total_bed_path = os.path.join(parent_folder_full_path, total_bed_name)

    counterparts_bed_name = label + "_counterparts_regions.bed"
    counterparts_bed_path = os.path.join(parent_folder_full_path, counterparts_bed_name)

    PC_folder_name = label
    PC_folder_full_path = os.path.join(parent_folder_full_path, PC_folder_name)

    os.makedirs(PC_folder_full_path, exist_ok=True)

    PC_bed_name = label + ".bed"
    PC_bed_path = os.path.join(PC_folder_full_path, PC_bed_name)

    return {"base_folder_path": parent_folder_full_path,
            "PC_bed": PC_bed_path,
            "All_region_bed": total_bed_path,
            "Counterparts_bed": counterparts_bed_path}

# - BED file manipulation using pybedtools - #

def sortBed_and_merge(bed_file, output=None, logger = logger):
    
    # Load your BED file
    bed = BedTool(bed_file)
    merged_bed = bed.sort().merge(s=True, c='4,5,6', o='first,first,first')
    
    # Inplace file sort and merge
    if not output:
        output = bed_file

    # You can save the results to a new file
    merged_bed.saveas(output)
    
    return None
    
def merge_bed_files(bed_files: List[str], logger = logger) -> BedTool:
    """Merges and sorts multiple BED files using pybedtools. Returns BedTool."""
    bed_files = list(dict.fromkeys(bed_files))  # Remove duplicates
    if not bed_files:
        return BedTool("", from_string=True)  # Return empty BedTool

    bedtools_list = [BedTool(f) for f in bed_files]
    merged_bedtool = bedtools_list[0]
    for bt in bedtools_list[1:]:
        merged_bedtool = merged_bedtool.cat(bt, postmerge=False)
    return merged_bedtool.sort().merge()

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


