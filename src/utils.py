import os
import subprocess
import logging
import uuid

from pybedtools import BedTool

def executeCmd(cmd, logger = logging.getLogger('SDrecall')) -> None:
    
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    logger.debug(f"Command run:\n{cmd}\nReturn code: {proc.returncode}\nError log:\n{proc.stderr.decode()}")
    code = proc.returncode
    
    if code != 0:
        cmd_lst = cmd.split(" ")
        if cmd_lst[1][0] != "-":
            raise RuntimeError("Error in " + " ".join(cmd_lst))
        else:
            raise RuntimeError("Error in {}".format(cmd_lst))
    
    return proc.stdout.decode()

def is_file_up_to_date(file_to_check, list_of_dependency_files):
    file_to_check_time = os.path.getmtime(file_to_check)
    return all(file_to_check_time > os.path.getmtime(dep_file) for dep_file in list_of_dependency_files)

def update_plain_file_on_md5(old_file, new_file, logger=logging.getLogger('SDrecall')):
    import hashlib
    
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
        logger.debug("The new file {} is different from the old one {} so the old one will be replaced with the new one.".format(new_file, old_file))
        shutil.move(new_file, old_file)
        executeCmd(f"ls -lht {old_file}", logger=logger)
        return True

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

def string_to_tuple(string):
    """Convert a string representation of a tuple into an actual tuple."""
    import ast
    try:
        return ast.literal_eval(string)
    except (ValueError, SyntaxError) as e:
        raise ValueError(f"Failed to convert string to tuple: {string}") from e

def perform_bedtools_sort_and_merge(bed_file,
                                    output_bed_file=None,
                                    logger = logging.getLogger('SDrecall')):
    # Load your BED file
    bed = BedTool(bed_file)

    # Sort the BED file
    sorted_bed = bed.sort()

    # Merge the BED file
    merged_bed = sorted_bed.merge(s=True, c='4,5,6', o='first,first,first')

    if not output_bed_file:
        output_bed_file = bed_file
        
    # You can save the results to a new file
    merged_bed.saveas(output_bed_file)

def filter_bed_by_interval_size(bed_obj, interval_size_cutoff):
    
    return bed_obj.filter(lambda x: len(x) > interval_size_cutoff).sort()

