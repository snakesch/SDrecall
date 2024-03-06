#!/usr/bin/env python3

import logging
logger = logging.getLogger("root")

def init_dtree(base_folder, label):
    import os
    
    parent_folder           = os.path.join(base_folder, label + "_related_homo_regions")
    total_bed_path          = os.path.join(parent_folder, label + "_related_homo_regions.bed")
    counterparts_bed_path   = os.path.join(parent_folder, label + "_counterparts_regions.bed")
    PC_folder               = os.path.join(parent_folder, label)
    PC_bed_path             = os.path.join(PC_folder, label + ".bed")

    folder_catalog = {"base_folder": parent_folder, "PC_bed": PC_bed_path, "All_region_bed": total_bed_path, "Counterparts_bed": counterparts_bed_path}

    for d in [parent_folder, PC_folder]:
        os.makedirs(d, exist_ok=True)

    return folder_catalog

def executeCmd(cmd) -> None:
    
    import subprocess
    
    proc = subprocess.run(cmd, shell=True)
    code = proc.returncode
    cmd_list = cmd.split(" ")
    if code != 0:
        message = "Error in {}:\n{}\n".format(" ".join(cmd_list), proc.stdout.decode()) if cmd_list[1][0] != "-" else "Error in {}:\n{}\n".format(cmd_list, proc.stdout.decode())
        raise RuntimeError(message)
    
    