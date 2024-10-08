import os
import re
import shutil
import sys
from itertools import repeat
from multiprocessing import Pool
import logging

from build_regions_from_PC_clusters import establish_beds_per_PC_cluster
from utils import executeCmd

logger = logging.getLogger("SDrecall")

def imap_establish(tup_args):
    return establish_beds_per_PC_cluster(*tup_args)

def convert_nodes_into_hierachical_beds(grouped_qnode_cnodes: list,
                                        sd_paralog_pairs: dict,
                                        output_folder,
                                        ref_genome,
                                        target_region_bed = "",
                                        nthreads = 12,
                                        length = None,
                                        avg_frag_size = 400,
                                        std_frag_size = 140):
    # Establish a series of labels for PC-counterparts pairs
    # For disconnected query nodes cluster, we simply name them collectively as PC0
    # For connected query nodes cluster, we name them as PC1, PC2, PC3, ...
    import shutil
    from itertools import repeat
    import re
    
    logger.info("Start to put beds into the rest of PC-counterparts paired beds into subfolders")
    # We kind of need to restructure the grouped_qnode_cnodes to pass the information of fc-nfc pairs to the establish_beds_per_PC_cluster function
    new_results = []
    
    for result in grouped_qnode_cnodes:
        print(result)
        new_result = {"PCs": {}, "SD_counterparts": {}}
        for i in range(0, len(result["PCs"])):
            fc_node = result["PCs"][i]
            new_result['PCs'][i] = [fc_node]
            new_result['SD_counterparts'][i] = sd_paralog_pairs[fc_node]
        new_results.append(new_result)
        print(new_result)
        sys.exit()

    # Load balancing
    new_results = sorted(new_results, key = lambda x: sum(v[0][2] - v[0][1] for k,v in x["PCs"].items()) * sum(len(v) for k,v in x["SD_counterparts"].items()), reverse=True)
    labels = [ "PC" + str(n) for n in range(0, len(grouped_qnode_cnodes))]
    
    pool = Pool(nthreads)
    results = pool.imap_unordered(imap_establish, zip(new_results,
                                                      repeat(output_folder),
                                                      labels,
                                                      repeat(target_region_bed),
                                                      repeat(ref_genome),
                                                      repeat(length),
                                                      repeat(avg_frag_size),
                                                      repeat(std_frag_size)))
    i = 0
    intrinsic_bams = []
    for success, result, logs in results:
        print(f"\n\n*********************************** {i}_subprocess_start ***************************************", file = sys.stderr)
        if not success:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
        else:
            intrinsic_bams.append(result)
        print(logs, file = sys.stderr)
        print(f"*********************************** {i}_subprocess_end ***************************************\n\n", file = sys.stderr)
        i+=1
    
    pool.close()
    if not all([t[0] for t in results]):
        raise RuntimeError("Error happened during the parallel execution of establish_beds_per_PC_cluster. So stop the total python script here.")

    execute = False
    valid_intrinsic_bams = " ".join(intrinsic_bams)
    total_intrinsic_bam = os.path.join(output_folder, "total_intrinsic_alignments.bam")
    test_cmd = f"bash {bash_utils_hub} check_bam_validity {total_intrinsic_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        logger.info(f"The merged intrinsic alignment BAM file {total_intrinsic_bam} is not ready")
        execute = True
    else:
        if all([os.path.getmtime(total_intrinsic_bam) > os.path.getmtime(vb) for vb in intrinsic_bams]) and \
           (os.path.getmtime(total_intrinsic_bam) > os.path.getmtime(os.path.abspath(__file__))):
            execute = False
        else:
            execute = True

    intrinsic_bam_header = total_intrinsic_bam.replace(".bam", ".bam.header")
    cmd = f"bash {bash_utils_hub} modify_bam_sq_lines {intrinsic_bams[0]} {ref_genome} {intrinsic_bam_header}"
    executeCmd(cmd, logger=logger)

    intrinsic_bam_list = total_intrinsic_bam.replace(".bam", ".bams.list.txt")
    with open(intrinsic_bam_list, "w") as f:
        f.write("\n".join(intrinsic_bams))

    cmd = f"samtools merge -@ {nthreads} -h {intrinsic_bam_header} -b {intrinsic_bam_list} -o - | \
            samtools sort -O bam -o {total_intrinsic_bam} && \
            samtools index {total_intrinsic_bam} && \
            ls -lht {total_intrinsic_bam} || \
            echo Failed to concatenate all the filtered realigned BAM files. It wont be a fatal error but brings troubles to debugging and variant tracing."
    if execute:
        executeCmd(cmd)
    
    # Third remove leftover PC folders
    subdir_gen = os.walk(output_folder)
    # Prepare all target folder names in this time's generation, the naming syntax should follow the function called construct_folder
    if len(grouped_qnode_cnodes) > 1:
        target_folder_names = set([f"PC{x}_related_homo_regions" for x in range(0, len(grouped_qnode_cnodes))])
    else:
        target_folder_names = ["PC0_related_homo_regions"]
    logger.info(f"This time, we only have these PC regions generated: {target_folder_names}")
    first_level_dirs = next(subdir_gen)[1]  # The first item should be a 3-item tuple containin: 1. Current dir full path 2. All subdir names 3. All subfile names
    for subdir in first_level_dirs:
        if subdir not in target_folder_names:
            # Remove the subdir forcely
            shutil.rmtree(os.path.join(output_folder, subdir))
            
    # Fourth, extract all PC*_related_homo_regions.bed file and concat them together.
    first_level_dirs = next(os.walk(output_folder))[1]
    # Build a temp file to store the names of all sub files, if directly put all the names following cat, it might cause argument list too long error
    tmp_lst = prepare_tmp_file(suffix=".txt").name
    with open(tmp_lst, "w") as tf:
        total_lines = []
        for bed_file in [os.path.join(output_folder, subdir, subdir + ".bed") for subdir in first_level_dirs]:
            with open(bed_file, "r") as bf:
                total_lines = total_lines + bf.readlines()
        tf.write("\n".join([l.strip("\n") for l in total_lines if len(l) > 1]))

    executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total}".format(all = tmp_lst,
                                                                                                        total = os.path.join(output_folder, "all_PC_related_homo_regions.bed")))
    
    executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total} && rm {all}".format(all = tmp_lst,
                                                                                                        total = os.path.join(output_folder, "all_PC_regions.bed")))
    
    
    # Fifth, extract all intrinsic vcfs and use bcftools to concat them together
    intrinsic_vcfs = []
    for root, dirs, files in os.walk(output_folder):
        for file in files:
            if re.search(r'PC[0-9]+\.raw\.vcf\.gz$', file):
                intrinsic_vcfs.append(os.path.join(root, file))
    
    vcf_list_file = os.path.join(output_folder, "intrinsic_vcf.lst")
    final_intrinsic_vcf = os.path.join(output_folder, "all_pc_region_intrinsic_variants.vcf.gz")
    tmp_intrinsic_vcf = os.path.join(output_folder, "all_pc_region_intrinsic_variants.tmp.vcf.gz")
    with open(vcf_list_file, "w") as f:
        for v in intrinsic_vcfs: f.write(v + "\n")
    
    executeCmd(f"bcftools concat -o {tmp_intrinsic_vcf} -a --no-version -Oz -f {vcf_list_file} && bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp -o {final_intrinsic_vcf} -Oz {tmp_intrinsic_vcf} && ls -lh {final_intrinsic_vcf} && rm {tmp_intrinsic_vcf}")
    return