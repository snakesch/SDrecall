import os
import sys
from itertools import repeat
from multiprocessing import Pool
import logging
import uuid
import subprocess
from tempfile import NamedTemporaryFile

from pybedtools import BedTool

from src.utils import executeCmd, is_file_up_to_date, construct_folder_struc, perform_bedtools_sort_and_merge, update_plain_file_on_md5
from src.const import *
from src.log import error_handling_decorator
from preparation.homoseq_region import HOMOSEQ_REGION
from preparation.genome import Genome
from intrinsic_variants import getIntrinsicVcf

logger = logging.getLogger("SDrecall")

def convert_nodes_into_hierachical_beds(grouped_qnode_cnodes: list,
                                        sd_paralog_pairs: dict,
                                        output_folder,
                                        ref_genome,
                                        nthreads = 12,
                                        avg_frag_size = 400,
                                        std_frag_size = 140):
    # Label SD-paralog pairs. Name disconnected qnodes as PC0, connected qnodes as PC1, PC2, ...
    from shutil import rmtree
    from itertools import repeat
    import re
    
    # We need to restructure the grouped_qnode_cnodes to pass the information of sd-paralog pairs to the establish_beds_per_PC_cluster function
    new_results = []
    
    for result in grouped_qnode_cnodes:
        new_result = {"PCs": {}, "SD_counterparts": {}}
        for i in range(0, len(result["PCs"])):
            fc_node = result["PCs"][i]
            new_result['PCs'][i] = [fc_node]
            new_result['SD_counterparts'][i] = sd_paralog_pairs[fc_node]
        new_results.append(new_result)

    # Load balancing
    new_results = sorted(new_results, key = lambda x: sum(v[0][2] - v[0][1] for k,v in x["PCs"].items()) * sum(len(v) for k,v in x["SD_counterparts"].items()), reverse=True)
    labels = [ "PC" + str(n) for n in range(0, len(grouped_qnode_cnodes))]
    
    pool = Pool(nthreads)
    results = pool.imap_unordered(imap_establish, zip(new_results,
                                                      repeat(output_folder),
                                                      labels,
                                                      repeat(ref_genome),
                                                      repeat(avg_frag_size),
                                                      repeat(std_frag_size),
                                                      repeat(logger)))
    i = 0
    intrinsic_bams = []
    for success, result, logs in results:
        logger.debug(f"*********************************** {i}_subprocess_start ***************************************")
        if not success:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
        else:
            intrinsic_bams.append(result)
        logger.debug(logs)
        logger.debug(f"*********************************** {i}_subprocess_end ***************************************")
        i+=1
    
    pool.close()
    if not all([t[0] for t in results]):
        raise RuntimeError("Error occurred during the parallel execution of establish_beds_per_PC_cluster. Exit.")

    ## total_intrinsic_alignments.bam is a merged alignment file for BILC model.
    total_intrinsic_bam = os.path.join(output_folder, "total_intrinsic_alignments.bam")
    test_cmd = f"source {shell_utils}; check_bam_validity {total_intrinsic_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        ## Create total_intrinsic_alignments.bam
        intrinsic_bam_header = total_intrinsic_bam.replace(".bam", ".bam.header")
        cmd = f"source {shell_utils}; modify_bam_sq_lines {intrinsic_bams[0]} {ref_genome} {intrinsic_bam_header}"
        executeCmd(cmd, logger=logger)

        intrinsic_bam_list = total_intrinsic_bam.replace(".bam", ".bams.list.txt")
        with open(intrinsic_bam_list, "w") as f:
            f.write("\n".join(intrinsic_bams))

        cmd = f"samtools merge -@ {nthreads} -h {intrinsic_bam_header} -b {intrinsic_bam_list} -o - | \
                samtools sort -O bam -o {total_intrinsic_bam} && \
                samtools index {total_intrinsic_bam} && \
                ls -lht {total_intrinsic_bam} || \
                echo Failed to concatenate all the filtered realigned BAM files."
        executeCmd(cmd)
    
    # Third remove leftover PC folders
    subdir_gen = os.walk(output_folder)
    # Prepare all target folder names in this time's generation, the naming syntax should follow the function called construct_folder
    if len(grouped_qnode_cnodes) > 1:
        target_folder_names = set([f"PC{x}_related_homo_regions" for x in range(0, len(grouped_qnode_cnodes))])
    else:
        target_folder_names = ["PC0_related_homo_regions"]
    logger.info(f"Generated PC regions: {sorted(target_folder_names)}")
    first_level_dirs = next(subdir_gen)[1]  # The first item should be a 3-item tuple containin: 1. Current dir full path 2. All subdir names 3. All subfile names
    for subdir in first_level_dirs:
        if subdir not in target_folder_names:
            # Remove the subdir forcely
            rmtree(os.path.join(output_folder, subdir))
            
    # Fourth, extract all PC*_related_homo_regions.bed file and concat them together.
    first_level_dirs = next(os.walk(output_folder))[1]

    with NamedTemporaryFile(dir = "/tmp") as fp:
        total_lines = []
        for bed_file in [os.path.join(output_folder, subdir, subdir + ".bed") for subdir in first_level_dirs]:
            with open(bed_file, "r") as bf:
                total_lines = total_lines + bf.readlines()
        encoded_string = "".join(total_lines).rstrip("\n").encode()
        fp.write(encoded_string)

        executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total}".format(all = fp.name,
                                                                                                            total = os.path.join(output_folder, "all_PC_related_homo_regions.bed")))
        
        executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total}".format(all = fp.name,
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
    
    executeCmd(f"bcftools concat -o {tmp_intrinsic_vcf} -a --no-version -Oz -f {vcf_list_file} && bcftools sort -o {final_intrinsic_vcf} -Oz {tmp_intrinsic_vcf} && ls -lh {final_intrinsic_vcf} && rm {tmp_intrinsic_vcf}")
    return

def imap_establish(tup_args):
    return establish_beds_per_PC_cluster(*tup_args)

@error_handling_decorator
def establish_beds_per_PC_cluster(cluster_dict={"PCs":{},
                                                "SD_counterparts":{}},
                                  base_folder = "",
                                  label = "PC0",
                                  ref_genome = "",
                                  avg_frag_size = 400,
                                  std_frag_size = 140,
                                  threads = 2,
                                  logger = logger):

    """
    input cluster dict now looks like this:
    {
        "PCs": {0: [], 1: [], 2: []}, 
        "SD_counterparts": {0: [], 1: [], 2: []}
    }
    """

    ## Initialize file structure
    paths = construct_folder_struc(base_folder=base_folder, label=label, logger=logger)

    # First convert the disconnected nodes to beds, each node is a 3-tuple (chr, start, end)
    tmp_id = str(uuid.uuid4())
    tmp_pc_bed = paths["PC_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_pc_bed = paths["PC_bed"].replace(".bed", ".raw.bed")
    tmp_counterparts_bed = paths["Counterparts_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_counterparts_bed = paths["Counterparts_bed"].replace(".bed", ".raw.bed")
    tmp_total_bed = paths["All_region_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_total_bed = paths["All_region_bed"].replace(".bed", ".raw.bed")
    
    with open(tmp_pc_bed, "w") as f:
        for idx, records in cluster_dict["PCs"].items():
            for record in records:
                if len(record) >= 3:
                    # Invoke __iter__ method instead of __getitem__
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3]]) + "\n")
                elif len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3]]) + "\n")
    
    subprocess.run(f"cp -f {tmp_pc_bed} {raw_pc_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_pc_bed, logger=logger)
    update_plain_file_on_md5(paths["PC_bed"], tmp_pc_bed, logger=logger)
            
    ## Create the counterparts region bed file
    with open(tmp_counterparts_bed, "w") as f:
        for idx, records in cluster_dict["SD_counterparts"].items():
            fc_node = cluster_dict["PCs"][idx][0]
            for record in records:
                assert isinstance(record, HOMOSEQ_REGION), "The record in the SD_counterparts list is not a HOMOSEQ_REGION object: {}".format(record)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)     
                ## Invoke __iter__ method instead of __getitem__
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3]]) + "\n")
    subprocess.run(f"cp -f {tmp_counterparts_bed} {raw_counterparts_bed}", shell=True)

    ## Remove PC regions from counterpart BEDs and force strandedness
    BedTool(tmp_counterparts_bed).subtract(BedTool(paths["PC_bed"]), s=True).saveas(tmp_counterparts_bed)
    perform_bedtools_sort_and_merge(tmp_counterparts_bed, logger=logger)
    update_plain_file_on_md5(paths["Counterparts_bed"], tmp_counterparts_bed, logger=logger)
    
    ## Create the total region bed file by combining PC BED and counterpart BED
    ## Invoke __iter__ method instead of __getitem__
    with open(tmp_total_bed, "w") as f:
        for idx, records in cluster_dict["PCs"].items():
            for record in records:
                if len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3], f"FC:{label}_{idx}"]) + "\n")
                elif len(record) >= 3:
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3], f"FC:{label}_{idx}"]) + "\n")
        for idx, records in cluster_dict["SD_counterparts"].items():
            fc_node = cluster_dict["PCs"][idx][0]
            for record in records:
                if not isinstance(record, HOMOSEQ_REGION):
                    logger.critical("Some elements of the SD_counterparts list are not HOMOSEQ_REGION objects: {}".format(record))
                    sys.exit(1)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)
                # The rela_start and rela_end are based on the corresponding PC region
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3], f"NFC:{label}_{idx}"]) + "\n")
    
    subprocess.run(f"cp -f {tmp_total_bed} {raw_total_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_total_bed, logger = logger)
    update_plain_file_on_md5(paths["All_region_bed"], tmp_total_bed, logger=logger)
    
    contig_sizes = ref_genome.replace(".fasta", ".fasta.fai")
    # Prepare masked genomes
    masked_genome = Genome(ref_genome).mask(paths["PC_bed"], avg_frag_size = avg_frag_size, std_frag_size=std_frag_size, genome=contig_sizes, logger=logger)
      
    # Call intrinsic variants
    bam_path = getIntrinsicVcf( pc_bed = paths["PC_bed"], 
                                all_homo_regions_bed = paths["All_region_bed"], 
                                pc_masked = masked_genome,
                                ref_genome = ref_genome,
                                avg_frag_size = avg_frag_size,
                                std_frag_size = std_frag_size,
                                threads = threads)
    return bam_path