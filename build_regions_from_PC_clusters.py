import uuid
import logging
import subprocess

from pybedtools import BedTool

from log import error_handling_decorator
from homoseq_region import HOMOSEQ_REGION
from genome import Genome
from intrinsic_variants import getIntrinsicVcf
from utils import construct_folder_struc, perform_bedtools_sort_and_merge, update_plain_file_on_md5

logger = logging.getLogger("SDrecall")

@error_handling_decorator
def establish_beds_per_PC_cluster(cluster_dict={"PCs":{},
                                                "SD_counterparts":{}},
                                  base_folder = "/paedyl01/disk1/yangyxt/indexed_genome/SD_priority_component_pairs",
                                  label = "PC0",
                                  target_region_bed = "",
                                  ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                  length = None,
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

    logger.info("The input cluster_dict is:\n{}\n".format(cluster_dict))
    assert isinstance(cluster_dict["SD_counterparts"][0][0], HOMOSEQ_REGION)
    assert isinstance(cluster_dict["PCs"][0][0], tuple), f"The first element in the PC cluster should be a tuple, but it is {type(cluster_dict['PCs'][0][0])} and it looks like {cluster_dict['PCs'][0][0]}"
    paths = construct_folder_struc(base_folder=base_folder, label=label, logger=logger)
    # logger.info("The PC bed file is {}, the counterparts bed file is {}, the total bed file is {}".format(paths["PC_bed"], paths["Counterparts_bed"], paths["All_region_bed"]))
    
    # First convert the disconnected nodes to beds, each node is a tuple consisting of three values (chr, start, end)
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
                    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3]]) + "\n")
                elif len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3]]) + "\n")
    
    subprocess.run(f"cp -f {tmp_pc_bed} {raw_pc_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_pc_bed, logger=logger)
    update_plain_file_on_md5(paths["PC_bed"], tmp_pc_bed, logger=logger)
            
    # Then compose the counterparts region bed file
    with open(tmp_counterparts_bed, "w") as f:
        for idx, records in cluster_dict["SD_counterparts"].items():
            fc_node = cluster_dict["PCs"][idx][0]
            for record in records:
                assert isinstance(record, HOMOSEQ_REGION), "The record in the SD_counterparts list is not a HOMOSEQ_REGION object: {}".format(record)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)     
                # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3]]) + "\n")
    subprocess.run(f"cp -f {tmp_counterparts_bed} {raw_counterparts_bed}", shell=True)

    # Now we need to subtract the PC from counterpart beds incase they have overlaps,  force strandness
    BedTool(tmp_counterparts_bed).subtract(BedTool(paths["PC_bed"]), s=True).saveas(tmp_counterparts_bed)
    perform_bedtools_sort_and_merge(tmp_counterparts_bed, logger=logger)
    
    update_plain_file_on_md5(paths["Counterparts_bed"], tmp_counterparts_bed, logger=logger)
    
    # Then concatenate the two bed file together to generate the total region bed file
    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
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
                assert isinstance(record, HOMOSEQ_REGION), "The record in the SD_counterparts list is not a HOMOSEQ_REGION object: {}".format(record)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)
                # The rela_start and rela_end are based on the corresponding FC interval
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3], f"NFC:{label}_{idx}"]) + "\n")
    
    subprocess.run(f"cp -f {tmp_total_bed} {raw_total_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_total_bed, logger = logger)
    update_plain_file_on_md5(paths["All_region_bed"], tmp_total_bed, logger=logger)
    
    contig_sizes = ref_genome.replace(".fasta", ".fasta.fai")
    # Alongside preparing the bed files, we can also perform masked genome preparation
    masked_genome, updated = Genome(ref_genome).mask(paths["PC_bed"], avg_frag_size = avg_frag_size, std_frag_size=std_frag_size, genome=contig_sizes, logger=logger)
    if updated:
        assert os.path.getmtime(masked_genome) > os.path.getmtime(paths["PC_bed"])
    
    # And we need to prepare the intrinsic VCF based on the masked genome and generated BED files
    bam_path = getIntrinsicVcf( pc_bed = paths["PC_bed"], 
                                all_homo_regions_bed = paths["All_region_bed"], 
                                counter_bed = paths["Counterparts_bed"],
                                pc_masked = masked_genome,
                                ref_genome = ref_genome,
                                avg_frag_size = avg_frag_size,
                                std_frag_size = std_frag_size,
                                threads = threads, 
                                logger = logger)
    return bam_path