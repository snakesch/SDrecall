import os
import sys
import logging

from pybedtools import BedTool
from pyfaidx import Fasta
import pysam

from src.utils import executeCmd, construct_folder_struc, sortBed_and_merge, combine_vcfs, merge_bed_files, update_plain_file_on_md5
from src.const import shell_utils
from src.log import error_handling_decorator, logger
from preparation.homoseq_region import HOMOSEQ_REGION
from preparation.masked_genome import create_masked_genome, get_reference_seq
from preparation.masked_realignment import create_minimap2_index, minimap2_align

import concurrent.futures


def build_beds_and_masked_genomes(grouped_qnode_cnodes: list,
                                  sd_paralog_pairs: dict,
                                  output_folder,
                                  ref_genome,
                                  nthreads = 12,
                                  avg_frag_size = 400,
                                  std_frag_size = 140):
    # Label SD-paralog pairs. Name disconnected qnodes as PC0, connected qnodes as PC1, PC2, ...
    all_PC_regions = []
    
    for r in grouped_qnode_cnodes:
        new_r = {"PCs": {}, "SD_counterparts": {}}
        for i in range(0, len(r["PCs"])):
            fc_node = r["PCs"][i]
            new_r['PCs'][i] = [fc_node]
            new_r['SD_counterparts'][i] = sd_paralog_pairs[fc_node]
        all_PC_regions.append(new_r)

    # Load balancing
    all_PC_regions = sorted(all_PC_regions, key = lambda x: sum(v[0][2] - v[0][1] for k,v in x["PCs"].items()) * sum(len(v) for k,v in x["SD_counterparts"].items()), reverse=True)
    labels = [ "PC" + str(n) for n in range(0, len(grouped_qnode_cnodes))]

    # output BAMs generated here are intrinsic BAM files in the manuscript
    output_bams = []
    output_beds = []

    ## total_intrinsic_alignments.bam is a merged alignment file for BILC model.
    total_intrinsic_bam = os.path.join(output_folder, "total_intrinsic_alignments.bam")
    
    ## Create total_intrinsic_alignments.bam
    intrinsic_bam_header = total_intrinsic_bam.replace(".bam", ".bam.header")
    cmd = f"source {shell_utils} && modify_bam_sq_lines {intrinsic_bams[0]} {ref_genome} {intrinsic_bam_header}"
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
               
    # Fourth, extract all PC*_related_homo_regions.bed file and concat them together.
    beds = glob(os.path.join(output_folder, "*/*_all/", "*_related_homo_regions.bed"))
    combined_bed = merge_bed_files(beds)
    combined_bed.saveas(os.path.join(output_folder, "all_PC_related_homo_regions.bed"))
    # combined_bed.saveas(os.path.join(output_folder, "all_PC_regions.bed"))
    
    # Fifth, extract all intrinsic vcfs and use bcftools to concat them together
    intrinsic_vcfs = []
    for root, dirs, files in os.walk(output_folder):
        for file in files:
            if re.search(r'PC[0-9]+\.raw\.vcf\.gz$', file):
                intrinsic_vcfs.append(os.path.join(root, file))
    
    final_intrinsic_vcf = os.path.join(output_folder, "all_pc_region_intrinsic_variants.vcf.gz")
    combine_vcfs(*intrinsic_vcfs, output=final_intrinsic_vcf)

    for future in concurrent.futures.as_completed(futures):
        try:
            bed_out, output_bam = future.result()
            output_beds.append(bed_out)
            output_bams.append(output_bam)
        except Exception as e:
            logger.error(f"Error processing a region: {e}")

    ## Merge all PC BED files
    all_PC_regions_out = os.path.join(outd, "all_PC.bed")
    merge_bed_files(output_beds).saveas(all_PC_regions_out)

    ## Reheader using the first BAM and merge all output BAMs with the new header
    merged_bam_out = os.path.join(outd, "total_intrinsic.bam")
    header = create_bam_header_from_fasta(output_bams[0], ref_genome)
    merge_bam_files(output_bams, merged_bam_out, header)

    return

def create_bam_header_from_fasta(input_bam: str, ref_fasta: str) -> pysam.AlignmentHeader:
    """
    Creates a new BAM header based on an input BAM file's header, replacing
    the @SQ lines with information derived from a reference FASTA file.

    Args:
        input_bam: Path to the input BAM file.
        ref_fasta: Path to the reference FASTA file.

    Returns:
        A pysam.AlignmentHeader object with updated @SQ lines.
    """

    fasta = Fasta(ref_fasta)

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        header_dict = bamfile.header.to_dict()

    header_dict["SQ"] = []

    sq_lines = []
    for name in fasta.keys():
        length = len(fasta[name])
        sq_lines.append({"SN": name, "LN": length})
    header_dict["SQ"] = sq_lines # type: ignore

    new_header = pysam.AlignmentHeader.from_dict(header_dict)  # type: ignore
    return new_header

def create_pc_bed(cluster_dict, bed_out = "", label = "PC0") -> None:
    
    # Tabulate query nodes and counterpart nodes in a BED for each PC region
    with open(bed_out, "w") as f:
        for idx, records in cluster_dict["PCs"].items():
            for record in records: ## use __iter__ here
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
    
    sortBed_and_merge(paths["All_region_bed"])
    
    contig_sizes = ref_genome.replace(".fasta", ".fasta.fai")
    
    # Prepare masked genomes
    masked_genome_path = os.path.join( os.path.dirname(paths["PC_bed"]), label + ".masked.fasta")
    masked_genome = Genome(ref_genome).mask(paths["PC_bed"], avg_frag_size = avg_frag_size, std_frag_size=std_frag_size, genome=contig_sizes, logger=logger, path=masked_genome_path)
    
    # Call intrinsic variants
    bam_path = getIntrinsicVcf( pc_bed = paths["PC_bed"], 
                                all_homo_regions_bed = paths["All_region_bed"], 
                                pc_masked = masked_genome,
                                ref_genome = ref_genome,
                                avg_frag_size = avg_frag_size,
                                std_frag_size = std_frag_size,
                                threads = threads)
    return bam_path
