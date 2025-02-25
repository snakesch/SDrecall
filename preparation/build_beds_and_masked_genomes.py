import os
import sys
import logging

from pybedtools import BedTool
from pyfaidx import Fasta
import pysam

from src.utils import sortBed_and_merge, merge_bed_files
from src.const import *

from preparation.homoseq_region import HOMOSEQ_REGION
from preparation.masked_genome import create_masked_genome, get_reference_seq
from preparation.masked_realignment import create_minimap2_index, minimap2_align

import concurrent.futures

logger = logging.getLogger("SDrecall")

def process_pc_region(pc_region, label, outd, ref_genome, avg_insert_size, std_insert_size, threads):
        out_base = os.path.join(outd, label)
        os.makedirs(out_base, exist_ok=True)
        bed_out = os.path.join(out_base, label + ".bed")
        masked_genome = os.path.join(out_base, label + ".masked.fasta")
        output_bam = os.path.join(out_base, label + ".realigned.bam")

        ## Prepare for masked realignment
        create_pc_bed(pc_region, bed_out = bed_out, label = label)
        create_masked_genome(ref_genome, bed_out, masked_genome, avg_insert_size, std_insert_size)

        ## Realign using minimap2
        index_file = masked_genome + ".mmi"
        masked_contig_sizes = masked_genome + ".fai"

        ### Prepare reference sequences for realignment
        ref_sequences = os.path.join(out_base, label + ".masked.fastq")
        get_reference_seq(bed_out, ref_sequences, masked_genome, padding = avg_insert_size + std_insert_size)

        if not os.path.exists(index_file):  # Check if index already exists
             create_minimap2_index(masked_genome, index_file, threads) # Use ref_genome directly
        minimap2_align(ref_sequences, None, index_file, "asm20", label, output_bam, masked_contig_sizes, threads)
        return bed_out, output_bam

def get_region_and_realign(grouped_qnode_cnodes: list,
                                 sd_paralog_pairs: dict,
                                 outd,
                                 ref_genome,
                                 avg_insert_size = 400,
                                 std_insert_size = 140,
                                 threads = 1): 
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

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_pc_region, pc_region, label, outd, ref_genome,
                                    avg_insert_size, std_insert_size, threads): (pc_region, label)
                    for pc_region, label in zip(all_PC_regions, labels)}

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
    
    sortBed_and_merge(bed_out)

def merge_bam_files(bam_files: list[str], output_bam_path, header: pysam.AlignmentHeader) -> None:

    from tempfile import NamedTemporaryFile

    temp_merged_bam = NamedTemporaryFile(suffix=".bam", delete=False)
    temp_merged_bam_path = temp_merged_bam.name
    temp_merged_bam.close()  # Close it; pysam will open it again

    with pysam.AlignmentFile(temp_merged_bam_path, "wb", header=header) as outfile:
        for bam_file in bam_files:
            with pysam.AlignmentFile(bam_file, "rb") as infile:
                for read in infile:
                    outfile.write(read)
    pysam.sort("-n", "-o", output_bam_path, temp_merged_bam_path)
    pysam.index(output_bam_path)
    return