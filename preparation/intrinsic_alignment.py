import os
import logging

from preparation.seq import getRawseq
from src.utils import executeCmd, is_file_up_to_date, prepare_tmp_file
from src.const import shell_utils
from src.log import logger


def getIntrinsicBam(rg_bed, 
                    all_homo_regions_bed, 
                    rg_masked,
                    ref_genome, 
                    intrinsic_bam,
                    rg_label = None,
                    avg_frag_size = 400,
                    std_frag_size = 120,
                    threads = 2):
    '''
    Generate an intrinsic VCF file from given genomic data.

    This function takes a masked genome and associated BED files to perform variant calling.
    It maps the counterpart reference sequences against the masked genome and generates a VCF 
    file containing the identified variants.

    Parameters:
    - rg_bed (str): Path to the BED file of the PC region.
    - all_homo_regions_bed (str): Path to the BED file of all homologous regions (SD counterparts) related to the PC region.
    - rg_masked (str): Path to the masked genome in FASTA format.
    - ref_genome (str): Path to the reference genome used for mapping.
    - avg_frag_size (int): Average fragment size of NGS reads. (400)
    - std_frag_size (int): Standard deviation of fragment size of NGS reads. (120)
    - threads (int): Number of threads to use for processing. (2)

    Returns:
    - str: Path to the generated BAM file.
    '''

    intrinsic_fastq = intrinsic_bam.replace(".bam", ".fastq")
    rg_label = rg_label or os.path.basename(rg_bed.replace(".bed",""))
        
    # Reference sequences of SD counterparts are extracted from reference genome
    getRawseq(all_homo_regions_bed, intrinsic_fastq, ref_genome, padding = avg_frag_size + std_frag_size)
    logger.info(f"Reference sequences for intrinsic alignment from {rg_bed} are written to: {intrinsic_fastq}")
    
    # Retrieved counterpart sequences are mapped against masked genomes using minimap2
    if not os.path.exists(intrinsic_bam) or not is_file_up_to_date(intrinsic_bam, [rg_masked, shell_utils, os.path.abspath(__file__)]):
        cmd = f"bash {shell_utils} \
                independent_minimap2_masked \
                -f {intrinsic_fastq} \
                -a {rg_masked} \
                -o {intrinsic_bam} \
                -s all_PC \
                -t {threads} \
                -i {rg_label} \
                -m asm20 \
                -g {ref_genome}"
        executeCmd(cmd, logger=logger)
        intrinsic_bam = filter_intrinsic_alignments(intrinsic_bam, logger=logger)

    return intrinsic_bam



def filter_intrinsic_alignments(bam_file, output_file=None, logger = logger):
    """
    Filter out alignments where the start position in the read name
    matches the alignment position in the reference. Meaning the reference sequence is mapped to its original genomic location
    
    Parameters:
    -----------
    bam_file : str
        Path to input BAM file
    output_file : str
        Path to output filtered BAM file
    """
    import pysam
    import re
    
    qname_regex = re.compile(r'(.*):(\d+)-(\d+)(.*)')
    tmp_output = prepare_tmp_file(suffix=".bam").name
    
    with pysam.AlignmentFile(bam_file, "rb") as input_bam:
        with pysam.AlignmentFile(tmp_output, "wb", header=input_bam.header) as output_bam:
            total_reads = 0
            filtered_reads = 0
            
            for read in input_bam:
                total_reads += 1
                
                # Extract start position from qname (format: chr:start-end)
                match = qname_regex.match(read.query_name)
                if match:
                    qname_start = int(match.group(2))
                    
                    # Compare with alignment position
                    if qname_start != read.reference_start + 1:  # pysam.read.reference_start is 0-based
                        output_bam.write(read)
                    else:
                        filtered_reads += 1
                else:
                    # If qname doesn't match the expected format, keep the read
                    output_bam.write(read)
            
            logger.info(f"Filtered out {filtered_reads} of {total_reads} reads where reference position matches qname start from {bam_file}")
    
    # Index the output BAM
    if output_file is None:
        cmd = f"mv {tmp_output} {bam_file} && samtools index {bam_file}"
        executeCmd(cmd, logger=logger)
        return bam_file
    else:
        cmd = f"mv {tmp_output} {output_file} && samtools index {output_file}"
        executeCmd(cmd, logger=logger)
        return output_file