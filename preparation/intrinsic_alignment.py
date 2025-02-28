import os
import logging

from preparation.seq import getRawseq
from src.utils import executeCmd, is_file_up_to_date
from src.const import shell_utils

logger = logging.getLogger("SDrecall")

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

    fastq_dir = os.path.dirname(all_homo_regions_bed)
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

    return intrinsic_bam