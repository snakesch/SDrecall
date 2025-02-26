import os
import logging

from preparation.seq import getRawseq
from src.utils import executeCmd, is_file_up_to_date
from src.const import *

logger = logging.getLogger("SDrecall")

def getIntrinsicBam(rg_bed, 
                    all_homo_regions_bed, 
                    rg_masked,
                    ref_genome, 
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
    vcf_dir = os.path.dirname(rg_bed)
    raw_fq_path = os.path.join(fastq_dir, os.path.basename(all_homo_regions_bed)[:-3] + "raw.fastq")
    bam_path = os.path.join(vcf_dir, os.path.basename(rg_bed)[:-3] + "raw.bam")
    vcf_path = bam_path.replace(".bam", ".vcf.gz")
    rg_label = os.path.basename(rg_bed.replace(".bed",""))
        
    # Reference sequences of SD counterparts are extracted from reference genome
    raw_fq_path = getRawseq(all_homo_regions_bed, raw_fq_path, ref_genome, padding = avg_frag_size + std_frag_size)
    logger.info(f"Reference sequences for calling intrinsic variants from {rg_bed} are written to: {raw_fq_path}")
    
    # Retrieved counterpart sequences are mapped against masked genomes using minimap2
    if not os.path.exists(bam_path) or not is_file_up_to_date(bam_path, [rg_masked, shell_utils, os.path.abspath(__file__)]):
        rg_masked_index = rg_masked.replace(".fasta", ".mmi")
        cmd = f"source {shell_utils} && independent_minimap2_masked \
                -f {raw_fq_path} \
                -a {rg_masked} \
                -o {bam_path} \
                -s all_PC \
                -t {threads} \
                -i {rg_label} \
                -m asm20 \
                -g {ref_genome}"
        executeCmd(cmd, logger=logger)

    ## Intrinsic variant calling
    os.makedirs(vcf_dir, exist_ok=True)
    tmp_vcf = vcf_path.replace(".vcf", ".tmp.vcf")

    cmd = f"export OPENBLAS_NUM_THREADS={threads} && \
            bcftools mpileup -a FORMAT/AD,FORMAT/DP -f {ref_genome} {bam_path} | " + \
          f"bcftools call -mv -P 0 -Oz -o {tmp_vcf} && tabix -p vcf {tmp_vcf} "
    executeCmd(cmd, logger=logger)
    cmd = f"export OPENBLAS_NUM_THREADS={threads} && bcftools norm -m -both -f {ref_genome} --multi-overlaps 0 --keep-sum AD -a {tmp_vcf} | \
            bcftools norm -d exact - | \
            bcftools view -i 'ALT!=\"*\"' | \
            bcftools sort -Oz -o {vcf_path} && tabix -f -p vcf {vcf_path} && rm {tmp_vcf}"
    executeCmd(cmd, logger=logger)
    
    logger.info(f"Writing intrinsic VCF to {vcf_path}")
    return bam_path