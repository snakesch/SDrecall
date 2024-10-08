import os
import logging

from preparation.seq import getRawseq
from src.utils import executeCmd, prepare_tmp_file

def getIntrinsicVcf(pc_bed, 
                    all_homo_regions_bed, 
                    counter_bed,
                    pc_masked,
                    ref_genome, 
                    avg_frag_size = 400,
                    std_frag_size = 120,
                    threads = 2,
                    logger = logging.getLogger('SDrecall')):
    # Here fp stands for file_path
    # This function requires an input of:
    # 1. masked genome
    # 2. BED file of PC region
    # 3. BED file of all homologous regions related to the PC region
    # 4. Reference genomef
    # 5. Length stands for the read length

    homo_sd = os.path.dirname(all_homo_regions_bed)
    fastq_dir = homo_sd
    bam_dir = os.path.dirname(pc_bed)
    vcf_dir = os.path.dirname(pc_bed)
    fq_path = os.path.join(fastq_dir, os.path.basename(all_homo_regions_bed)[:-3] + "fastq")
    raw_fq_path = os.path.join(fastq_dir, os.path.basename(all_homo_regions_bed)[:-3] + "raw.fastq")
    bam_path = os.path.join(bam_dir, os.path.basename(pc_bed)[:-3] + "raw.bam")
    vcf_path = bam_path.replace(".bam", ".vcf.gz")
    pc_label = os.path.basename(pc_bed.replace(".bed",""))
        
    # Get Fastq files, note that these reference genome sequences are extracted to single-end sequences intead of paired end sequences.
    # First predetermine the path of the fastq file
    raw_fq_path = getRawseq(all_homo_regions_bed, 
                            raw_fq_path, 
                            ref_genome, 
                            padding = avg_frag_size + std_frag_size)
    logger.info(f"The paired fastq files to get the intrinsic VCF for {pc_bed} is {raw_fq_path}")
    
    # After getting the fastq file it should be used to map against the masked genome and perform straight forward variants calling
    if not os.path.exists(bam_path) or \
       (os.path.getmtime(bam_path) < os.path.getmtime(pc_masked)) or \
       (os.path.getmtime(bam_path) < os.path.getmtime(bash_utils_hub)) or \
       (os.path.getmtime(bam_path) < os.path.getmtime(os.path.abspath(__file__))):
        pc_masked_index = pc_masked.replace(".fasta", ".mmi")
        cmd = f"source {os.path.dirname(__file__)}/shell_utils.sh; independent_minimap2_masked \
               -f {raw_fq_path} \
               -a {pc_masked} \
               -o {bam_path} \
               -s all_PC \
               -t {threads} \
               -i {pc_label} \
               -m asm20 \
               -g {ref_genome}"
        executeCmd(cmd, logger=logger)

    os.makedirs(vcf_dir, exist_ok=True)
    tmp_vcf = vcf_path.replace(".vcf", ".tmp.vcf")

    cmd = f"export OPENBLAS_NUM_THREADS={threads} && \
            bcftools mpileup -a FORMAT/AD,FORMAT/DP -f {ref_genome} {bam_path} | " + \
          f"bcftools call -mv -P 0 -Oz -o {tmp_vcf} && tabix -p vcf {tmp_vcf} "
    executeCmd(cmd, logger=logger)
    cmd = f"export OPENBLAS_NUM_THREADS={threads} && bcftools norm -m -both -f {ref_genome} --multi-overlaps 0 --keep-sum AD -a {tmp_vcf} | \
            bcftools norm -d exact - | \
            bcftools view -i 'ALT!=\"*\"' | \
            bcftools sort --temp-dir /tmp -Oz -o {vcf_path} && tabix -f -p vcf {vcf_path} && rm {tmp_vcf}"
    executeCmd(cmd, logger=logger)
    
    logger.info(f"Writing intrinsic VCF to {vcf_path}")
    return bam_path