


from src.log import logger, log_command
from src.const import shell_utils
from src.utils import executeCmd, merge_bed_files
from realign_recall.read_extraction import bam_to_fastq_biobambam


def imap_process_masked_bam(tup_args):
    return process_masked_bam(*tup_args)


@log_command
def process_masked_bam( rg_tag,
                        rg_bed,
                        rg_fc_beds,
                        rg_nfc_beds,
                        original_bam,
                        raw_masked_bam,
                        masked_genome,
                        sample_ID,
                        sd_freads,
                        sd_rreads,
                        threads=1,
                        ref_genome="",
                        logger = logger ):

    # First merge the FC and NFC regions and masked genomes  
    rg_fc_beds = list(dict.fromkeys(rg_fc_beds))
    rg_nfc_beds = list(dict.fromkeys(rg_nfc_beds))

    rg_fc_bed = rg_bed.replace(".bed", ".targeted.bed")
    rg_nfc_bed = rg_bed.replace(".bed", ".counterparts_regions.targeted.bed")

    rg_fc_bed = merge_bed_files(rg_fc_beds, output_bed=rg_fc_bed)
    rg_nfc_bed = merge_bed_files(rg_nfc_beds, output_bed=rg_nfc_bed)

    raw_masked_vcf = raw_masked_bam.replace(".bam", ".vcf.gz")

    cmd = f"bash {shell_utils} quick_check_bam_validity {raw_masked_bam} && \
            [[ {raw_masked_bam} -nt {sd_freads} ]] && \
            [[ {raw_masked_bam} -nt {sd_rreads} ]] && \
            bash {shell_utils} check_vcf_validity {raw_masked_vcf} 1 && \
            [[ {raw_masked_vcf} -nt {raw_masked_bam} ]] && \
            [[ {sd_freads} -nt {rg_bed} ]] && \
            [[ {sd_freads} -nt {original_bam} ]]"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.info(f"The masked bam {raw_masked_bam} is not valid, so we need to generate it.")
        execute = True
    else:
        logger.info(f"The masked bam {raw_masked_bam} is already up-to-date. No need to map again.")
        execute = False

    if execute:
        # Extract the reads directly from the original bam
        fc_sd_freads, fc_sd_rreads = sd_freads.replace(".fastq", ".fc.fastq"), sd_rreads.replace(".fastq", ".fc.fastq")
        nfc_sd_freads, nfc_sd_rreads = sd_freads.replace(".fastq", ".nfc.fastq"), sd_rreads.replace(".fastq", ".nfc.fastq")
        fc_sd_freads, fc_sd_rreads = bam_to_fastq_biobambam(original_bam, 
                                                            rg_fc_bed, 
                                                            fc_sd_freads, 
                                                            fc_sd_rreads, 
                                                            threads=threads, 
                                                            logger=logger)
        
        nfc_sd_freads, nfc_sd_rreads = bam_to_fastq_biobambam(original_bam, 
                                                            rg_nfc_bed, 
                                                            nfc_sd_freads, 
                                                            nfc_sd_rreads, 
                                                            multi_aligned=True,
                                                            threads=threads, 
                                                            logger=logger)
        
        cmd = f"cat {fc_sd_freads} {nfc_sd_freads} > {sd_freads} && \
                cat {fc_sd_rreads} {nfc_sd_rreads} > {sd_rreads}"
        executeCmd(cmd, logger=logger)
        
        # Now perform the mapping
        cmd = f"bash {shell_utils} \
                independent_minimap2_masked \
                -a {masked_genome} \
                -s {sample_ID} \
                -f {sd_freads} \
                -r {sd_rreads} \
                -g {ref_genome} \
                -o {raw_masked_bam} \
                -t {threads} \
                -i {rg_tag} && ls -lh {raw_masked_bam}"
        try:
            executeCmd(cmd, logger=logger)
        except RuntimeError:
            return f"{rg_tag},{rg_bed},NaN,NaN"
        else:
            # Now perform the variants calling
            # First determine the names of the output VCFs
            logger.info(f"Now we have prepared the masked BAM file {raw_masked_bam} for sample {sample_ID} for {rg_bed}, we need to generate the VCF file.")
            
            sub_running_log = raw_masked_vcf.replace(".vcf.gz", ".log")
            cmd = f"bash {shell_utils} bcftools_call_per_RG \
                    -m {ref_genome} \
                    -c {threads} \
                    -o {raw_masked_vcf} \
                    -p {rg_tag} \
                    -b {raw_masked_bam}"
                # We need to allow this function to be failed sometimes
            try:
                executeCmd(cmd, logger=logger)
            except RuntimeError:
                logger.error(f"Failed to generate {raw_masked_vcf}, running log in {sub_running_log}")
                cmd = f"bash {shell_utils} check_vcf_validity {raw_masked_vcf} 1 && [[ {raw_masked_vcf} -nt {script_path} ]]"
                try:
                    executeCmd(cmd, logger=logger)
                except RuntimeError:
                    return f"{rg_tag},{rg_bed},{raw_masked_bam},NaN"
                else:
                    logger.warning(f"Though the process reported error. The VCF file {raw_masked_vcf} is still valid and updated")
                    return f"{rg_tag},{rg_bed},{raw_masked_bam},{raw_masked_vcf}"
            else:
                logger.info(f"Succesfully generate the vcf file {raw_masked_vcf} for region {rg_bed}")
                return f"{rg_tag},{rg_bed},{raw_masked_bam},{raw_masked_vcf}"
    else:
        logger.info(f"The masked bam {raw_masked_bam} is already up-to-date. No need to map again.")
        return f"{rg_tag},{rg_bed},{raw_masked_bam},{raw_masked_vcf}"