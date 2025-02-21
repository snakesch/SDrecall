## Realign --> realignment + variant calling
import logging
import os
from typing import List, Set, Tuple, Dict

import pysam
import pybedtools as pb

from src.utils import executeCmd, prepare_tmp_file

logger = logging.getLogger('SDrecall')


## Helper functions

def merge_bed_files(bed_files: List[str], logger: logging.Logger = logger) -> pb.BedTool:
    """Merges and sorts multiple BED files using pybedtools. Returns BedTool."""
    bed_files = list(dict.fromkeys(bed_files))  # Remove duplicates
    if not bed_files:
        return pb.BedTool("", from_string=True)  # Return empty BedTool

    bedtools_list = [pb.BedTool(f) for f in bed_files]
    merged_bedtool = bedtools_list[0]
    for bt in bedtools_list[1:]:
        merged_bedtool = merged_bedtool.cat(bt, postmerge=False)
    return merged_bedtool.sort().merge()


def select_reads_by_qnames(input_bam: str,
                                   qname_set: Set[str],
                                   output_freads: str,
                                   output_rreads: str,
                                   threads: int = 1,
                                   logger: logging.Logger = logger) -> Tuple[str, str]:
    """Creates forward and reverse FASTQ files from a **qname-sorted** BAM based on QNAMEs.
       Uses pysam to identify read pairs.
    """
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as bamfile:
        with open(output_freads, "w") as f_out, open(output_rreads, "w") as r_out:
            for read in bamfile:
                if read.query_name in qname_set:
                    if read.is_read1:
                        f_out.write(read.to_string(format=pysam.FASTQ) + "\n")
                    elif read.is_read2:
                        r_out.write(read.to_string(format=pysam.FASTQ) + "\n")
                    else:  # Reads could be unpaired.  Handle as needed (e.g., log).
                        logger.warning(f"Read {read.query_name} is neither read1 nor read2.")
    return output_freads, output_rreads




def process_masked_bam( pc_tag,
                        pc_fc_beds,
                        pc_nfc_beds,
                        original_bam,
                        pc_bed,
                        sample_ID,
                        total_freads,
                        total_rreads,
                        max_varno=5,
                        threads=1,
                        mq_cutoff=20,
                        ref_genome="",
                        logger = logger ):


    
    # First merge the FC and NFC regions and masked genomes
    mq_cutoff = int(mq_cutoff)
    pc_dir = os.path.dirname(pc_bed)
    pc_fc_beds = list(dict.fromkeys(pc_fc_beds))
    pc_nfc_beds = list(dict.fromkeys(pc_nfc_beds))
    genome_name = os.path.basename(ref_genome).replace(".fasta", f".sub.{pc_tag}.masked.fasta")
    masked_genome = os.path.join(pc_dir, genome_name)
    pc_fc_bed = pc_bed.replace(".bed", ".targeted.bed")
    pc_nfc_bed = pc_bed.replace(".bed", ".counterparts_regions.targeted.bed")
    cmd = f"cat {' '.join(pc_fc_beds)} | bedtools sort -i - | bedtools merge -i - > {pc_fc_bed}"
    executeCmd(cmd, logger=logger)
    cmd = f"cat {' '.join(pc_nfc_beds)} | bedtools sort -i - | bedtools merge -i - > {pc_nfc_bed}"
    executeCmd(cmd, logger=logger)

    # Now select the reads from the total fastq files
    per_sample_pc_dir = os.path.join(os.path.dirname(original_bam), f"{sample_ID}_ref_SD_recall")
    os.makedirs(per_sample_pc_dir, exist_ok = True)

    sd_freads = os.path.join(per_sample_pc_dir, os.path.basename(total_freads.replace(".fastq", f".{pc_tag}_sdrecall.fastq")))
    sd_rreads = os.path.join(per_sample_pc_dir, os.path.basename(total_rreads.replace(".fastq", f".{pc_tag}_sdrecall.fastq")))

    assert sd_freads != total_freads, f"The output fastq file {sd_freads} should not be the same as the input fastq file {total_freads}"
    assert sd_rreads != total_rreads, f"The output fastq file {sd_rreads} should not be the same as the input fastq file {total_rreads}"

    masked_bam = os.path.join(per_sample_pc_dir, os.path.basename(original_bam.replace(".bam", f".only_{pc_tag}.bam")))
    raw_masked_bam = masked_bam.replace(".bam", ".raw.bam")
    raw_masked_vcf = masked_bam.replace(".bam", ".raw.vcf.gz")

    # Extract the qnames
    pc_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                 pc_fc_bed, 
                                                 output_qnames = "Yes",
                                                 threads=threads,
                                                 multi_aligned=False,
                                                 logger=logger )

    counterpart_qname_lst = bam_reads_selection_by_region(original_bam, 
                                                          pc_nfc_bed, 
                                                          output_qnames = "Yes",
                                                          threads=threads,
                                                          multi_aligned=True,
                                                          mq_cutoff=mq_cutoff,
                                                          logger=logger )

    ## TODO: 1. Sort the input BAM, 2. Select aligned reads by query names, 3. minimap2
    qname_lst_f = prepare_tmp_file().name
    qname_lst_r = prepare_tmp_file().name
    cmdf = f"cat {pc_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_f}"
    cmdr = f"cat {pc_qname_lst} {counterpart_qname_lst} | sort - | uniq > {qname_lst_r}"
    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)

    sd_freads, sd_rreads = extract_reads_by_qnames( qname_lst_f, qname_lst_r,
                                                    total_freads, total_rreads, 
                                                    sd_freads, sd_rreads, 
                                                    logger=logger )

    # Now perform the mapping
    cmd = f"bash {bash_utils_hub} independent_minimap2_masked \
            -b {original_bam} \
            -a {masked_genome} \
            -s {sample_ID} \
            -f {sd_freads} \
            -r {sd_rreads} \
            -g {ref_genome} \
            -o {raw_masked_bam} \
            -t {threads} \
            -i {pc_tag} \
            -c {max_varno} && ls -lh {raw_masked_bam}"
    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        return f"{pc_tag},{pc_bed},NaN,NaN"
    else:
        # Now perform the variants calling
        # First determine the names of the output VCFs
        logger.info(f"Now we have prepared the masked BAM file {raw_masked_bam} for sample {sample_ID} for {pc_bed}, we need to generate the VCF file.")

        sub_running_log = raw_masked_vcf.replace(".vcf.gz", ".log")
        cmd = f"bash {bash_utils_hub} call_polyploidy_per_PC \
                -a {original_bam} \
                -m {ref_genome} \
                -c {threads} \
                -o {raw_masked_vcf} \
                -p {pc_tag} \
                -b {raw_masked_bam} > {sub_running_log} 2>&1"
            # We need to allow this function to be failed sometimes
        try:
            executeCmd(cmd, logger=logger)
        except RuntimeError:
            logger.error(f"Failed to generate {raw_masked_vcf}, running log in {sub_running_log}")
            cmd = f"bash {bash_utils_hub} check_vcf_validity {raw_masked_vcf} 1 && [[ {raw_masked_vcf} -nt {script_path} ]]"
            try:
                executeCmd(cmd, logger=logger)
            except RuntimeError:
                return f"{pc_tag},{pc_bed},{raw_masked_bam},NaN"
            else:
                logger.warning(f"Though the process reported error. The VCF file {raw_masked_vcf} is still valid and updated")
                return f"{pc_tag},{pc_bed},{raw_masked_bam},{raw_masked_vcf}"
        else:
            logger.info(f"Succesfully generate the vcf file {raw_masked_vcf} for region ${pc_bed}")
            return f"{pc_tag},{pc_bed},{raw_masked_bam},{raw_masked_vcf}"

def imap_process_masked_bam(tup_args):
    return process_masked_bam(*tup_args)

def bam_reads_selection_by_region(input_bam: str, 
                                  region_bed: str, 
                                  output_bam = None,
                                  # output_freads = None,
                                  # output_rreads = None,
                                  # input_freads = None,
                                  # input_rreads = None,
                                  multi_aligned = False,
                                  threads = 1,
                                  mq_cutoff = 20, 
                                  gatk = f"bash {bash_utils_hub} gatk_wrapper",
                                  output_qnames = "",
                                  logger=logger):
    if os.path.exists(os.path.dirname(output_qnames)):
        tmp_file = output_qnames
    else:
        tmp_file = prepare_tmp_file().name
    
    if multi_aligned:
        cmd = f"""sambamba view -q -t {threads} -L {region_bed} {input_bam} | \
                  mawk -F '\\t' '($0 !~ /SA:Z:/) && (($5 < 60) || ($0 ~ /XA:Z:/)) {{print $1;}}' | tail -n +1 | sort - | uniq - > {tmp_file}"""
    else:
        cmd = f"""sambamba view -q -t {threads} -L {region_bed} {input_bam} | cut -f 1 | tail -n +1 | sort - | uniq - > {tmp_file}"""
    executeCmd(cmd, logger=logger)
    
    if stat_file_rows(tmp_file) < 1:
        logger.error(f"The region {region_bed} does not have any poor aligned reads in BAM ${input_bam}")
        if output_qnames:
            return None
        elif output_bam:
            return None
        else:
            return None, None
    
    if output_qnames:
        return tmp_file

    if output_bam:
        cmd = f"""{gatk} FilterSamReads -I ${input_bam} -O ${output_bam} --FILTER includeReadList -RLF ${tmp_file} -SO coordinate"""
        executeCmd(cmd, logger=logger)
        os.remove(tmp_file)
        return output_bam

    if output_freads and output_rreads and input_freads and input_rreads:
        cmdf = f"seqtk subseq {input_freads} {tmp_file} > {output_freads}"
        cmdr = f"seqtk subseq {input_rreads} {tmp_file} > {output_rreads}"
        executeCmd(cmdf, logger=logger)
        executeCmd(cmdr, logger=logger)
        # Make sure the output fastq files are consistent regarding the read qnames
        cmd = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
        executeCmd(cmd, logger=logger)
 
        return output_freads, output_rreads

def stat_file_rows(file):
    with open(file, "r") as f:
        lines = f.readlines()
    return len(lines)

def extract_reads_by_qnames(qname_lst_f, qname_lst_r,
                            input_freads, input_rreads, 
                            output_freads, output_rreads, 
                            logger = logger):

    cmdf = f"seqtk subseq {input_freads} {qname_lst_f} > {output_freads}"
    cmdr = f"seqtk subseq {input_rreads} {qname_lst_r} > {output_rreads}"

    executeCmd(cmdf, logger=logger)
    executeCmd(cmdr, logger=logger)
    # Make sure the output fastq files are consistent regarding the read qnames
    cmds = f"bash {bash_utils_hub} sync_fastq -1 {output_freads} -2 {output_rreads}"
    executeCmd(cmds, logger=logger)

    return output_freads, output_rreads
