import os
import re
import pysam
import intervaltree

from collections import defaultdict

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
                    tmp_dir = "/tmp",
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
    getRawseq(all_homo_regions_bed, intrinsic_fastq, ref_genome, tmp_dir = tmp_dir, padding = avg_frag_size + std_frag_size)
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



def compute_interval_status(bam_file, logger):
    """
    Computes allowed status for each distinct interval extracted from the qname.
    An interval is represented as a tuple (chrom, start, end). If an interval 
    is enclosed by a larger interval on the same chromosome, it is marked False.
    For identical intervals, only one copy will be considered allowed.
    """
    # Updated regex: the :label part is optional.
    qname_regex = re.compile(r'(.*):(\d+)-(\d+):(.+)')
    intervals_by_chr = defaultdict(set)
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            # Only consider reads with a query sequence and mapped
            if read.is_unmapped or read.query_sequence is None:
                continue
            m = qname_regex.match(read.query_name)
            if m:
                chrom = m.group(1)
                start_val = int(m.group(2))
                end_val = int(m.group(3))
                intervals_by_chr[chrom].add((start_val, end_val))
    
    allowed_status = {}
    # Process intervals per chromosome
    for chrom, intervals in intervals_by_chr.items():
        # Sort by start ascending; for equal starts, sort by end descending.
        sorted_intervals = sorted(list(intervals), key=lambda x: (x[0], -x[1]))
        max_end = -1
        for i, (s, e) in enumerate(sorted_intervals):
            # The first interval is always allowed.
            if i == 0:
                allowed_status[(chrom, s, e)] = True
                max_end = e
            else:
                # If current interval's end is <= max_end, it is enclosed by a previous interval.
                if e <= max_end:
                    allowed_status[(chrom, s, e)] = False
                else:
                    allowed_status[(chrom, s, e)] = True
                    max_end = e
    logger.info(f"Computed interval allowed statuses for {bam_file}")
    return allowed_status



def filter_intrinsic_alignments(bam_file, output_file=None, logger=logger):
    """
    Filter out alignments where the start position in the read name matches 
    the alignment position in the reference (i.e. mapping to its original 
    genomic location) and now also filter out reads whose intervals (from the 
    qname) are either enclosed by larger intervals or are duplicates.
    
    For reads with qnames matching the pattern: chr:start-end(:label)?,
    only one copy per distinct interval (chr, start, end) will be kept 
    if that interval is not enclosed by a larger one.
    
    Parameters:
    -----------
    bam_file : str
        Path to input BAM file
    output_file : str, optional
        Path to output filtered BAM file
    """
    # Updated regex: the :label part is not required.
    qname_regex = re.compile(r'(.*):(\d+)-(\d+):(.+)')
    tmp_output = prepare_tmp_file(suffix=".bam").name
    
    # Compute allowed status for each distinct interval.
    allowed_status = compute_interval_status(bam_file, logger)
    seen_intervals = set()  # To track intervals already allowed (one copy per distinct interval)
    
    with pysam.AlignmentFile(bam_file, "rb") as input_bam:
        with pysam.AlignmentFile(tmp_output, "wb", header=input_bam.header) as output_bam:
            total_reads = 0
            primary_align_origin_qnames = set()
            sec_to_pri_qnames = set() 
            buffer_sec_aligns = {}
            
            for read in input_bam:
                total_reads += 1
                
                # Skip unmapped or reads with no query sequence.
                if read.is_unmapped or read.query_sequence is None:
                    continue
                
                qname = read.query_name
                # Check the interval filtering only if the qname matches the expected pattern.
                m = qname_regex.match(qname)
                if m:
                    chrom = m.group(1)
                    start_val = int(m.group(2))
                    end_val = int(m.group(3))
                    interval = (chrom, start_val, end_val)
                    # If the interval is in our computed dictionary and is enclosed, skip.
                    if not allowed_status.get(interval, True):
                        continue
                    # If the interval has already been seen, skip duplicate.
                    if interval in seen_intervals:
                        continue
                    # Otherwise, mark this interval as seen.
                    seen_intervals.add(interval)
                
                    # Now follow the original intrinsic alignment filtering.
                    if qname in primary_align_origin_qnames:
                        if read.is_supplementary:
                            continue
                        elif read.is_secondary:
                            sec_to_pri_qnames.add(qname)
                            # Remove the secondary alignment flag (256)
                            read.flag = read.flag - 256
                            output_bam.write(read)
                            continue
                    elif read.is_supplementary:
                        continue
                    elif read.is_secondary:
                        buffer_sec_aligns[qname] = read
                        continue
                    # Compare with alignment position (pysam.read.reference_start is 0-based)
                    if start_val != read.reference_start + 1:
                        output_bam.write(read)
                    else:
                        primary_align_origin_qnames.add(qname)
                else:
                    # If qname doesn't match the expected pattern, keep the read.
                    output_bam.write(read)
            
            # Process buffered secondary alignments that have a primary counterpart.
            primary_align_origin_qnames = primary_align_origin_qnames - sec_to_pri_qnames
            for qname in primary_align_origin_qnames:
                if qname in buffer_sec_aligns:
                    sec_read = buffer_sec_aligns[qname]
                    sec_read.flag = sec_read.flag - 256
                    output_bam.write(sec_read)
                    sec_to_pri_qnames.add(qname)

            logger.info(
                f"Filtered out {len(primary_align_origin_qnames)} reads (intrinsic alignments) "
                f"and duplicates/enclosed intervals (duplicate intervals: {len(seen_intervals)} kept) "
                f"out of {total_reads} reads in {bam_file}. "
                f"Converted {len(sec_to_pri_qnames)} secondary alignments to primary alignments."
            )
    
    # Index the output BAM.
    if output_file is None:
        cmd = f"mv {tmp_output} {bam_file} && samtools index {bam_file}"
        executeCmd(cmd, logger=logger)
        return bam_file
    else:
        cmd = f"mv {tmp_output} {output_file} && samtools index {output_file}"
        executeCmd(cmd, logger=logger)
        return output_file

