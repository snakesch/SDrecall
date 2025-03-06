import os
import sys
import multiprocessing as mp
ctx = mp.get_context("spawn")

from itertools import repeat

from src.utils import executeCmd, prepare_tmp_file
from src.log import logger, log_command
from src.const import shell_utils


def imap_slice_bam_per_bed(tup_args):
    return slice_bam_per_bed(*tup_args)


@log_command
def slice_bam_per_bed(bed, bam, ref_genome, threads = 4, tmp_dir = "/tmp", logger = logger):
    # assert os.path.basename(bed).split(".")[-2] != chunk_id, f"The bed file {bed} should not have the same chunk_id as the previous bed file"
    chunk_id = os.path.basename(bed).split(".")[-2]
    cov_bam = bam.replace(".bam", f".{chunk_id}.bam")
    cov_bam_header = cov_bam.replace(".bam", ".header")
    bam_index = bam + ".bai"
    if not os.path.exists(bam_index):
        execute = True
    elif os.path.getmtime(bam) > os.path.getmtime(bam_index):
        execute = True
    else:
        execute = False

    if execute:
        tmp_index = prepare_tmp_file(suffix=".bai", tmp_dir = tmp_dir).name
        cmd = f"samtools index -b -o {tmp_index} {bam} && \
                [[ -f {bam_index} ]] && \
                [[ {bam_index} -nt {bam} ]] || \
                mv -f {tmp_index} {bam_index}"
        executeCmd(cmd, logger = logger)
        # We need to make sure the index update is atomic

    cmd = f"samtools view -@ {threads} -u -L {bed} {bam} | \
            samtools sort -@ {threads} -T {tmp_dir}/{chunk_id}_{os.path.basename(cov_bam)} -o {cov_bam} -O bam - && \
            bash {shell_utils} modify_bam_sq_lines {cov_bam} {ref_genome} {cov_bam_header} && \
            samtools reheader {cov_bam_header} {cov_bam} > {cov_bam}.tmp && \
            mv {cov_bam}.tmp {cov_bam} && \
            samtools index {cov_bam}"
    try:
        executeCmd(cmd, logger = logger)
    except RuntimeError:
        logger.warning(f"Failed to slice the bam file {bam} by bed {bed} and generate a {cov_bam}")
        return "NaN"
    else:
        return cov_bam


def split_bed_by_size(bed_file, 
                      chunksize=100000, 
                      delimiter=1000,
                      logger = logger):
    import pybedtools
    # Read the BED file
    bed = pybedtools.BedTool(bed_file)
    
    # Sort the BED file by chromosome and start position
    bed = bed.sort().merge(d=delimiter)
    
    # Initialize variables
    chunk_num = 1
    chunk_size = 0
    chunk_intervals = []
    
    # Iterate over the intervals
    for interval in bed:
        interval_size = interval.end - interval.start
        
        if chunk_size + interval_size > chunksize:
            # Write the current chunk to a new file
            if len(chunk_intervals) > 0:
                chunk_bed = pybedtools.BedTool(chunk_intervals)
                chunk_bed.saveas(bed_file.replace(".bed", f".chunk{chunk_num}.bed"))
            
                # Reset variables for the next chunk
                chunk_size = interval_size
                chunk_intervals = [interval]
            else:
                logger.warning("The interval size is larger than the chunk size, so we output the interval as a single chunk: {}".format(interval))
                chunk_bed = pybedtools.BedTool([interval])
                chunk_bed.saveas(bed_file.replace(".bed", f".chunk{chunk_num}.bed"))
                chunk_size = 0
                chunk_intervals = []
            
            chunk_num += 1
        else:
            # Add the interval to the current chunk
            chunk_size += interval_size
            chunk_intervals.append(interval)
    
    # Write the last chunk to a file
    if len(chunk_intervals) > 0:
        chunk_bed = pybedtools.BedTool(chunk_intervals)
        chunk_bed.saveas(bed_file.replace(".bed", f".chunk{chunk_num}.bed"))

    return [bed_file.replace(".bed", f".chunk{i}.bed") for i in range(1, chunk_num + 1)]


def split_bam_by_cov(bam, 
                     chunksize=10000, 
                     delimiter_size=1000,
                     base_qual = 0,
                     map_qual = 0,
                     beds = [],
                     threads = 4,
                     ref_genome = "",
                     tmp_dir = "/tmp",
                     logger = logger):
    
    if len(beds) == 0:
        cov_bed = bam.replace(".bam", f".bed")
        cmd = f"bash {shell_utils} samtools_bam_coverage \
                -i {bam} \
                -d {delimiter_size} \
                -o {cov_bed} \
                -b {base_qual} \
                -m {map_qual}"
        executeCmd(cmd, logger = logger)

        beds = split_bed_by_size(cov_bed, 
                                chunksize = chunksize, 
                                logger = logger)

        logger.info(f"Getting the coverage bed file {cov_bed} and split it into {len(beds)} parts")
    else:
        logger.info(f"Using the provided coverage bed files({len(beds)} beds) to split the BAM file {bam}")

    with ctx.Pool(threads - 1) as pool:
        result_records = pool.imap(imap_slice_bam_per_bed, zip( beds, 
                                                                repeat(bam),
                                                                repeat(ref_genome), 
                                                                repeat(1), 
                                                                repeat(tmp_dir)))

        splitted_bams = []
        i=0
        for success, result, log_contents in result_records:
            i+=1
            print(f"\n************************************{i}_subprocess_start_for_slice_bam************************************\n", file=sys.stderr)
            if success:
                splitted_bams.append(result)
                print(f"Successfully slice the bam file {bam}. The log info are:\n{log_contents}\n", file=sys.stderr)
            else:
                error_mes, tb_str = result
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")

            print(f"\n************************************{i}_subprocess_end_for_slice_bam************************************\n", file=sys.stderr)
     
    logger.info(f"Successfully splitted the BAM file {bam} into {len(splitted_bams)} parts")
    return splitted_bams, beds