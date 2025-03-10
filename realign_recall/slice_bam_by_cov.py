import os
import sys
import re
import multiprocessing as mp
ctx = mp.get_context("spawn")

from itertools import repeat
import pybedtools as pb
    
from src.utils import executeCmd, prepare_tmp_file
from src.log import logger, log_command
from src.const import shell_utils


def imap_slice_bam_per_bed(tup_args):
    return slice_bam_per_bed(*tup_args)


@log_command
def slice_bam_per_bed(bed, bam, ref_genome, chunk_id, threads = 4, tmp_dir = "/tmp", logger = logger):
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

    if os.path.exists(bed):
        # Coordinates in bed file are 0-indexed and half-open
        cmd = f"samtools view -@ {threads} -u -P -L {bed} {bam} | \
                samtools sort -@ {threads} -T {tmp_dir}/{chunk_id}_{os.path.basename(cov_bam)} -o {cov_bam} -O bam - && \
                bash {shell_utils} modify_bam_sq_lines {cov_bam} {ref_genome} {cov_bam_header} && \
                samtools reheader {cov_bam_header} {cov_bam} > {cov_bam}.tmp && \
                mv {cov_bam}.tmp {cov_bam} && \
                samtools index {cov_bam}"
    elif re.match(r"^chr.*\d+$", bed):
        region = bed # When input the region str, both start and end are 1-indexed and inclusive
        cmd = f"samtools view -@ {threads} -u -P {bam} {region} | \
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


def extract_depth_blocks(depth_file, min_depth=5, output_file=None):
    """
    Extract continuous blocks with at least the specified minimum depth.
    
    Args:
        depth_file (str): Path to the samtools depth output file
        min_depth (int): Minimum depth threshold (default: 5)
        output_file (str): Output file path (default: stdout)
    """
    out = open(output_file, 'w') if output_file else sys.stdout
    
    current_block = None
    
    with open(depth_file, 'r') as f:
        for line in f:
            chrom, pos, depth = line.strip().split()
            pos = int(pos)
            depth = int(depth)
            
            if depth >= min_depth:
                # We're in a valid region
                if current_block is None:
                    # Start a new block
                    current_block = [chrom, pos, pos]
                elif current_block[0] == chrom and current_block[2] == pos - 1:
                    # Extend current block
                    current_block[2] = pos
                else:
                    # Output previous block and start a new one
                    out.write(f"{current_block[0]}\t{current_block[1]}\t{current_block[2] + 1}\n")
                    current_block = [chrom, pos, pos]
            else:
                # We're not in a valid region
                if current_block is not None:
                    # Output previous block and reset
                    out.write(f"{current_block[0]}\t{current_block[1]}\t{current_block[2] + 1}\n")
                    current_block = None
    
    # Don't forget the last block if we have one
    if current_block is not None:
        out.write(f"{current_block[0]}\t{current_block[1]}\t{current_block[2] + 1}\n")
    
    if output_file:
        out.close()
        return output_file


def process_target_regions_with_coverage(target_bed, 
                                         cov_bed, 
                                         min_interval=2000, 
                                         max_interval=10000,
                                         delimiter_size = 500,
                                         ref_genome = "",
                                         tmp_dir = "/tmp",
                                         logger=logger):
    """
    Process target regions based on coverage data:
    1. Large intervals (>10000bp) are split using coverage data
    2. Small intervals (<2000bp) are merged if they fall within a coverage island (2000-5000bp)
    
    Args:
        target_bed (str): BED file with target regions
        cov_bed (str): BED file with coverage data
        min_interval (int): Minimum interval size to consider (default: 2000bp)
        max_interval (int): Maximum interval size before splitting (default: 10000bp)
        
    Returns:
        list: Processed intervals for BAM splitting
    """
    # Load the target and coverage BEDs
    target_regions = pb.BedTool(target_bed).sort().merge(d = delimiter_size)
    coverage_regions = pb.BedTool(cov_bed)
    
    # Track our processed intervals
    processed_intervals = []
    small_intervals = []
    
    # Process each target region
    for region in target_regions:
        region_size = region.end - region.start
        
        # Case 1: Large interval that needs potential splitting
        if region_size > max_interval:
            logger.debug(f"Large region found: {region.chrom}:{region.start}-{region.end} ({region_size}bp)")
            
            # tmp_bed = prepare_tmp_file(suffix=".bed", tmp_dir = tmp_dir).name
            # Create a single-region BED for this large region
            temp_region = pb.BedTool([region])
            
            # Intersect with coverage data to find covered subregions
            covered_subregions = temp_region.intersect(coverage_regions)
            
            if len(covered_subregions) > 1:
                # We found covered subregions - use these instead
                logger.debug(f"  - Split into {len(covered_subregions)} covered subregions")
                for subregion in covered_subregions:
                    processed_intervals.append(subregion)
            else:
                # No coverage data - keep original but log warning
                logger.debug(f"  - No coverage data for large region {region.chrom}:{region.start}-{region.end}")
                processed_intervals.append(region)
        
        # Case 2: Small interval that might need merging
        elif region_size < min_interval:
            logger.debug(f"Small region found: {region.chrom}:{region.start}-{region.end} ({region_size}bp)")
            small_intervals.append(region)
        # Case 3: Medium interval that's already appropriate size
        else:
            processed_intervals.append(region)


    # Deal with small intervals in batch
    # Intersect small bed with coverage bed and pick the coverage islands that covered multiple small intervals
    covered_small_intervals = pb.BedTool(small_intervals).intersect(coverage_regions, wao = True).to_dataframe( disable_auto_names = True, 
                                                                                                                names = ["chrom_target", "start_target", "end_target", 
                                                                                                                         "chrom_cov", "start_cov", "end_cov", 
                                                                                                                         "overlap_len"] )
    logger.debug(f"These small intervals are covered by reads:\n{covered_small_intervals.to_string(index=False)}\n")
    groupby_cov_island = covered_small_intervals.groupby(["chrom_cov", "start_cov", "end_cov"], as_index = False)
    for cov_island, df in groupby_cov_island:
        if df.overlap_len.values[0] == 0:
            logger.warning(f"These small targeting regions are not covered by sufficient reads:\n{df[['chrom_target', 'start_target', 'end_target']].to_string(index=False)}\n")
            continue
        # First check whether the cov_island is appropriately large
        if cov_island[2] - cov_island[1] >= min_interval and cov_island[2] - cov_island[1] <= max_interval:
            # Directly input the cov_island into processed_intervals
            processed_intervals.append((cov_island[0], cov_island[1], cov_island[2]))
        elif cov_island[2] - cov_island[1] < min_interval:
            # cov_island is too small, just let pad the overlapping small intervals by delimiter_size and merge
            processed_intervals.append((cov_island[0], cov_island[1] - delimiter_size, cov_island[2] + delimiter_size))
        else:
            # cov_island is too large, just let pad the overlapping small intervals by delimiter_size and merge
            cov_small_regions = df[["chrom_target", "start_target", "end_target"]]
            cov_small_regions = pb.BedTool.from_dataframe(cov_small_regions)
            cov_small_regions = cov_small_regions.sort().slop(b = delimiter_size, g = f"{ref_genome}.fai").merge(d = delimiter_size)
            for region in cov_small_regions:
                processed_intervals.append((region.chrom, region.start, region.end))
            
    result_bed = pb.BedTool(processed_intervals).sort().merge()
    
    return result_bed


def split_bam_by_cov(bam, 
                     target_bed = None,
                     min_depth = 3,
                     min_interval = 2000,
                     max_interval = 10000,
                     delimiter_size = 1000,
                     base_qual = 0,
                     map_qual = 0,
                     beds = [],
                     threads = 4,
                     ref_genome = "",
                     tmp_dir = "/tmp",
                     logger = logger):
    """
    Split a BAM file based on coverage and target regions using the new strategy:
    1. Start with target regions
    2. Process large regions (>10000bp) by splitting with coverage data
    3. Process small regions (<2000bp) by merging with coverage islands
    4. Split the BAM file based on the processed regions
    
    Args:
        bam (str): Input BAM file
        target_bed (str, optional): BED file with target regions
        min_depth (int): Minimum depth for coverage regions (default: 5)
        min_interval (int): Minimum interval size (default: 2000)
        max_interval (int): Maximum interval size (default: 10000)
        delimiter_size (int): Size for merging nearby regions (default: 1000)
        base_qual (int): Minimum base quality for coverage calculation
        map_qual (int): Minimum mapping quality for coverage calculation
        beds (list): Pre-defined BED files to use instead of computing intervals
        threads (int): Number of threads to use
        ref_genome (str): Reference genome path
        tmp_dir (str): Temporary directory
        
    Returns:
        tuple: (list of split BAM files, list of BED files used for splitting)
    """
    if len(beds) == 0 and target_bed is None:
        raise ValueError("Either beds or target_bed must be provided")
    
    if len(beds) == 0:
        # First, get coverage information
        cov_depth_file = bam.replace(".bam", ".depth")
        cov_bed = bam.replace(".bam", ".cov.bed")
        
        # Use samtools to get coverage depth
        cmd = f"[[ {cov_depth_file} -nt {bam} ]] || samtools depth -q {base_qual} -Q {map_qual} {bam} > {cov_depth_file}"
        executeCmd(cmd, logger=logger)
        
        # Extract blocks with minimum depth
        extract_depth_blocks(cov_depth_file, min_depth=min_depth, output_file=cov_bed)
        logger.info(f"Extracted coverage regions to {cov_bed}")
        
        # Process target regions with coverage data
        processed_intervals = process_target_regions_with_coverage(
            target_bed=target_bed,
            cov_bed=cov_bed,
            min_interval=min_interval,
            max_interval=max_interval,
            delimiter_size=delimiter_size,
            ref_genome=ref_genome,
            tmp_dir=tmp_dir,
            logger=logger
        )
        
        # Save processed regions to a BED file
        processed_bed = bam.replace(".bam", ".processed.bed")
        processed_intervals.saveas(processed_bed)
        for interval in processed_intervals:
            beds.append(interval.chrom + ":" + str(interval.start + 1) + "-" + str(interval.end))
        logger.info(f"Saved processed regions to {processed_bed} and it has {len(processed_intervals)} intervals")
    else:
        logger.info(f"Using provided coverage bed files ({len(beds)} beds) to split the BAM file {bam}")

    # Slice the BAM file based on the BED files
    with ctx.Pool(threads - 1) as pool:
        result_records = pool.imap(imap_slice_bam_per_bed, zip(
            beds, 
            repeat(bam),
            repeat(ref_genome), 
            [i for i in range(1, len(beds) + 1)],
            repeat(1), 
            repeat(tmp_dir)))

        splitted_bams = []
        i = 0
        for success, result, log_contents in result_records:
            i += 1
            print(f"\n************************************{i}_subprocess_start_for_slice_bam************************************\n", file=sys.stderr)
            if success:
                splitted_bams.append(result)
                print(f"Successfully slice the bam file {bam}. The log info are:\n{log_contents}\n", file=sys.stderr)
            else:
                error_mes, tb_str = result
                logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")

            print(f"\n************************************{i}_subprocess_end_for_slice_bam************************************\n", file=sys.stderr)
     
    logger.info(f"Successfully split the BAM file {bam} into {len(splitted_bams)} parts")
    return splitted_bams, beds