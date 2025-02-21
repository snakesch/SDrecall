#!/usr/bin/env python3

from glob import glob
import pandas as pd
import os
from typing import List, Tuple, Optional
import logging
from collections import defaultdict

logger = logging.getLogger('SDrecall')

import pybedtools as pb

pb.set_bedtools_path("/share1/bedtools/2.30.0/bin/")
work_dir = "/home/snakesch/work/SDrecall/test/HG002/test_run_output"
reference = "/home/snakesch/work/SDrecall/test/ucsc.hg19.fasta"

## TODO: Get fragment size distribution from file
avg_frag_size, std_frag_size = 570.4, 150.7

def enumerate_PCs(work_dir):
    # Prepare arguments for downstream imap_prepare_masked_align_region_per_PC
    beds = glob(f"{work_dir}/PC*_related_homo_regions/PC*_all/PC*_related_homo_regions.bed")
    
    # - For load balancing: high coverage PCs top on the list - #
    beds = sorted(beds, key = lambda f: pb.BedTool(f).total_coverage(), reverse=True)

    # Extract FC regions from the total bed
    pc_names, pc_units = [], []
    for bed_path in beds:
        pc_names.append(os.path.basename(bed_path).split("_")[0])

        fc_regions = []
        for interval in pb.BedTool(bed_path):
            if interval[-1].startswith("FC:"):
                fc_regions.append(interval)
        
        # Sort FC regions by their tag (last field). 
        fc_regions.sort(key=lambda interval: interval[-1])

        pc_units.append(tuple(range(len(fc_regions))))
    print(pc_names, pc_units)
    return pc_names, pc_units

def prepare_masked_align_region_per_PC(
    pc_tag: str,
    pc_subids: Tuple[str, ...],
    sd_map_dir: str,
    target_region_bed: str,
    ref_genome: str = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
    logger: logging.Logger = logger,
) -> Tuple[str, ...]:

    whole_region_raw_bed = os.path.join(
        sd_map_dir, f"{pc_tag}_related_homo_regions", f"{pc_tag}_related_homo_regions.raw.bed"
    )
    logger.info(f"The raw bed file for this {pc_tag} is {whole_region_raw_bed}")

    # Efficiently filter the raw BED file using pybedtools *directly*
    fc_bed = pb.BedTool(whole_region_raw_bed).filter(lambda interval: interval[6].startswith("FC")).saveas()

    target_regions = pb.BedTool(target_region_bed)
    ref_genome_fai = ref_genome + ".fai"

    # Directly create and use a temporary file for the targeted FC regions
    target_fc_bed = (
        fc_bed.intersect(target_regions)
        .slop(b=500, g=ref_genome_fai)
        .sort()
        .merge()
        .saveas()  # saveas() returns a BedTool object.
    )
    target_fc_size = target_fc_bed.total_coverage()

    records = []
    for subgroup_id in pc_subids:
        logger.info(
            f"Start to fetch the masked align region for NFC regions for PC {pc_tag} subgroup {subgroup_id}"
        )

        # Filter for NFC and FC regions directly using pybedtools.
        nfc_bed = pb.BedTool(whole_region_raw_bed).filter(
            lambda interval: interval[6] == f"NFC:{pc_tag}_{subgroup_id}"
        ).saveas()

        fc_region_bed = pb.BedTool(whole_region_raw_bed).filter(
            lambda x: x.name == f"FC:{pc_tag}_{subgroup_id}"
        ).saveas()


        # Directly use the result of extract_and_pad.
        targeted_nfc_regions_bed, target_nfc_regions_bed_size = extract_and_pad_segments(
            fc_region_bed,
            nfc_bed,
            target_regions,
            padding=600,
            logger=logger,
        )

        record = (
            f"{pc_tag},{subgroup_id},{target_fc_bed.fn},{targeted_nfc_regions_bed.fn},"
            f"{target_fc_size},{target_nfc_regions_bed_size}"
        )
        logger.info(f"The returned record for {pc_tag} subgroup {subgroup_id} is {record}")
        records.append(record)

    return tuple(records)


def extract_and_pad_segments(
    single_interval_bed: pb.BedTool,
    multiple_intervals_bed: pb.BedTool,
    target_regions_bed: pb.BedTool,
    padding: int = 600,
    logger: logging.Logger = logger,
) -> Tuple[pb.BedTool, int]:
    """
    Extracts, pads, and merges segments from multiple_intervals_bed that overlap
    target regions within a single interval in single_interval_bed.
    """

    overlapping_segments = single_interval_bed.intersect(target_regions_bed)
    main_interval = single_interval_bed[0]
    main_start, main_end, main_strand = main_interval.start, main_interval.end, main_interval.strand

    # 1. Merge overlapping segments *after* padding.  Crucially, we now keep track of
    #    the *original* start/end for merging purposes, and use the *padded*
    #    start/end for BedTool creation.
    merged_segments = []
    for segment in overlapping_segments:
        padded_start = max(segment.start - padding, 0)
        padded_end = segment.end + padding
        merged_segments.append((segment.start, segment.end, padded_start, padded_end))  # (orig_start, orig_end, padded_start, padded_end)

    merged_segments.sort(key=lambda x: x[0])  # Sort by original start

    if not merged_segments:
        return pb.BedTool("", from_string=True), 0
    
    final_merged = [merged_segments[0]]
    for cur_orig_start, cur_orig_end, cur_padded_start, cur_padded_end in merged_segments[1:]:
        last_orig_start, last_orig_end, last_padded_start, last_padded_end = final_merged[-1]
        if cur_orig_start <= last_orig_end:  # Use *original* coords for overlap
            final_merged[-1] = (
                last_orig_start,
                max(last_orig_end, cur_orig_end),
                min(last_padded_start, cur_padded_start),
                max(last_padded_end, cur_padded_end)
            )
        else:
            final_merged.append((cur_orig_start, cur_orig_end, cur_padded_start, cur_padded_end))
            

    # 2. Extract corresponding segments, adjusting for strand and relative coords.
    corresponding_segments = []
    for interval in multiple_intervals_bed:
        interval_strand = interval.strand
        rel_start_interval = int(interval[3])
        rel_end_interval = int(interval[4])


        for orig_start, orig_end, padded_start, padded_end in final_merged:

            # Convert padded start/end to relative coords
            rel_start = padded_start - main_start
            rel_end = padded_end - main_start

            # Calculate overlap *after* converting to relative coords
            overlap_start = max(rel_start, rel_start_interval)
            overlap_end = min(rel_end, rel_end_interval)

            if overlap_start < overlap_end:  # Only proceed if there's actual overlap
                if interval_strand == main_strand:
                    rel_small_start = overlap_start - rel_start_interval
                    rel_small_end = overlap_end - rel_start_interval
                else:
                    rel_small_start = rel_end_interval - overlap_end
                    rel_small_end = rel_end_interval - overlap_start

                abs_start = max(interval.start, interval.start + rel_small_start)
                abs_end = min(interval.end, interval.start + rel_small_end)
                if abs_start < abs_end:
                    corresponding_segments.append(
                        pb.Interval(interval.chrom, abs_start, abs_end, strand=interval_strand, name=interval.name)
                    )

    return_bed = pb.BedTool(corresponding_segments).sort().merge()
    return_bed_size = return_bed.total_coverage()
    return return_bed, return_bed_size

prepared_arguments_df = enumerate_PCs(work_dir)
uniq_pc_names = prepared_arguments_df.loc[:, 0].drop_duplicates().tolist()

print(uniq_pc_names)
pc_subids_tup_list = [tuple(prepared_arguments_df.loc[prepared_arguments_df[0] == pc, 1].drop_duplicates().tolist()) for pc in uniq_pc_names]
print(pc_subids_tup_list)
# pool = ctx.Pool(threads, initializer=pool_init)
# prepared_beds = pool.imap_unordered(imap_prepare_masked_align_region_per_PC, zip(uniq_pc_names,
#                                                                                 pc_subids_tup_list,
#                                                                                 repeat(all_PC_folder),
#                                                                                 repeat(target_region),
#                                                                                 repeat(ref_fasta)))