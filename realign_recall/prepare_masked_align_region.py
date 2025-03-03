import pybedtools as pb
import pandas as pd
import os


from src.log import logger, log_command
from src.utils import prepare_tmp_file

def extract_and_pad_segments(single_interval_bed, 
                             multiple_intervals_bed, 
                             target_regions_bed, 
                             output_bed = None,
                             padding=600, 
                             logger = logger):
    """
    Extracts and pads segments from genomic intervals in multiple_intervals_bed that correspond to 
    the parts of the interval in single_interval_bed which overlap with target regions. Merges overlapping
    segments after padding and adjusts for strand orientation.

    :param single_interval_bed: Path to the BED file containing the single genomic interval.
    :param multiple_intervals_bed: Path to the BED file containing multiple genomic intervals.
    :param target_regions_bed: Path to the BED file containing target regions.
    :param padding: Number of base pairs to pad to each side of the segment.
    :return: A list of pybedtools.Interval objects representing the corresponding segments.
    """

    # Load the BED files
    if type(single_interval_bed) == str:
        single_interval = pb.BedTool(single_interval_bed)
    else:
        single_interval = single_interval_bed
        
    if type(multiple_intervals_bed) == str:
        multiple_intervals = pb.BedTool(multiple_intervals_bed)
    else:
        multiple_intervals = multiple_intervals_bed
        
    if type(target_regions_bed) == str:
        target_regions = pb.BedTool(target_regions_bed)
    else:
        target_regions = target_regions_bed
    
    # Intersect to get overlapping segments with target regions
    overlapping_segments = single_interval.intersect(target_regions)
    
    # Assuming there's only one interval in single_interval_bed
    main_interval = single_interval[0]
    main_interval_start, main_interval_end, main_interval_strand = main_interval.start, main_interval.end, main_interval.strand

    # Calculate and pad relative coordinates of overlapping segments
    relative_segments = []
    for segment in overlapping_segments:
        # Padding and ensuring the coordinates do not become negative
        start = max(segment.start - padding, 0)
        end = segment.end + padding
        relative_start = start - main_interval_start
        relative_end = end - main_interval_start
        # Store the padded relative coordinates with original start and end for sorting and merging
        relative_segments.append((relative_start, relative_end, start, end))
        
    # logger.info(relative_segments)
    
    # Sort by the original start and merge overlapping segments
    relative_segments.sort(key=lambda x: x[2])  # Sort by original start coordinate
    merged_segments = [relative_segments[0]]
    for current in relative_segments[1:]:
        last = merged_segments[-1]
        # If current segment overlaps with the last, merge them
        if current[2] <= last[3]:
            merged_end = max(last[3], current[3])
            new_segment = (last[0], merged_end - main_interval_start, last[2], merged_end)
            merged_segments[-1] = new_segment  # Update the last segment
        else:
            merged_segments.append(current)
    
    # logger.info(merged_segments)

    # Use relative coordinates to slice out corresponding parts from multiple_intervals_bed
    corresponding_segments = []
    for interval in multiple_intervals:
        interval_strand = interval.strand
        # Extract the relative start and end from the 4th and 5th columns
        rel_start_interval = int(interval[3])  # assuming this is the 4th column in your bed file
        rel_end_interval = int(interval[4])  # assuming this is the 5th column in your bed file
        
        # Calculate the subset of relative segments that overlap with the smaller interval
        overlapping_with_smaller = [
            (max(rel_start, rel_start_interval), min(rel_end, rel_end_interval))
            for rel_start, rel_end, _, _ in merged_segments
            if rel_start < rel_end_interval or rel_end > rel_start_interval
        ]
        
        # logger.info(overlapping_with_smaller)
        
        for rel_start, rel_end in overlapping_with_smaller:
            # Be careful that the relative coordinates are corresponding to the main_interval
            if rel_start >= rel_end:
                continue
            
            if interval_strand == main_interval_strand:
                # Calculate the relative coordinates based on the iterating smaller interval
                rel_small_start = rel_start - rel_start_interval
                rel_small_end = rel_end - rel_start_interval
            else:
                # Calculate the relative coordinates based on the iterating smaller interval
                rel_small_start = rel_end_interval - rel_end
                rel_small_end = rel_end_interval - rel_start
                
            # Calculate absolute coordinates for slicing based on the interval start
            abs_start = interval.start + rel_small_start
            abs_end = interval.start + rel_small_end
            
            # Adjust if coordinates exceed the interval's boundaries
            abs_start = max(interval.start, abs_start)
            abs_end = min(interval.end, abs_end)
            
            # Add to results if the segment is within the interval's bounds
            if abs_start < abs_end:  # Ensure start is less than end
                corresponding_segment = pb.Interval(interval.chrom, abs_start, abs_end, strand=interval_strand, name=interval.name)
                corresponding_segments.append(corresponding_segment)

    if output_bed:
        return_bed = pb.BedTool(corresponding_segments).sort().merge()
        return_bed_size = return_bed.total_coverage()
        if os.path.exists(output_bed) and \
           os.path.getmtime(output_bed) > os.path.getmtime(target_regions_bed) and \
           os.path.getmtime(output_bed) > os.path.getmtime(multiple_intervals_bed) and \
           os.path.getmtime(output_bed) > os.path.getmtime(single_interval_bed):
            ori_output_bed = pb.BedTool(output_bed)
            ori_bed_size = ori_output_bed.total_coverage()
            if ori_bed_size == return_bed_size:
                logger.info(f"Reuse the existing output bed file {output_bed}")
                return output_bed, ori_bed_size
            
        return_bed.saveas(output_bed)
        return output_bed, return_bed_size 
    else:
        return_bed = pb.BedTool(corresponding_segments).sort().merge()
        return_bed_size = return_bed.total_coverage()
        return return_bed, return_bed_size




def imap_prepare_masked_align_region_per_RG(tup_args):
    return prepare_masked_align_region_per_RG(*tup_args)
    


@log_command
def prepare_masked_align_region_per_RG( rg_label: str,
                                        rg_subids: tuple,
                                        tmp_dir: str,
                                        target_region_bed: str,
										whole_region_bed: str,
                                        ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                        logger = logger ):
    # The rg_subid_dict should only contain one PC tag
    rg_label = os.path.basename(whole_region_bed).split("_")[0]
    logger.info(f"The raw bed file for this {rg_label} is {whole_region_bed}")
    whole_region_bedf = pd.read_table(whole_region_bed, header=None, names=["chrom", "start", "end", "col4", "col5", "strand", "tag"])
    fc_region_bedf = whole_region_bedf.loc[whole_region_bedf["tag"].str.contains(r"^FC.*"), :]
    fc_bed_path = prepare_tmp_file(tmp_dir = tmp_dir, suffix=".bed").name
    fc_region_bedf.to_csv(fc_bed_path, sep="\t", header=False, index=False)
    fc_bed = pb.BedTool(fc_bed_path)

    target_region_df = pd.read_table(target_region_bed, header=None).iloc[:, :3]
    target_region_bed_obj = pb.BedTool.from_dataframe(target_region_df)

    targeted_fc_region_bed = prepare_tmp_file(tmp_dir = tmp_dir, suffix=".bed").name
    ref_genome_fai = ref_genome.replace(".fasta", ".fasta.fai")
    target_fc_bed = fc_bed.intersect(target_region_bed_obj).slop(b=500, g=ref_genome_fai).sort().merge()
    target_fc_size = target_fc_bed.total_coverage()
    target_fc_bed.saveas(targeted_fc_region_bed)

    records = []
    for subgroup_id in rg_subids:
        logger.info("Start to fetch the masked align region for NFC regions for PC {} subgroup {}".format(rg_label, subgroup_id))
        nfc_region_bedf = whole_region_bedf.loc[whole_region_bedf["tag"] == f"NFC:{rg_label}_{subgroup_id}", :]
        fc_region_bedf = whole_region_bedf.loc[whole_region_bedf["tag"] == f"FC:{rg_label}_{subgroup_id}", :]
        record = prepare_masked_align_region_per_RG_subgroup(fc_region_bedf,
                                                             nfc_region_bedf,
                                                             whole_region_bed,
                                                             rg_label, 
                                                             subgroup_id, 
                                                             tmp_dir,
                                                             target_region_bed, 
                                                             ref_genome, 
                                                             logger = logger)
        record = record.replace("target_fc_region_bed", f"{targeted_fc_region_bed}")
        record = record.replace("target_fc_region_size", f"{target_fc_size}")
        logger.info(f"The returned record for {rg_label} subgroup {subgroup_id} is {record}")
        records.append(record)

    return tuple(records)



def prepare_masked_align_region_per_RG_subgroup(fc_region_bedf,
                                                nfc_region_bedf,
                                                whole_region_bed: str,
                                                rg_label: str,
                                                subgroup_id: str,
                                                tmp_dir: str,
                                                target_region_bed: str,
                                                ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                                logger = logger):
    # Parse to generate new variables
    
    # fc_region_bedf = whole_region_bedf.loc[whole_region_bedf["tag"] == f"FC:{rg_label}_{subgroup_id}", :]
    # nfc_region_bedf = whole_region_bedf.loc[whole_region_bedf["tag"] == f"NFC:{rg_label}_{subgroup_id}", :]
    assert len(fc_region_bedf) == 1, f"The FC region record {whole_region_bed} for subgroup {subgroup_id} is not unique"
    assert len(nfc_region_bedf) >= 1, f"The NFC region record {whole_region_bed} for subgroup {subgroup_id} does not exist"

    fc_bed_path = prepare_tmp_file(tmp_dir = tmp_dir, suffix=".bed").name
    nfc_bed_path = prepare_tmp_file(tmp_dir = tmp_dir, suffix=".bed").name

    fc_region_bedf.iloc[:, :3].to_csv(fc_bed_path, sep="\t", header=False, index=False)
    nfc_region_bedf.to_csv(nfc_bed_path, sep="\t", header=False, index=False)

    # target_region_bed_obj = pb.BedTool(target_region_bed)
    
    # target_fc_bed = fc_bed.intersect(target_region_bed_obj).slop(b=500, g=ref_genome_fai).sort().merge()
    # target_fc_size = target_fc_bed.sort().total_coverage()
    # target_fc_bed.saveas(targeted_fc_region_bed)
    # logger.info(f"Target FC region bed for {rg_label} subgroup {subgroup_id} is saved at {targeted_fc_region_bed}")

    targeted_nfc_regions_bed = prepare_tmp_file(tmp_dir = tmp_dir, suffix=".bed").name
    targeted_nfc_regions_bed, target_nfc_regions_bed_size = extract_and_pad_segments(fc_bed_path,
                                                                                     nfc_bed_path,
                                                                                     target_region_bed,
                                                                                     targeted_nfc_regions_bed,
                                                                                     logger = logger)

    return f"{rg_label},{subgroup_id},target_fc_region_bed,{targeted_nfc_regions_bed},target_fc_region_size,{target_nfc_regions_bed_size}"

