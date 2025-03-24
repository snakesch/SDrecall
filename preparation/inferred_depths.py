import os

import pysam
import pandas as pd
from pybedtools import BedTool

from src.utils import logger


def filter_and_process_read(read, min_mapq, filter_tags, filter_logic):
    """Filters a single read and returns processed data or None if filtered."""
    if filter_tags:
        if filter_logic not in ["or", "and", "not"]:
            raise ValueError("filter_logic must be 'and', 'or', or 'not'")
        filter_func = {"or": any, "and": all, "not": lambda x: not any(x)}[filter_logic]

        filter_pass = filter_func([
            read.get_tag("AS") - read.get_tag("XS") <= 5 if tag == "XS" and read.has_tag("AS") and read.has_tag("XS")
            else read.has_tag(tag) for tag in filter_tags
        ])
    else:
        filter_pass = True

    if (read.mapping_quality >= min_mapq and not read.is_unmapped and not read.is_duplicate and
            not read.is_secondary and not read.is_supplementary and not read.is_qcfail and filter_pass):

        chrom, start, end, read_id = read.reference_name, read.reference_start, read.reference_end, f"{read.query_name}:{read.flag}"
        return [chrom, start, end, read_id]
    return None


def create_genome_dict(fai_file):
    """Creates a genome dictionary from a FASTA index file."""
    genome_dict = {}
    with open(fai_file, 'r') as f:
        for line in f:
            chrom, length, *rest = line.strip().split('\t')
            genome_dict[chrom] = (0, int(length))
    return genome_dict


def calculate_coverage(bam_file, min_mapq=10, filter_tags: list[str] = None,
                       target_region="", target_tag="FCRs", genome_file="",
                       filter_logic="and"):
    """Calculates coverage from a BAM file."""
    region_arg = f"samtools view -P -L {target_region} {bam_file}" if target_region else ""
    filter_arg = f"-F \"[{filter_tags[0]}] != null\"" if filter_tags else ""
    output_part = f"> {output_depth}" if output_depth else ""

    if region_arg and filter_arg:
        cmd = f"{region_arg} | sambamba view -q -f bam -h -t {threads} {filter_arg} /dev/stdin | samtools depth -J -a -q 13 -s -Q {MQ_threshold} -b {target_region} - {output_part}"
    elif region_arg:
        cmd = f"{region_arg} | samtools depth -J -a -q 13 -s -Q {MQ_threshold} -b {target_region} - {output_part}"
    else:
        cmd = f"samtools depth -J -a -q 13 -s -Q {MQ_threshold} {bam_file} {output_part}"
    
    depth_str = executeCmd(cmd, stdout_only=True, logger=logger)
    if return_df:
        return pd.read_table(StringIO(depth_str), header=None, names=["chrom", "pos", "depth"])
    elif output_depth:
        with open(output_depth, "r") as of:
            of.write(depth_str)
        return output_depth
    
    

def calculate_inferred_coverage(bam_file, min_mapq=10, filter_tags: list[str] = None,
                                target_region="", target_tag="FCRs", genome_file="",
                                filter_logic="and"):
    """Calculates inferred coverage from a BAM file."""

    logger.info(f"Running (MQ {min_mapq}; filter {filter_tags})")

    if not os.path.isfile(genome_file):
        raise FileNotFoundError(f"Genome index not found for reference {genome_file}")

    name_filter = f".{filter_logic}_" + "_".join(filter_tags) if filter_tags else ""
    target_suffix = f".{target_tag or (os.path.basename(target_region).split('.')[0] if target_region else '')}"
    output_tsv_base = bam_file.replace(".bam", f".infercov.minMQ_{min_mapq}{name_filter}{target_suffix}")
    output_tsv = output_tsv_base + ".tsv"

    self_script = os.path.abspath(__file__)

    if os.path.exists(output_tsv) and \
       os.path.getmtime(output_tsv) > os.path.getmtime(bam_file) and \
       os.path.getmtime(output_tsv) > os.path.getmtime(self_script):
        logger.info(f"Inferred coverage file {output_tsv} already exists. Skipping calculation.")
        return pd.read_csv(output_tsv, sep="\t", header=None, names=["chrom", "pos", "depth"])

    reads_list = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        if target_region:
            bed_obj = BedTool(target_region).sort().merge()
            df = bed_obj.to_dataframe()
            if not df.empty:
                bed_str = df.iloc[:, :3].apply(lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1).to_list()
                logger.debug(f"Target region(s): {bed_str[:5]} ... (first 5 shown)")
                for region_str in bed_str:
                    for read in bam.fetch(region=region_str):
                        processed_read = filter_and_process_read(read, min_mapq, filter_tags, filter_logic)
                        if processed_read:
                            reads_list.append(processed_read)
            else:
                logger.warning("Target BED file is empty after merging.")

        else:
            for read in bam.fetch():
                processed_read = filter_and_process_read(read, min_mapq, filter_tags, filter_logic)
                if processed_read:
                    reads_list.append(processed_read)


    if not reads_list:
        logger.warning("No reads passed filtering. Creating empty output files.")
        coverage_df = pd.DataFrame(columns=['chrom', 'pos', 'depth'])
    else:
        bed_df = pd.DataFrame(reads_list, columns=['chrom', 'start', 'end', 'read_id']).dropna().astype({'start': int, 'end': int})
        logger.debug(f"BED DataFrame (first 5 rows):\n{bed_df.head(5).to_string(index=False)}")

        # Create the genome dictionary from the .fai file
        genome_dict = create_genome_dict(genome_file)
        
        coverages = []
        for interval in BedTool.from_dataframe(bed_df).set_chromsizes(genome_dict).genome_coverage(bg=True):
            for i in range(interval.start + 1, interval.end + 1):
                coverages.append((interval.chrom, i, int(interval[3])))
        coverage_df = pd.DataFrame(coverages, columns=['chrom', 'pos', 'depth'])
    
    logger.debug(coverage_df.head(5))
    coverage_df.to_csv(output_tsv, sep="\t", header=False, index=False)
    logger.debug(f"Coverage df written to {output_tsv}")

    logger.debug(f"Coverage table (first 5 rows):\n{coverage_df.head(5).to_string(index=False)}")
    logger.debug(f"Shape of table [MQ {min_mapq}; filter {filter_tags}]: {coverage_df.shape}")

    return coverage_df