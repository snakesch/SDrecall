import os
import uuid
import sys
import logging
import multiprocessing as mp

import pandas as pd
import pysam
import pybedtools as pb
from src.utils import executeCmd

logger = logging.getLogger("SDrecall")

def get_read_metadata(read):
    return read.reference_name, read.reference_start, read.reference_end, read.query_name

def common_read_filter(read, min_mapq):
    read_based = not (read.is_duplicate or read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail)
    good_quality = read.mapping_quality >= min_mapq
    return read_based and good_quality

def get_filtered_read_qname(read, min_mapq: int, filter_tags: list[str]) -> str:
    '''
    Selects reads mapped with ambiguity for SDrecall

    Returns:
    bed_feature (str): chrom, start, end and query name in BED format
    '''
    passed_filter = all([read.get_tag("AS") - read.get_tag("XS") <= 5 if tag == "XS" else read.has_tag(tag) for tag in filter_tags]) if len(filter_tags) > 1 else True
    if not common_read_filter(read, min_mapq) or not passed_filter:
        return None
    if (not read.mate_is_unmapped) and read.is_proper_pair:
        # Ensure we process each pair only once
        if read.is_read1:
            chrom, start, end, read_query_name = get_read_metadata(read)
            if read.next_reference_name == chrom:
                bed_feature = f"{chrom}\t{min(start, read.next_reference_start)}\t{max(end, read.next_reference_start + 150)}\t{read_query_name}\n"
            else:
                bed_feature = f"{chrom}\t{start}\t{end}\t{read_query_name}\n"
        else:
            return None
    elif (not read.is_proper_pair) and read.is_paired:
        chrom, start, end, read_query_name = get_read_metadata(read)
        bed_feature = f"{chrom}\t{start}\t{end}\t{read_query_name}\n"
    elif not read.is_paired:
        chrom, start, end, read_query_name = get_read_metadata(read)
        bed_feature = f"{chrom}\t{start}\t{end}\t{read_query_name}\n"
    
    return bed_feature

def calculate_inferred_coverage(bam_file,
                              min_mapq,
                              filter_tags: str,
                              target_region,
                              genome_file = "",
                              output_bed = "",
                              logger = logger):
    '''
    Multiple filter tags are supported as comma-delimited strings. E.g. XA,XS
    '''
    logger.debug(f"Start calculating coverage with MAPQ = {min_mapq}; TAGS = {filter_tags}")

    filter_tags = filter_tags.split(",") if filter_tags != "" else []

    target_tag = os.path.basename(target_region).split(".")[0] if target_region else None

    output_prefix = os.path.join(os.path.dirname(output_bed), os.path.basename(bam_file).replace(".bam", ""))

    if len(filter_tags) > 0:
        output_tsv = output_prefix + f".infercov.minMQ_{min_mapq}.{','.join(filter_tags)}.tsv"
    else:
        output_tsv = output_prefix + f".infercov.minMQ_{min_mapq}.tsv"
    
    ## Avoid concurrent writing
    output_bed = output_tsv.replace(".tsv", ".bed")
    tmp_output_bed = output_bed.replace(".bed", f".{uuid.uuid4()}.bed")
    
    # Adjacent intervals might be counted twice, duplicates have to be dropped
    if os.path.exists(output_bed) and (os.path.getmtime(output_bed) > os.path.getmtime(bam_file)):
        logger.info("The output bed file {} is newer than the input bam file {}, skipping coverage calculation".format(output_bed, bam_file))
    else:
        bed_obj = pb.BedTool(target_region).sort().merge() if target_region else None

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            with open(tmp_output_bed, "w") as ob:
                if bed_obj:
                    for interval in bed_obj:
                        bam_iter = bam.fetch(interval.chrom, interval.start, interval.end)
                        for read in bam_iter:
                            bed_feature = get_filtered_read_qname(read, min_mapq, filter_tags)
                            if bed_feature is not None:
                                ob.write(bed_feature)
                else:
                    bam_iter = bam.fetch()
                    for read in bam_iter:
                        bed_feature = get_filtered_read_qname(read, min_mapq, filter_tags)
                        if bed_feature is not None:
                            ob.write(bed_feature)

        # Sometimes the adjacent interval might fetch the same pair twice. To prevent double counting in the depth calculation. We need to drop the duplicates
        executeCmd(f"echo \"First 5 qualifying reads are:\" && head -n5 {tmp_output_bed}")
        executeCmd(f"sort -uV {tmp_output_bed} > {output_bed} && rm {tmp_output_bed}")
    
    logger.debug(f"Start calculating bedtools genomic coverage with MAPQ = {min_mapq}; TAGS = {filter_tags}")
    cmd = f"bedtools genomecov -bg -i {output_bed} -g {genome_file} | \
            sort -k 1,1 -k 2,2n - | \
            mawk -F '\\t' '{{ for (i=$2+1;i<=$3;i++) {{printf \"%s\\t%i\\t%i\\n\", $1, i, $4;}} }}' - > {output_tsv}"
    executeCmd(cmd)
    if os.path.exists(output_tsv):
        cmd = f"bgzip -f {output_tsv} && tabix -f -s 1 -b 2 -e 2 {output_tsv}.gz"
        executeCmd(cmd)
    else:
        logger.error(f"Unable to compute genomic coverage. Probably due to corrupted BED files. ")
        sys.exit(1)

    logger.debug(f"Finished calculating bedtools genomic coverage with MAPQ = {min_mapq}; TAGS = {filter_tags}")

    depth_df = pd.read_table(output_tsv + '.gz', compression='gzip', header=None, names=['chrom', 'pos', 'depth'], low_memory=False)
    logger.debug("First 5 rows of coverage table {} are: \n{}".format(output_tsv + '.gz', depth_df[:5].to_string(index=False)))
    logger.debug("Inferred coverage of {} bases:\n{}".format(depth_df.shape[0], depth_df.head(5).to_string(index=False)))
    return depth_df

def pick_multialigned_regions(input_bam: str,
                output_bed: str,
                MQ_threshold=41,
                high_quality_depth=10, 
                minimum_depth=3,
                target_region=None,
                multialign_frac = 0.7,
                threads=10,
                genome_file="",
                logger=logger):
    '''
    Consolidate read depths computed with 4 different settings and derive a set of recall regions
    '''
    # Compute the inferred coverage of individual bases within target region by different qualities
    parallel_args = [ (input_bam, 0,            "",     target_region, genome_file, output_bed, logger), 
                      (input_bam, MQ_threshold, "",     target_region, genome_file, output_bed, logger),
                      (input_bam, 0,            "XA",   target_region, genome_file, output_bed, logger),
                      (input_bam, 0,            "XS",   target_region, genome_file, output_bed, logger) ]

    with mp.Pool(min(1, threads)) as pool:
        raw_depth, high_mq_depth, raw_xa_depth, raw_xs_depth = pool.starmap(calculate_inferred_coverage, parallel_args)
    
    raw_depth = raw_depth.rename(columns={"depth": "raw_depth"})
    high_mq_depth = high_mq_depth.rename(columns={"depth": "high_MQ_depth"})
    raw_xa_depth = raw_xa_depth.rename(columns={"depth": "raw_XA_depth"})
    raw_xs_depth = raw_xs_depth.rename(columns={"depth": "raw_XS_depth"})

    logger.info(f"Base coverage:\nAll bases: {raw_depth.shape[0]}\nBases covered by high quality reads: {high_mq_depth.shape[0]}\nBases covered by XA reads: {raw_xa_depth.shape[0]}")
    merged_depth = raw_depth.merge(high_mq_depth, how="left", on=["chrom", "pos"]).merge(raw_xa_depth, how="left", on=["chrom", "pos"]).merge(raw_xs_depth, how="left", on=["chrom", "pos"]).drop_duplicates()
    merged_depth = merged_depth.fillna(0)
    logger.debug("After consolidating multiple inferred depths, {} bases remained for filtration:\n{}\n".format(merged_depth.shape[0],
                                                                                                     merged_depth.head(5).to_string(index=False)))

    # - Identify candidate regions for short variant recall - #
    # Bases not satisfying any of the below criteria are excluded
    # 1. minimum coverage of 5 (irrespective of MAPQ)
    # 2. 60% of the overlapping reads tagged with XA or XS
    # 3. covered by no more than ten MQ >= 50 reads

    min_depth = merged_depth["raw_depth"] >= minimum_depth
    most_xa = (merged_depth["raw_XA_depth"]/merged_depth["raw_depth"]) >= multialign_frac
    most_xs = (merged_depth["raw_XS_depth"]/merged_depth["raw_depth"]) >= multialign_frac
    not_enough_evidence = merged_depth["high_MQ_depth"] <= high_quality_depth

    merged_depth = merged_depth.loc[min_depth & (most_xa | most_xs) & not_enough_evidence, :].drop_duplicates()
    logger.debug("After filtration, {} bp remained for recall.".format(merged_depth.shape[0]))
    merged_depth["start"] = merged_depth.loc[:, "pos"] - 1
    merged_depth["end"] = merged_depth.loc[:, "pos"]
    bed_str = "\n".join(["\t".join([str(x) for x in row.values]) for ind, row in merged_depth.loc[:, ["chrom", "start", "end"]].drop_duplicates().iterrows()])
    bed_obj = pb.BedTool(bed_str, from_string=True).sort().merge()

    if target_region:
        target_bed = pb.BedTool(target_region)
        bed_obj = bed_obj.intersect(target_bed).sort().merge()
    
    bed_obj.saveas(output_bed)

    return output_bed

