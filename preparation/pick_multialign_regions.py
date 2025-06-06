import multiprocessing as mp
from pybedtools import BedTool
from src.log import logger
from preparation.inferred_depths import calculate_inferred_coverage


def pick_multialigned_regions(input_bam, 
                              output_bed = None,
                              MQ_threshold=41,
                              high_quality_depth=10, 
                              minimum_depth=3,
                              multialign_frac = 0.5,
                              target_region=None,
                              target_tag = "FCRs",
                              threads=4,
                              genome_file=""):
    
    ## genome_file should be the FAI index of the reference genome

    ## no need to specify filter tags since we are using XA and XS separately
    parallel_args = [ (input_bam, 0, [], target_region, target_tag, genome_file), 
                      (input_bam, MQ_threshold, [], target_region,  target_tag, genome_file),
                      (input_bam, 0, ["XA"], target_region, target_tag, genome_file),
                      (input_bam, 0, ["XS"], target_region, target_tag, genome_file) ]

    with mp.Pool(min(4, threads)) as pool:
        raw_depth, high_mq_depth, raw_xa_depth, raw_xs_depth = pool.starmap(calculate_inferred_coverage, parallel_args)

    raw_depth = raw_depth.rename(columns={"depth": "raw_depth"})
    high_mq_depth = high_mq_depth.rename(columns={"depth": "high_MQ_depth"})
    raw_xa_depth = raw_xa_depth.rename(columns={"depth": "raw_XA_depth"})
    raw_xs_depth = raw_xs_depth.rename(columns={"depth": "raw_XS_depth"})

    merged_depth = raw_depth.merge(high_mq_depth, how="left", on=["chrom", "pos"]).merge(raw_xa_depth, how="left", on=["chrom", "pos"]).merge(raw_xs_depth, how="left", on=["chrom", "pos"]).drop_duplicates()
    merged_depth = merged_depth.fillna(0)
    logger.debug("After merging, the table has shape {}".format(merged_depth.shape))

    # Now we filter the bases that do not fulfill the standards below: (AND logic)
    # 1. Has at least 5 reads covered. (no MQ considered)
    # 2. Has 60% of the overlapping reads with XA tag
    # 3. Has no more than 10 MQ >= 50 reads covered

    min_depth = merged_depth["raw_depth"] >= minimum_depth
    most_xa = (merged_depth["raw_XA_depth"]/merged_depth["raw_depth"]) >= multialign_frac
    most_xs = (merged_depth["raw_XS_depth"]/merged_depth["raw_depth"]) >= multialign_frac
    not_enough_evidence = merged_depth["high_MQ_depth"] <= high_quality_depth

    merged_depth = merged_depth.loc[min_depth & (most_xa | most_xs) & not_enough_evidence, :].drop_duplicates()
    logger.info("After filtration, the multialign BED covers {}bp.".format(merged_depth.shape[0]))
    merged_depth["start"] = merged_depth.loc[:, "pos"] - 1
    merged_depth["end"] = merged_depth.loc[:, "pos"]

    bed_obj = BedTool.from_dataframe(merged_depth.loc[:, ["chrom", "start", "end"]]).sort().merge()
    
    if target_region:
        target_bed = BedTool(target_region)
        bed_obj = bed_obj.intersect(target_bed).sort().merge()
    
    if output_bed:
        bed_obj.saveas(output_bed)
        return output_bed
    
    return bed_obj


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Pick multialigned regions from BAM file")
    parser.add_argument("--input_bam", type=str, required=True, help="Input BAM file")
    parser.add_argument("--output_bed", type=str, required=True, help="Output BED file")
    parser.add_argument("--MQ_threshold", type=int, default=41, help="MQ threshold for high-quality reads")
    parser.add_argument("--high_quality_depth", type=int, default=10, help="High-quality depth threshold")
    parser.add_argument("--minimum_depth", type=int, default=3, help="Minimum depth threshold")
    parser.add_argument("--multialign_frac", type=float, default=0.5, help="Multialign fraction threshold")
    parser.add_argument("--target_region", type=str, help="Target region BED file")
    parser.add_argument("--target_tag", type=str, default="exome", help="Target tag")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--genome_file", type=str, required=True, help="Genome file")

    args = parser.parse_args()

    output_bed = pick_multialigned_regions( args.input_bam, 
                                            args.output_bed, 
                                            args.MQ_threshold, 
                                            args.high_quality_depth, 
                                            args.minimum_depth, 
                                            args.multialign_frac, 
                                            args.target_region, 
                                            args.target_tag, 
                                            args.threads, 
                                            args.genome_file )
