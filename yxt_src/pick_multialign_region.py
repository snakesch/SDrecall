import pandas as pd
import numpy as np
import multiprocessing as mp
import subprocess
import logging
import pybedtools as pb
import argparse as ap
import uuid
from io import StringIO
from python_utils import convert_input_value, calculate_inferred_coverage


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def init_logger():
    logger = logging.getLogger(f"Process-{mp.current_process().name}")
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter('%(asctime)s:%(name)s:%(message)s'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def log_decorator(func):
    def wrapper(*args, **kwargs):
        logger = init_logger()
        log_stream = StringIO()
        ch = logging.StreamHandler(log_stream)
        logger.addHandler(ch)
        result = func(*args, logger = logger, **kwargs)

        log_contents = log_stream.getvalue()
        log_stream.close()
        
        print(log_contents, file = sys.stderr)
        return result
    return wrapper



def executeCmd(cmd, stdout_only = True, logger = logger) -> None:
    import subprocess
    from subprocess import PIPE, DEVNULL
    
    if stdout_only:
        result = subprocess.run(cmd, shell=True, capture_output=True)
    else:
        result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=PIPE)
        
    # logger.info(f"Running the following shell command inside python:\n{cmd}\nAnd the output goes like this:\n{result.stdout.decode()}\n\n")
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        if cmd_lst[1][0] != "-":
            raise RuntimeError("Error in " + " ".join(cmd_lst))
        else:
            raise RuntimeError("Error in {}".format(cmd_lst))
    
    return result.stdout.decode()



def stat_depth_over_region(input_bam, 
                           MQ_threshold, 
                           filter_tag,
                           target_region,
                           inferred_coverage = False,
                           return_df = True,
                           target_tag = "FCRs",
                           output_depth = None, 
                           threads=1,
                           logger=logger,
                           genome_file = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.contigsize.genome"):
    
    assert (return_df or output_depth), "Neither specified the output depth file or directly return dataframe in memory. Quit with error"
    
    if inferred_coverage:
        if filter_tag:
            filter_tags = [filter_tag]
        else:
            filter_tags = []

        depth_df = calculate_inferred_coverage( input_bam, 
                                                min_mapq = MQ_threshold, 
                                                filter_tags = filter_tags, 
                                                logger=logger, 
                                                return_df=True, 
                                                genome_file = genome_file,
                                                target_region = target_region,
                                                target_tag = target_tag )

        logger.info("The inferred coverage table has shape of {} and it looks like:\n{}\n\n".format(depth_df.shape,
                                                                                                    depth_df[:5].to_string(index=False)))
            
        if return_df:
            return depth_df
        elif output_depth:
            depth_df.to_csv(output_depth, sep="\t", index=False, header=False)
            return output_depth
    else:
        region_arg = f"sambamba slice -q -L {target_region} {input_bam}" if target_region else ""
        filter_arg = f"-F \"[{filter_tag}] != null\"" if filter_tag else ""
        output_part = f"> {output_depth}" if output_depth else ""

        if region_arg and filter_arg:
            cmd = f"{region_arg} | sambamba view -q -f bam -h -t {threads} {filter_arg} /dev/stdin | samtools depth -J -a -q 13 -s -Q {MQ_threshold} -b {target_region} - {output_part}"
        elif region_arg:
            cmd = f"{region_arg} | samtools depth -J -a -q 13 -s -Q {MQ_threshold} -b {target_region} - {output_part}"
        else:
            cmd = f"samtools depth -J -a -q 13 -s -Q {MQ_threshold} {input_bam} {output_part}"
        
        depth_str = executeCmd(cmd, stdout_only=True, logger=logger)
        if return_df:
            return pd.read_table(StringIO(depth_str), header=None, names=["chrom", "pos", "depth"])
        elif output_depth:
            with open(output_depth, "r") as of:
                of.write(depth_str)
            return output_depth


def main_func_pick_region(input_bam, 
                          output_bed=None,
                          MQ_threshold=50,
                          depth_threshold=10, 
                          minimum_depth=5,
                          target_region=None,
                          target_tag = "FCRs",
                          filter_tag="XA",
                          multialign_frac = 0.7,
                          threads=10,
                          inferred_coverage = True,
                          logger=logger):
    parallel_args = [ (input_bam, 0, None, target_region, inferred_coverage, True, target_tag), 
                      (input_bam, MQ_threshold, None, target_region, inferred_coverage, True, target_tag),
                      (input_bam, 0, filter_tag, target_region, inferred_coverage, True, target_tag) ]

    with mp.Pool(3) as pool:
        raw_depth, high_mq_depth, raw_xa_depth = pool.starmap(stat_depth_over_region, parallel_args)
    
    raw_depth = raw_depth.rename(columns={"depth": "raw_depth"})
    high_mq_depth = high_mq_depth.rename(columns={"depth": "high_MQ_depth"})
    raw_xa_depth = raw_xa_depth.rename(columns={"depth": "raw_XA_depth"})

    logger.info("Before merging, the Raw depth table has shape of {}, the high MQ depth table has shape of {}, while the raw Multialigned depth table has shape of {}".format(raw_depth.shape,
                                                                                                                                                                              high_mq_depth.shape,
                                                                                                                                                                              raw_xa_depth.shape))
    merged_depth = raw_depth.merge(high_mq_depth, how="left", on=["chrom", "pos"]).merge(raw_xa_depth, how="left", on=["chrom", "pos"]).drop_duplicates()
    merged_depth = merged_depth.fillna(0)
    logger.info("After merging, the merged table has shape of {} and it looks like :\n{}\n\n".format(merged_depth.shape,
                                                                                                     merged_depth[:5].to_string(index=False)))

    # Now we filter the bases that do not fulfill the standards below: (AND logic)
    # 1. Has at least 5 reads covered. (no MQ considered)
    # 2. Has 60% of the overlapping reads with XA tag
    # 3. Has no more than 10 MQ >= 50 reads covered

    min_depth = merged_depth["raw_depth"] >= minimum_depth
    most_xa = (merged_depth["raw_XA_depth"]/merged_depth["raw_depth"]) >= multialign_frac
    not_enough_evidence = merged_depth["high_MQ_depth"] <= depth_threshold

    merged_depth = merged_depth.loc[min_depth & most_xa & not_enough_evidence, :].drop_duplicates()
    logger.info("After the filtration, the merged depth has covered {} bp.".format(merged_depth.shape[0]))
    merged_depth["start"] = merged_depth.loc[:, "pos"] - 1
    merged_depth["end"] = merged_depth.loc[:, "pos"]
    bed_str = "\n".join(["\t".join([str(x) for x in row.values]) for ind, row in merged_depth.loc[:, ["chrom", "start", "end"]].drop_duplicates().iterrows()])
    bed_obj = pb.BedTool(bed_str, from_string=True).sort().merge()
    if output_bed:
        bed_obj.saveas(output_bed)
        return output_bed
    else:
        return bed_obj



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--function", type=str, help="The function name", required=True)
    parser.add_argument("-a", "--arguments", type=str, help="The function's input arguments, delimited by semi-colon ;", required=False, default=None)
    parser.add_argument("-k", "--key_arguments", type=str, help="Keyword arguments for the function, delimited by semi-colon ;", required=False, default=None)
    
    args = parser.parse_args()
    try:
        fargs = [ convert_input_value(a) for a in args.arguments.split(";") ] if type(args.arguments) == str else []
        fkwargs = { t.split("=")[0]: convert_input_value(t.split("=")[1]) for t in args.key_arguments.split(";") } if type(args.key_arguments) == str else {}
        logger.info("Running function: {}, input args are {}, input kwargs are {}".format(args.function, fargs, fkwargs))
    except Exception as e:
        logger.error("Input argument does not meet the expected format, encounter Parsing error {}, Let's check the input:\n-f {}, -a {}, -k {}".format(
            e,
            args.function,
            args.arguments,
            args.key_arguments
        ))
        raise e
    globals()[args.function](*fargs, **fkwargs)
