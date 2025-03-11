import os

import gc
import re
import sys

import pandas as pd
import numpy as np
import multiprocessing as mp
ctx = mp.get_context("spawn")

from io import StringIO
from datetime import datetime
from src.utils import executeCmd, configure_parallelism, merge_bams
from src.log import logger
from src.const import shell_utils
from src.suppress_warning import *

from realign_recall.slice_bam_by_cov import split_bam_by_cov
from fp_control.realign_filter_per_cov import imap_filter_out
from realign_recall.cal_edge_NM_values import calculate_NM_distribution_poisson

def gt_filter_init(threads, cache_dir=None):
    # Set thread count for all common numerical libraries
    # Ensure threading layer is configured before any numba call
    os.environ["NUMBA_THREADING_LAYER"] = "omp"

    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["TBB_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)  # Add Intel MKL
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)  # Add OpenBLAS
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)  # Add Apple's Accelerate
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)  # Add numexpr

    # Configure cache directory if provided
    if cache_dir:
        os.environ["NUMBA_CACHE_DIR"] = cache_dir

    # Set NumPy threading (if using recent NumPy versions)
    try:
        np.config.threading_layer = 'threaded'  # Use Python threads
    except:
        pass
        



def eliminate_misalignments(input_bam,
                            output_bam,
                            intrinsic_bam,
                            original_bam,
                            ref_genome,
                            target_regions,
                            avg_frag_size = 400,
                            threads = 12,
                            numba_threads = 4,
                            conf_level = 0.01,
                            stat_sample_size = 1000000,
                            mapq_cutoff = 20,
                            basequal_median_cutoff = 15,
                            edge_weight_cutoff = 0.201,
                            cache_dir = "/tmp",
                            logger = logger):

    self_path = os.path.abspath(__file__)
    cmd = f"bash {shell_utils} quick_check_bam_validity {output_bam} && \
            bash {shell_utils} quick_check_bam_validity {input_bam} && \
            [[ {output_bam} -nt {input_bam} ]] && \
            [[ {output_bam} -nt {self_path} ]]"

    try:
        executeCmd(cmd, logger=logger)
    except RuntimeError:
        logger.info(f"Start to post process the bam file {input_bam}")
        splitted_bams, splitted_beds = split_bam_by_cov(input_bam, 
                                                        target_bed = target_regions,
                                                        delimiter_size=int(np.ceil(avg_frag_size * 1.5)), 
                                                        logger = logger, 
                                                        threads = threads,
                                                        ref_genome = ref_genome,
                                                        tmp_dir = cache_dir)
        splitted_intrin_bams, _ = split_bam_by_cov( intrinsic_bam, 
                                                    target_bed = target_regions,
                                                    beds = splitted_beds, 
                                                    logger = logger, 
                                                    threads = threads,
                                                    ref_genome = ref_genome,
                                                    tmp_dir = cache_dir)
        
        # There is a possiblity that the splitted_intrin_bams might be empty, so we need to check it
        # Sort the list for load balancing
        tup_bams = list(zip(splitted_bams, splitted_intrin_bams, splitted_beds))
        assert len(splitted_intrin_bams) == len(splitted_bams), f"The splitted intrin BAMs number is not equal to splitted BAMs number"
        tb_removed_indices = set()
        for i in range(len(splitted_bams)):
            if splitted_bams[i] == "NaN":
                logger.warning(f"The BAM file {splitted_bams[i]} for {splitted_beds[i]} is not valid, so we need to remove it from the list")
                tb_removed_indices.add(i)
            elif splitted_intrin_bams[i] == "NaN":
                logger.warning(f"The intrinsic BAM file {splitted_intrin_bams[i]} for {splitted_beds[i]} is not valid, so we need to remove it from the list")
                tb_removed_indices.add(i)

        tup_bams = [tup_bams[i] for i in range(len(tup_bams)) if i not in tb_removed_indices]
        tup_bams = sorted(tup_bams, key=lambda x: os.path.getsize(x[0]), reverse=True)
        raw_bams = [cb[0] for cb in tup_bams]
        intrinsic_bams = [cb[1] for cb in tup_bams]
        raw_bam_regions = [cb[2] for cb in tup_bams]
        clean_bams = [ re.sub(r"\.raw\.", ".", cb) for cb in raw_bams]
        
        # Filter out the misaligned_reads per chromosome
        job_num, _ = configure_parallelism(threads, numba_threads)
        nm_cutoff, _ = calculate_NM_distribution_poisson(original_bam, conf_level, stat_sample_size, logger=logger)

        # Create a dedicated log directory
        log_dir = os.path.dirname(output_bam)
        os.makedirs(log_dir, exist_ok=True)

        with ctx.Pool(job_num, initializer=gt_filter_init, initargs=(numba_threads, cache_dir)) as pool:
            # Pass the log_dir parameter to the function
            result_records = pool.starmap(
                imap_filter_out, 
                [(
                    (raw_bam, 
                     clean_bam, 
                     intrinsic_bam, 
                     raw_bam_region,
                     nm_cutoff,
                     mapq_cutoff,
                     basequal_median_cutoff,
                     edge_weight_cutoff,
                     numba_threads,
                     cache_dir,
                     i), 
                    log_dir
                ) for i, (raw_bam, clean_bam, intrinsic_bam, raw_bam_region) in enumerate(zip(raw_bams, clean_bams, intrinsic_bams, raw_bam_regions))]
            )
            
            result_record_strs = []
            i = 0
            for result in result_records:
                i += 1
                print(f"\n************************************{i}_subprocess_start_for_filtering************************************\n", file=sys.stderr)
                
                result_parts = result.split(",")
                raw_bam = result_parts[0]
                
                # Extract the log file path (now at position 2)
                log_file = result_parts[2] if len(result_parts) >= 3 else "No log file"
                
                # For result processing, we need to add placeholders to match the expected 5 fields
                # The original expected: raw_bam, masked_bam, clean_masked_bam, noise_masked_bam, masked_vcf
                if len(result_parts) == 3 and result_parts[1] != "NaN":  # Successful result: raw_bam, output_bam, log_file
                    standard_result = f"{result_parts[0]},{result_parts[1]}"
                else:  # Error result: raw_bam, NaN, log_file
                    standard_result = f"{result_parts[0]},NaN"
                    
                result_record_strs.append(standard_result)
                
                # Print a reference to the log file instead of all log content
                print(f"Subprocess {i} for {raw_bam} complete. Full log available at: {log_file}", file=sys.stderr)
                
                # Display brief status based on result
                if len(result_parts) >= 2 and result_parts[1] != "NaN":
                    print(f"  Status: SUCCESS - Created output: {result_parts[1]}", file=sys.stderr)
                else:
                    # For errors, show the last few lines of the log file to help with debugging
                    try:
                        with open(log_file, 'r') as f:
                            log_tail = "".join(f.readlines()[-5:])  # Last 5 lines
                        print(f"  Status: ERROR - Check log file for details", file=sys.stderr)
                        print(f"  Error summary:\n{log_tail}", file=sys.stderr)
                    except:
                        print(f"  Status: ERROR - Could not read log file", file=sys.stderr)
                    
                print(f"{datetime.now()}: ************************************{i}_subprocess_end_for_filtering_{raw_bam}************************************", file=sys.stderr)

            pool.close()
            gc.collect()
            
            # Continue with existing result processing...
            logger.info("Start to compose the parallel returned results into a dataframe")
            post_result_tab = os.path.join(os.path.dirname(raw_bams[0]), "post_result.tsv")
            post_result_df = pd.read_csv(StringIO("\n".join(result_record_strs)), header=None, names=["raw_masked_bam", "masked_bam"], na_values="NaN")
            post_result_df.to_csv(post_result_tab, sep="\t", index=False)
            
            # Log summary statistics about the results
            success_count = len(post_result_df.dropna(subset=["masked_bam"]))
            total_count = len(post_result_df)
            logger.info(f"Processing completed: {success_count}/{total_count} regions successfully processed")
            
            # Log output file location but not the entire content
            logger.info(f"Detailed results saved to {post_result_tab}")
            
            # Show a summary table with just the counts by status
            summary = post_result_df.notna().sum()
            logger.info(f"Result summary:\n{summary}")
            
            # Now we merge the masked bams and vcfs across different chromosomes
            merge_bams(post_result_df.loc[:, "masked_bam"].dropna().drop_duplicates().tolist(), 
                        output_bam, 
                        ref_fasta = ref_genome, 
                        threads=threads - 1, 
                        tmp_dir = cache_dir, 
                        logger=logger)
            merge_bams(post_result_df.loc[:, "raw_masked_bam"].dropna().drop_duplicates().tolist(), 
                        input_bam, 
                        ref_fasta = ref_genome, 
                        threads=threads - 1, 
                        tmp_dir = cache_dir, 
                        logger=logger)
            return input_bam, output_bam
            
    else:
        logger.info(f"The post processed masked bam and vcf files are already up-to-date. No need to run the post processing again.")
        return input_bam, output_bam


if __name__ == "__main__":
    import argparse as ap
    parser = ap.ArgumentParser(description="Eliminate the misalignments from the input bam file")
    parser.add_argument("--input_bam", type=str, required=True, help="The input bam file")
    parser.add_argument("--output_bam", type=str, required=True, help="The output bam file")
    parser.add_argument("--intrinsic_bam", type=str, required=True, help="The intrinsic bam file")
    parser.add_argument("--ref_genome", type=str, required=True, help="The reference genome file")
    parser.add_argument("--avg_frag_size", type=int, required=True, help="The average fragment size")
    parser.add_argument("--threads", type=int, required=True, help="The number of threads")
    parser.add_argument("--numba_threads", type=int, required=True, help="The number of threads for numba")
    parser.add_argument("--logger", type=str, required=True, help="The logger")
    args = parser.parse_args()

    eliminate_misalignments(args.input_bam,
                            args.output_bam,
                            args.intrinsic_bam,
                            args.ref_genome,
                            args.avg_frag_size,
                            args.threads,
                            args.numba_threads,
                            args.logger)


