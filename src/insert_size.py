#!/usr/bin/env python3 

## This script randomly samples reads from a BAM file and computes the insert size distribution.

import pysam
import numpy as np
import random
from src.log import logger

def get_insert_size_distribution(bam_file, num_runs=10, num_samples=2000, max_template_length=8000, sample_prob=0.001):
    """
    Runs insert size estimation multiple times, averages statistics, and handles filtering.

    Args:
        bam_file (str): Path to the BAM file.
        num_runs (int): Number of times to run the estimation.
        num_samples (int): Number of samples per run.
        max_template_length (int): Maximum template length.
        sample_prob (float): Sampling probability.

    Returns:
        tuple: (Average Mean, Average Median, Average Std), or (None, None, None) if no data.
    """
    all_stats = []
    for _ in range(num_runs):
        insert_sizes = []
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam:
                    if (read.flag & 0x2 and not (read.flag & 0x90C) and
                        read.reference_id == read.next_reference_id and
                        abs(read.template_length) <= max_template_length and
                        random.random() < sample_prob):

                        insert_sizes.append(abs(read.template_length))
                        if len(insert_sizes) >= num_samples:
                            break  # Stop when num_samples reached

            if insert_sizes:
                threshold = np.percentile(insert_sizes, 99)
                filtered_sizes = [size for size in insert_sizes if size <= threshold][:num_samples]
                if filtered_sizes:
                  stats = np.array([np.mean(filtered_sizes), np.median(filtered_sizes), np.std(filtered_sizes)])
                  all_stats.append(stats)

        except (FileNotFoundError, ValueError, Exception) as e:
            logger.error(f"Error during run {_ + 1}: {e}")

    return tuple(np.mean(np.stack(all_stats), axis=0)) if all_stats else (None, None, None)

# print(get_insert_size_distribution(bam_file, num_samples=2000))