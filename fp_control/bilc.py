import numba
import numpy as np
import pandas as pd
from numba import prange, types
from highspy import Highs
from collections import defaultdict

from numba_operators import numba_sum
from bam_ncls import overlapping_reads_generator





@numba.njit(types.int32[:](types.int32[:,:], types.float32[:,:], types.int32[:,:]), fastmath=True)
def assemble_consensus(seq_arrays, qual_arrays, reads):
    # Every read in the reads can have different length
    # Find the start and end position of the consensus sequence
    start_pos = reads[:, 0].min()
    end_pos = reads[:, 1].max()

    # logger.info(", ".join([f"{r.reference_start}-{r.reference_end}" for r in reads]))

    # Initialize the consensus sequence and quality scores with zeros
    consensus_seq = np.ones(end_pos - start_pos, dtype=np.int32)
    consensus_qual = np.full(end_pos - start_pos, 0.1, dtype=np.float32)

    for i in prange(len(seq_arrays)):
        '''
        For every iteration,
        seq and qual should have the same length

        Across iterations,
        seq and qual are not bound with the same length
        '''
        seq = seq_arrays[i]
        qual = qual_arrays[i]
        read = reads[i]
        start = read[0]
        # assert seq.size == qual.size, f"Found the hap_vector {list(seq)} ({seq.size}bp) and qual_vector {list(qual)} ({qual.size}bp) have different lengths for read {read} at position {start}"
        # Filter out NaN values
        non_na_values = numba_sum(seq >= -8)

        nona_seq = np.empty(non_na_values, dtype=np.int32)
        nona_qual = np.empty(non_na_values, dtype=np.float32)

        for j in prange(non_na_values):
            nona_seq[j] = seq[j]
            nona_qual[j] = qual[j]

        seq = nona_seq
        qual = nona_qual

        # Calculate the relative start and end positions of the current sequence
        rel_start = start - start_pos
        rel_end = rel_start + non_na_values

        # Create a mask for positions where the current sequence has higher quality
        # assert rel_end - rel_start == len(qual), f"Found the relative start and end positions (determined by {len(seq)}) of the current sequence are not equal to the quality vector size {len(qual)} for read {read} at position {start}"
        # qual here is actually err prob, so smaller the better
        mask = qual <= consensus_qual[rel_start:rel_end]
        mask = mask & (qual < 0.1)

        # Update the consensus sequence and quality scores where the current sequence has higher quality
        consensus_seq[rel_start:rel_end][mask] = seq[mask]
        consensus_qual[rel_start:rel_end][mask] = qual[mask]

    return consensus_seq



def extract_continuous_regions_dict(reads):
    # Sort reads by reference start position
    reads = sorted(reads, key=lambda read: read.reference_start)

    continuous_regions = {}
    current_region_reads = []
    current_start = None
    current_end = None

    for read in reads:
        read_start = read.reference_start
        read_end = read.reference_end

        if current_start is None:
            # Initialize the first region
            current_start = read_start
            current_end = read_end
            current_region_reads.append(read)
        else:
            if read_start <= current_end:
                # Extend the current region
                current_end = max(current_end, read_end)
                current_region_reads.append(read)
            else:
                # Save the current region and start a new one
                continuous_regions[(current_start, current_end)] = current_region_reads
                current_start = read_start
                current_end = read_end
                current_region_reads = [read]

    # Append the last region
    if current_region_reads:
        continuous_regions[(current_start, current_end)] = current_region_reads

    return continuous_regions



def lp_solve_remained_haplotypes(total_record_df,
                                 logger):
    # total_record_df = total_record_df.loc[np.logical_not(total_record_df["extreme_vard"]), :]
    hap_id_coefficients = total_record_df.groupby(["hap_id"])["coefficient"].mean().reset_index()
    hap_id_coefficients = hap_id_coefficients.merge(total_record_df.loc[:, ["hap_id", "ref_genome_similarities"]], how="left", on="hap_id")
    hap_id_coefficients.loc[:, "coefficient"] = hap_id_coefficients["coefficient"] + hap_id_coefficients["ref_genome_similarities"]
    hap_id_coefficients.drop_duplicates(inplace=True)

    logger.info(f"The hap_id_coefficients looks like \n{hap_id_coefficients.to_string(index=False)}\n")
    # rank_3_count = total_record_df[total_record_df["rank"] <= 3].shape[0]
    # rank_1_count = total_record_df[total_record_df["rank"] <= 1].shape[0]
    hap_no = hap_id_coefficients.shape[0]
    '''
    Note that here the hap_id_coefficient table is a dataframe with columns ["hap_id", "coefficient"]
    The row index is used to represent the hap_id in the constraint matrix below, but the the value does not necessarily equal to hap_id
    Therefore, we need two dicts to map from row_index to hap_id and from hap_id to row_index
    '''

    assert len(hap_id_coefficients["hap_id"].unique()) == hap_no, f"The hap_id_coefficients table has {hap_no} rows, but the hap_id column has {len(hap_id_coefficients['hap_id'].unique())} unique values."
    index_to_hapid = dict(zip(range(hap_no), hap_id_coefficients["hap_id"].to_numpy()))
    hapid_to_index = dict(zip(hap_id_coefficients["hap_id"].to_numpy(), range(hap_no)))

    logger.info(f"Now the index_to_hapid dict looks like \n{index_to_hapid}\n")
    logger.info(f"Now the hapid_to_index dict looks like \n{hapid_to_index}\n")

    by_region = total_record_df.groupby(["chrom", "start", "end"])
    inspect_region_no = by_region.ngroups
    # Initialize the HiGHS model
    highs = Highs()

    # Set the output_flag to False to suppress output
    highs.setOptionValue('output_flag', False)
    # Create a constraint matrix
    logger.info(f"The constraint matrix will have a shape of {(inspect_region_no, hap_no)}")

    # Then we need to fill the constraint matrix A properly to correctly specify start indices, row_indices, and values in the future function calls.
    lower_bounds = np.zeros(hap_no, dtype=np.int32)
    upper_bounds = np.ones(hap_no, dtype=np.int32)

    # Set the output_flag to False to suppress output
    highs.setOptionValue('output_flag', True)

    col_costs = - hap_id_coefficients["coefficient"].to_numpy().astype(np.double)

    # Add Vars
    highs.addVars(hap_no, lower_bounds, upper_bounds)

    # Change column costs
    highs.changeColsCost(hap_no, np.arange(hap_no, dtype=np.int32), col_costs)

    # Set variable types to integer
    integrality = np.ones(hap_no, dtype=np.int32)  # 1 for integer variables
    highs.changeColsIntegrality(hap_no, np.arange(hap_no, dtype=np.int32), integrality)

    # Adding rows
    for name, group in by_region:
        included_hapids = group["hap_id"].unique()
        if included_hapids.size <= 2:
            continue
        rank_2_count = group[group["rank"] <= 2].shape[0]
        # Column index extraction
        # hapid_indices = [hap_id_coefficients[hap_id_coefficients["hap_id"] == hapid].index[0] for hapid in included_hapids]
        # hapid_indices = find_indices(hap_id_coefficients["hap_id"].to_numpy(), included_hapids)
        hapid_indices = [hapid_to_index[hapid] for hapid in included_hapids]
        lower_bound = included_hapids.size - rank_2_count
        upper_bound = max(included_hapids.size - 2, lower_bound)
        status = highs.addRow(0,
                              lower_bound,
                              included_hapids.size,
                              np.array(hapid_indices, dtype=np.int32),
                              np.ones(included_hapids.size, dtype=np.double))
        logger.info(f"The addrow status is {status}, the group included hap_ids are {included_hapids}, the corresponding hap_id indices are {hapid_indices} the lower bound is {lower_bound}, the upper bound is {upper_bound}.")

    # Run solver
    highs.run()
    solution = highs.getSolution()
    info = highs.getInfo()
    num_var = highs.getNumCol()
    model_status = highs.getModelStatus()
    basis = highs.getBasis()

    # Access and print the constraint matrix
    lp = highs.getLp()
    logger.info(f"column coefficients: {lp.col_cost_}")
    logger.info(f"Model status = {highs.modelStatusToString(model_status)}")
    logger.debug(f"Optimal objective = {info.objective_function_value}")
    logger.info(f"Iteration count = {info.simplex_iteration_count}")
    logger.info(f"Primal solution status = {highs.solutionStatusToString(info.primal_solution_status)}")

    # Get solution
    solution = np.array(highs.getSolution().col_value)
    # print(f"Solution is {solution}")
    select_hap_indices = np.where(solution == 0)[0]
    drop_hap_indices = np.where(solution == 1)[0]
    select_hap_ids = hap_id_coefficients.iloc[select_hap_indices, 0]
    drop_hap_ids = hap_id_coefficients.iloc[drop_hap_indices, 0]
    return set(select_hap_ids), set(drop_hap_ids), highs.modelStatusToString(model_status)









