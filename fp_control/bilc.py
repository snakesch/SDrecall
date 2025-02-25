#! /usr/bin/env python3

import logging
import numba
import numpy as np
from highspy import Highs

from src.log import logger

"""
This module implements a Binary Integer Linear Constraint (BILC) Programming approach to select
optimal haplotypes from a set of candidate haplotypes.

It uses the HiGHS solver (open source) to solve the BILC problem.

The main function, lp_solve_remained_haplotypes, takes a DataFrame of haplotype
records and constructs a BILC problem to select the best subset of haplotypes.
The objective is to maximize the sum of haplotype coefficients while satisfying
constraints on the number of selected haplotypes per genomic region.

Key Features:
- Constructs a BILC problem from haplotype records
- Uses HiGHS solver for efficient solving of large-scale BILC problems
- Handles constraints on haplotype selection per genomic region
- Every haplotype is treated as a binary variable, and an coefficient is assigned to it which represents the general likelihood that the haplotype is misaligned

Functions:
- find_indices: Helper function to find indices (currently deprecated)
- lp_solve_remained_haplotypes: Main function to solve the haplotype selection problem

Dependencies:
- numpy
- pandas
- highspy (HiGHS solver Python interface)
- numba (for potential performance optimizations)

Usage:
This script only contains the HiGHs solver part used in a bigger workflow to identify misaligned haplotypes.

Note:
This implementation assumes that the input DataFrame has specific columns including
'hap_id', 'coefficient', 'ref_genome_similarities', 'chrom', 'start', 'end', and 'rank'.
Ensure that the input data is properly formatted before using this module.
"""



@numba.njit
def find_indices(hap_ids, included_hapids):
    '''
    Deprecated for now
    '''
    hapid_indices = np.empty(len(included_hapids), dtype=np.int32)
    for i in range(len(included_hapids)):
        hapid = included_hapids[i]
        indices = np.where(hap_ids == hapid)[0]
        if len(indices) > 0:
            hapid_indices[i] = indices[0]
    return hapid_indices



def lp_solve_remained_haplotypes(total_record_df,
                                 logger = logger):
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
    # num_var = highs.getNumCol()
    model_status = highs.getModelStatus()
    # basis = highs.getBasis()

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









