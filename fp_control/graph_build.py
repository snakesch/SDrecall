import numba
import graph_tool.all as gt
import numpy as np
import pandas as pd

from collections import defaultdict
from numba import types
from numba.typed import Dict

from src.utils import executeCmd
from src.suppress_warning import *
from fp_control.numba_operators import any_false_numba, numba_sum
from fp_control.bam_ncls import overlap_qname_idx_iterator
from fp_control.pairwise_read_inspection import determine_same_haplotype
from src.log import logger


def stat_ad_to_dict(bam_file, empty_dict, reference_genome, base_dict = {"A": np.int8(0), "T": np.int8(1), "C": np.int8(2), "G": np.int8(3), "N": np.int8(4)}, logger = logger):
    # Build up an AD query dict by bcftools mpileup
    bam_ad_file = f"{bam_file}.ad"
    cmd = f"""bcftools mpileup -Ou --fasta-ref {reference_genome} -a FORMAT/AD --indels-2.0 -q 10 -Q 10 {bam_file} | \
              bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n' - > {bam_ad_file}"""
    executeCmd(cmd, logger = logger)
    ad_table = pd.read_table(bam_ad_file, header = None, sep = "\t", names = ["chrom", "pos", "ref", "alt", "ad"], dtype = {"chrom": str, "pos": int, "ref": str, "alt": str, "ad": str}, na_values=["", "<*>"])
    if ad_table.dropna(subset = ["alt"]).shape[0] < 1:
        logger.warning(f"No ALT allele found in this BAM file. Skip this entire script")
        logger.warning(f"Check the {bam_ad_file} or take a look at the original dataframes now: \n{ad_table.to_string(index=False)}\n")
        return None
    ad_table.dropna(subset = ["alt"], inplace = True)
    ad_expanded = ad_table["ad"].str.split(",", expand=True).replace({None: np.nan, "": np.nan, "0": np.nan}).astype(float)
    # Fix: Instead of using subset, we'll filter columns directly
    if ad_expanded.shape[1] > 1:
        # Keep column 0 (reference) and any non-alt columns that have at least one non-NA value
        cols_to_keep = [0]  # Always keep reference column
        for col in ad_expanded.columns[1:]:
            if not ad_expanded[col].isna().all():
                cols_to_keep.append(col)
        ad_expanded = ad_expanded[cols_to_keep]
    
    alt_expanded = ad_table["alt"].str.rstrip(",<*>").str.split(",", expand=True).replace({None: np.nan, "": np.nan}).dropna(axis = 1, how = "all")

    logger.debug("The AD expanded table looks like: \n{}\nThe ALT expanded table looks like: \n{}\n".format(ad_expanded.loc[~ad_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False),
                                                                                                            alt_expanded.loc[~alt_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False)))

    # Initialize the nested dictionary
    nested_ad_dict = {chrom: {} for chrom in ad_table["chrom"].unique()}
    column_width = alt_expanded.shape[1]

    # Not intended to handle insertions and deletions because they are rarely are due to low quality single bases. 
    base_dict["DP"] = np.int8(5)

    # Iterate over the rows of dfA and dfB
    for i in range(len(ad_table)):
        chrom = ad_table.iloc[i, 0]
        position = ad_table.iloc[i, 1] - 1 # position, output of bcftools mpileup is 1-based, now needs to be 0-based
        ref_allele = ad_table.iloc[i, 2]
        alt_alleles = alt_expanded.iloc[i, :] # first allele
        ref_depth = np.int16(ad_expanded.iloc[i, 0])
        alt_depths = np.int16(ad_expanded.iloc[i, 1:])

        if pd.isna(ref_depth):
            logger.warning(f"The reference depth is NA for the row {i} at {chrom}:{position} of the ad table which is: \n{ad_table.iloc[i, :].to_string(index=False)}\n")
            continue

        if pd.isna(alt_depths).all():
            logger.warning(f"The alt depths are all NA for the row {i} at {chrom}:{position} of the ad table which is: \n{ad_table.iloc[i, :].to_string(index=False)}\n")
            continue

        total_dp = ad_expanded.iloc[i, :].dropna().sum()

        if total_dp == 0 or pd.isna(total_dp):
            logger.warning(f"The total depth is 0 for the row {i} at {chrom}:{position} of the ad table which is: \n{ad_table.iloc[i, :].to_string(index=False)}\n")
            continue

        # Initialize the inner dictionary if the outer key is not present
        if position not in nested_ad_dict[chrom]:
            nested_ad_dict[chrom][position] = empty_dict.copy()

        nested_ad_dict[chrom][position][base_dict["DP"]] = np.int16(total_dp)
    
        # Check for the second pair of inner key-value
        for c in range(0, column_width):
            alt_allele = alt_alleles[c]
            alt_depth = alt_depths[c]
            if pd.isna(alt_allele):
                continue
            if c >= alt_depths.size:
                continue
            if pd.isna(alt_depth):
                continue
            if len(ref_allele) > len(alt_allele):
                continue
            if len(ref_allele) < len(alt_allele):
                continue
            if len(ref_allele) == len(alt_allele) and len(alt_allele) > 1:
                diff_indx = [i for i in range(len(ref_allele)) if ref_allele[i] != alt_allele[i]]
                assert len(diff_indx) == 1, f"The reference allele and the alternate allele are of the same length but there are more than one difference at {chrom}:{position} of the ad table which is: \n{ad_table.iloc[i, :].to_string(index=False)}\n"
                logger.warning(f"The reference allele and the alternate allele are of the same length but there are more than 1 bases at the allele and they should differ at {chrom}:{position} of the ad table which is: \n{ad_table.iloc[i, :].to_string(index=False)}\n")
                diff_indx = diff_indx[0]
                alt_allele = alt_allele[diff_indx]

            # Now the remaining ref-alt pair is a SNV
            nested_ad_dict[chrom][position][base_dict[alt_allele]] = np.int16(alt_depth)
            # logger.debug(f"The allele depth for base {alt_allele} at {chrom}:{position} is {alt_depth} and the allele depth looks like {nested_ad_dict[chrom][position]}")

    return nested_ad_dict


def check_edge(u, v, adj_set):
    '''
    1 means not sure
    2 means accept same haplotype
    -1 means reject same haplotype

    Cant use 0 here because python treats 0 as False
    '''
    return adj_set.get((u, v), False) or adj_set.get((v, u), False)


def get_overlap_intervals(read_pair1, read_pair2, logger = logger):
    '''
    This function is to extract the overlapping intervals between two pairs of reads
    '''
    overlap_intervals = {}
    # logger.debug(f"Trying to find overlap intervals between the read pair {[get_read_id(x) for x in read_pair1]} and {[get_read_id(x) for x in read_pair2]}")

    for r1 in read_pair1:
        for r2 in read_pair2:
            r1_start = r1.reference_start
            r2_end = r2.reference_end

            if r2_end <= r1_start:
                # logger.debug(f"Read {get_read_id(r1)} and {get_read_id(r2)} have no overlap")
                continue

            r1_end = r1.reference_end
            r2_start = r2.reference_start

            if r1_end <= r2_start:
                # logger.debug(f"Read {get_read_id(r1)} and {get_read_id(r2)} have no overlap")
                continue

            overlap_start = max(r1_start, r2_start)
            overlap_end = min(r1_end, r2_end)

            # logger.debug(f"Found the overlap interval {overlap_start}-{overlap_end} between the read pair {get_read_id(r1)} and {get_read_id(r2)}")
            overlap_intervals[(overlap_start, overlap_end)] = (r1, r2)

    return overlap_intervals



class FastIntervals:
    def __init__(self, max_size = 10):
        self.starts = np.empty(max_size, dtype=np.int32)
        self.ends = np.empty(max_size, dtype=np.int32)
        self.size = 0
    
    def add(self, start, end):
        """Add new interval maintaining sorted order"""
        self.starts[self.size] = start
        self.ends[self.size] = end
        self.size += 1

    def clear(self):
        """Reset arrays"""
        self.starts = np.array([], dtype=np.int32)
        self.ends = np.array([], dtype=np.int32)



@numba.njit(types.int32[:, :](types.int32[:], types.int32[:], types.int32, types.int32), fastmath=True)
def find_uncovered_regions_numba(existing_starts, existing_ends, new_start, new_end):
    n = len(existing_starts)
    if n == 0:
        return np.array([[new_start, new_end]], dtype=np.int32)
        
    # Find overlapping intervals using binary search
    overlaps = (existing_starts <= new_end) & (existing_ends >= new_start)
    
    if numba_sum(overlaps) == 0:
        return np.array([[new_start, new_end]], dtype=np.int32)
    
    # Pre-allocate result array
    result = np.empty((len(existing_starts) + 1, 2), dtype=np.int32)
    n_results = 0
    current_start = new_start
    
    for i in range(len(existing_starts)):
        if overlaps[i]:
            if existing_starts[i] > current_start:
                result[n_results, 0] = current_start
                result[n_results, 1] = existing_starts[i]
                n_results += 1
            current_start = max(current_start, existing_ends[i])
    
    if current_start < new_end:
        result[n_results, 0] = current_start
        result[n_results, 1] = new_end
        n_results += 1
    
    return result[:n_results]




def build_phasing_graph(bam_file,
                        intrinsic_bam,
                        ncls_dict,
                        ncls_read_dict,
                        ncls_qname_dict,
                        mean_read_length,
                        total_lowqual_qnames,
                        reference_genome,
                        base_dict = {"A": np.int8(0), "T": np.int8(1), "C": np.int8(2), "G": np.int8(3), "N": np.int8(4)},
                        logger = logger):
    '''
    Construct a phasing graph from BAM data for efficient haplotype identification.

    This function builds a graph where each vertex represents a read pair and each edge
    represents the confidence that two read pairs originated from the same haplotype.
    The graph is used for subsequent haplotype identification through clique finding.

    Parameters:
    -----------
    bam_file : str
        Path to the input BAM file.

    ncls_dict : dict
        Dictionary of NCLS objects, keyed by chromsome names.
        The ncls_dict looks like:
                                        { chromsome1: ncls_chr1,
                                          chromsome2: ncls_chr2,
                                          ...
                                          chromsomeN: ncls_chrN  }

    ncls_read_dict : dict
        Dictionary of read objects, keyed by query name indices (which are also the interval indices stored in the NCLS object)
        The ncls_read_dict looks like:
                                        { qname_idx1: [read1, read2],
                                          qname_idx2: [read1, read2]  }

    ncls_qname_dict : dict
        Dictionary mapping query name indices to query names. (qname_indices -> qnames)

    ncls_qname_idx_dict : dict
        Dictionary mapping query names to query name indices. (qnames -> qname_indices)

    mean_read_length : float
        Average read length in the BAM file.
    edge_weight_cutoff : float, optional
        Minimum edge weight threshold for including an edge in the graph (default is 0.201, a heuristic cutoff).
    logger : logging.Logger, optional
        Logger object for output messages.

    Returns:
    --------
    tuple
        A tuple containing:
        - phased_graph : graph_tool.Graph
            The constructed phasing graph.
        - weight_matrix : numpy.ndarray
            Matrix of edge weights between read pairs. (adjacency matrix)
        - qname_to_node : dict
            Dictionary mapping query names to graph vertices. (qnames -> graph_vertices_indices)
        - total_readhap_vector : dict
            Dictionary of read haplotype vectors. (read_ids -> haplotype_vectors), please refer read_ids to function get_read_id

    Notes:
    ------
    - The function uses the graph-tool library for efficient graph operations.
    - Edge weights are calculated based on two factors:
        1. The total overlapping span between two read pairs
        2. The number of variants (SNVs and small indels) shared by two read pairs
    - The resulting graph is used for subsequent haplotype identification through clique finding (Greedy-Clique-Expansion algorithm).

    See Also:
    ---------
    migrate_bam_to_ncls : For creating the NCLS data structures.
    find_cliques_in_components : For identifying haplotypes in the constructed graph.
    '''

    logger.info(f"There are totally {len(ncls_read_dict)} pair of reads, mean read length is {mean_read_length}. with adequate mapping or base quality which can be used to build the graph")
    # Create an empty graph
    g = gt.Graph(directed = False)
    g.set_fast_edge_removal(fast = True)

    # Create a property map to store the query names for each node
    qname_prop = g.new_vertex_property("string")
    weight = g.new_edge_property("float")

    # Create a dictionary to map query names to their corresponding nodes
    qname_to_node = {}

    # Create a set to store the hap_vectors corresponding to each read
    read_hap_vectors = {}
    read_error_vectors = {}
    # Create a dictionary to store the start and end positions for each read pair
    read_ref_pos_dict = {}

    # Create a dictionary to store Allele Depth for each position
    empty_dict = Dict.empty(key_type=types.int8, value_type=types.int16)
    # Make sure every value is 0 to initiate the empty dict
    for key in empty_dict.keys():
        empty_dict[key] = np.int16(0)
    nested_ad_dict = stat_ad_to_dict(bam_file, empty_dict, reference_genome, base_dict, logger = logger)
    if nested_ad_dict is None:
        logger.warning(f"No ALT allele found in this BAM file. Skip this entire script")
        return None, None, None, None, None, None, total_lowqual_qnames
    
    intrinsic_ad_dict = stat_ad_to_dict(intrinsic_bam, empty_dict, reference_genome, base_dict, logger = logger)
    if intrinsic_ad_dict is None: intrinsic_ad_dict = {}

    qname_check_dict = {}
    # Use NCLS read dict to iterate through all the reads to build a graph. One good thing is that every key: value stores a pair of read objects
    total_qname_num = len(ncls_read_dict)
    weight_matrix = np.eye(total_qname_num, dtype=np.float32)
    score_arr = np.array([mean_read_length + mean_read_length * i for i in range(50)])

    for qname_idx, paired_reads in ncls_read_dict.items():
        qname = ncls_qname_dict[qname_idx]
        assert len(dict.fromkeys(r.query_name for r in paired_reads)) == 1, f"The query names of the reads in the paired_reads are not the same: {paired_reads}"
        # Check if the query name already exists in the graph
        if qname not in qname_to_node:
            # Add a new node to the graph
            qv = g.add_vertex()
            qname_prop[qv] = qname
            qname_to_node[qname] = int(qv)
            # logger.debug(f"Added a new node {v} for qname {qname} to the graph. The current vertices are {g.get_vertices()}")
        else:
            qv = g.vertex(qname_to_node[qname])

        chrom = paired_reads[0].reference_name
        start = min(r.reference_start for r in paired_reads)
        end = max(r.reference_end for r in paired_reads)
        qidx_iter = overlap_qname_idx_iterator(ncls_dict, 
                                               chrom, start, end)

        # Iterate through the reads
        for qidx in qidx_iter:
            other_reads = ncls_read_dict[qidx]
            # Get the query name
            other_qname = other_reads[0].query_name

            if qname == other_qname:
                continue
                
            # Check if the query name already exists in the graph
            if not other_qname in qname_to_node:
                # Add a new node to the graph
                oqv = g.add_vertex()
                qname_prop[oqv] = other_qname
                qname_to_node[other_qname] = int(oqv)
                # logger.debug(f"Added a new node {v} for qname {qname} to the graph. The current vertices are {g.get_vertices()}")
            else:
                oqv = g.vertex(qname_to_node[other_qname])

            inspect_res = check_edge(int(qv), int(oqv), qname_check_dict)
            if inspect_res:
                continue

            qname_check_dict[(int(qv), int(oqv))] = True

            if qname in total_lowqual_qnames:
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1
                # logger.debug(f"The read pair {qname} is low quality, skip inspecting this pair of fragments")
                continue

            if other_qname in total_lowqual_qnames:
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1
                # logger.debug(f"The read pair {other_qname} is low quality, skip inspecting this pair of fragments")
                continue

            # Inspect the overlap between the two pairs of reads
            overlap_intervals = get_overlap_intervals(paired_reads, other_reads, logger = logger)
            # logger.debug(f"Found these overlap intervals between the read pair {qname} and {other_qname}: {overlap_intervals}")

            if len(overlap_intervals) == 0:
                continue

            # logger.info("These are the overlapping intervals: {}".format("\n".join([f"{oi}: {[get_read_id(r) for r in (read1, read2)]}" for oi, (read1, read2) in overlap_intervals.items()])))

            qname_bools = np.zeros(4, dtype=np.int32)
            n = 0
            pair_weight = None
            # Clear all the contents inside the inspected_overlaps
            inspected_overlaps = FastIntervals(max_size = len(overlap_intervals))

            for (overlap_start, overlap_end), (read1, read2) in overlap_intervals.items():
                uncovered_overlaps = find_uncovered_regions_numba(inspected_overlaps.starts[:inspected_overlaps.size], inspected_overlaps.ends[:inspected_overlaps.size], 
                                                                  overlap_start, overlap_end)
                # logger.debug(f"Found the uncovered regions {uncovered_overlaps} for the reads {read1.query_name} and {read2.query_name} in the overlap region ({chrom}:{overlap_start}-{overlap_end})")
                for row_ind in range(uncovered_overlaps.shape[0]):
                    if other_qname in total_lowqual_qnames:
                        # logger.debug(f"The read pair {other_qname} is low quality, skip this pair of reads")
                        break
                    if qname in total_lowqual_qnames:
                        # logger.debug(f"The read pair {qname} is low quality, skip this pair of reads")
                        break
                    uncovered_start, uncovered_end = uncovered_overlaps[row_ind, 0], uncovered_overlaps[row_ind, 1]
                    bool_res, read_ref_pos_dict, read_hap_vectors, read_error_vectors, total_lowqual_qnames, read_weight = determine_same_haplotype(read1, read2,
                                                                                                                                                    uncovered_start, uncovered_end,
                                                                                                                                                    score_arr,
                                                                                                                                                    read_hap_vectors = read_hap_vectors,
                                                                                                                                                    read_error_vectors = read_error_vectors,
                                                                                                                                                    nested_ad_dict = nested_ad_dict[chrom],
                                                                                                                                                    read_ref_pos_dict = read_ref_pos_dict,
                                                                                                                                                    total_lowqual_qnames = total_lowqual_qnames,
                                                                                                                                                    mean_read_length = mean_read_length,
                                                                                                                                                    empty_dict = empty_dict,
                                                                                                                                                    intrinsic_ad_dict = intrinsic_ad_dict.get(chrom, {}),
                                                                                                                                                    base_dict = base_dict,
                                                                                                                                                    logger = logger)
                    
                    if np.isnan(bool_res):
                        qname_bools[n] = 0
                    elif bool_res:
                        qname_bools[n] = 1
                        # logger.info(f"Weight is {read_weight}. Found the reads {read1.query_name} ({read1.reference_name}:{read1.reference_start}-{read1.reference_end}) and {read2.query_name} ({read2.reference_name}:{read2.reference_start}-{read2.reference_end}) are in the same haplotype, overlap region ({overlap_start}-{overlap_end})")
                        # assert read_weight is not None, f"The read weight {read_weight} is not calculated for the reads {read1.query_name} and {read2.query_name}"
                        # logger.info(f"Found the reads {read1.query_name} ({read1.reference_name}:{read1.reference_start}-{read1.reference_end}) and {read2.query_name} ({read2.reference_name}:{read2.reference_start}-{read2.reference_end}) are in the same haplotype, overlap region ({overlap_start}-{overlap_end})")
                    else:
                        qname_bools[n] = -1
                        # logger.info(f"Found the reads {read1.query_name} ({read1.reference_name}:{read1.reference_start}-{read1.reference_end}) and {read2.query_name} ({read2.reference_name}:{read2.reference_start}-{read2.reference_end}) are in different haplotypes, overlap region ({overlap_start}-{overlap_end})")

                    if read_weight is not None:
                        read_weight = read_weight if read_weight > 0 else 0
                        norm_weight = read_weight/(mean_read_length * 10)
                        if pair_weight is None:
                            pair_weight = norm_weight
                        else:
                            pair_weight += norm_weight
                    else:
                        norm_weight = None
                    
                inspected_overlaps.add(overlap_start, overlap_end)
                n += 1

            if other_qname in total_lowqual_qnames or qname in total_lowqual_qnames:
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1
                continue

            qname_bools = qname_bools[:n]
            pair_weight = 0 if pair_weight is None else pair_weight
            pair_weight = pair_weight if pair_weight > 0 else 1e-4
            
            if any_false_numba(qname_bools):
                # logger.info(f"Qname_bools are {qname_bools}, Found two pairs {qname_prop[qv]} and {qname_prop[oqv]} are in different haplotypes, Removing the edge with the biggest weight")
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1
            else:
                # logger.debug(f"Between {qname_prop[qv]} and {qname_prop[oqv]}, the pair weight is {pair_weight}")
                weight_matrix[int(qv), int(oqv)] = pair_weight
                weight_matrix[int(oqv), int(qv)] = pair_weight
                e = g.add_edge(qv, oqv)
                weight[e] = pair_weight

    # Set the query name property for the graph
    g.vertex_properties["qname"] = qname_prop
    g.edge_properties["weight"] = weight

    logger.info(f"Now we finished building up the edges in the graph. There are currently {g.num_vertices()} vertices and {g.num_edges()} edges in the graph")
    return g, weight_matrix, qname_to_node, read_hap_vectors, read_error_vectors, read_ref_pos_dict, total_lowqual_qnames