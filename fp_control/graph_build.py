import logging
import graph_tool.all as gt
import numpy as np
import pandas as pd

from intervaltree import IntervalTree

from shell_cmds import executeCmd
from numba_operators import any_false_numba, custom_all_numba
from bam_ncls import overlapping_qname_idx_generator
from pairwise_read_inspection import determine_same_haplotype


logger = logging.getLogger('SDrecall')



def stat_ad_dict(bam_file, logger = logger):
    '''
    This function is to extract the allele depth (AD) information from the BAM file
    The AD info is stored in a nested dictionary, which is then used to calculate the edge weights in the graph
    '''
    bam_ad_file = f"{bam_file}.ad"
    cmd = f"""bcftools mpileup -Ou --no-reference -a FORMAT/AD --indels-2.0 -q 10 -Q 15 {bam_file} | \
              bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n' - > {bam_ad_file}"""
    executeCmd(cmd, logger = logger)
    ad_table = pd.read_table(bam_ad_file, header = None, sep = "\t", names = ["chrom", "pos", "alt", "ad"], dtype = {"chrom": str, "pos": int, "alt": str, "ad": str}, na_values=["", "<*>"])
    ad_expanded = ad_table["ad"].str.split(",", expand=True).replace({None: np.nan, "": np.nan, "0": np.nan}).dropna(axis = 1, how="all").astype(float)
    alt_expanded = ad_table["alt"].str.rstrip(",<*>").str.split(",", expand=True).replace({None: np.nan, "": np.nan}).dropna(axis = 1, how="all")

    if alt_expanded.shape[1] <= 1 or ad_expanded.shape[1] <= 1:
        logger.warning(f"No ALT allele found in this BAM file. Skip this entire script")
        logger.warning(f"Look at the two dataframes now: \n{alt_expanded.to_string(index=False)}\n{ad_expanded.to_string(index=False)}\n")
        return None

    logger.info("\n{}\n{}\n".format(ad_expanded.loc[~ad_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False),
                                    alt_expanded.loc[~alt_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False)))

    # After getting the table, we seek to organize the AD info by chromsome, pos, and alt into a nested dictionary
    # The dict structure is like:
    # nested_ad_dict[chrom][pos][alt] = ad
    # nested_ad_dict[chrom][pos]["DP"] = total_dp
    # Initialize the nested dictionary
    nested_ad_dict = {chrom: {} for chrom in ad_table["chrom"].unique()}
    column_width = alt_expanded.shape[1]

    # Iterate over the rows of ad_table to build the nested_ad_dict
    for i in range(len(ad_table)):
        chrom = ad_table.iloc[i, 0]
        outer_key = ad_table.iloc[i, 1]
        inner_key = alt_expanded.iloc[i, 0]
        value = ad_expanded.iloc[i, 0]

        if pd.isna(value):
            continue

        total_dp = value

        # Initialize the inner dictionary if the outer key is not present
        if outer_key not in nested_ad_dict[chrom]:
            nested_ad_dict[chrom][outer_key] = {}

        # Add the first pair of inner key-value
        nested_ad_dict[chrom][outer_key][inner_key] = value

        # Check for the second pair of inner key-value
        for c in range(1, column_width):
            if not pd.isna(ad_expanded.iloc[i, c]) and not pd.isna(alt_expanded.iloc[i, c]):
                nested_ad_dict[chrom][outer_key][alt_expanded.iloc[i, c]] = ad_expanded.iloc[i, c]
                total_dp += ad_expanded.iloc[i, c]

        nested_ad_dict[chrom][outer_key]["DP"] = total_dp

    return nested_ad_dict




def check_edge(u, v, adj_set):
    '''
    1 means not sure
    2 means accept same haplotype
    -1 means reject same haplotype

    Cant use 0 here because python treats 0 as False
    '''
    return adj_set.get((u, v), False) or adj_set.get((v, u), False)





def get_overlap_intervals(read_pair1, read_pair2):
    '''
    This function is to extract the overlapping intervals between two pairs of reads
    '''
    overlap_intervals = {}

    for r1 in read_pair1:
        for r2 in read_pair2:
            r1_start = r1.reference_start
            r2_end = r2.reference_end

            if r2_end - r1_start <= 0 or r2_end - r1_start >= 300:
                continue

            r1_end = r1.reference_end
            r2_start = r2.reference_start

            if r1_end - r2_start <= 0:
                continue

            overlap_start = max(r1_start, r2_start)
            overlap_end = min(r1_end, r2_end)

            overlap_intervals[(overlap_start, overlap_end)] = (r1, r2)

    return overlap_intervals




def find_uncovered_regions(existing_intervals, new_interval):
    # Find overlapping intervals
    overlapping_intervals = existing_intervals.overlap(new_interval[0], new_interval[1])

    # Sort the overlapping intervals by start time
    sorted_intervals = sorted(overlapping_intervals, key=lambda x: x.begin)

    # Determine the uncovered regions
    uncovered_regions = []
    current_start = new_interval[0]

    if len(overlapping_intervals) == 0:
        return [new_interval]

    for interval in sorted_intervals:
        if interval.begin > current_start:
            uncovered_regions.append((current_start, interval.begin))
        current_start = max(current_start, interval.end)

    if current_start < new_interval[1]:
        uncovered_regions.append((current_start, new_interval[1]))

    return uncovered_regions




def build_phasing_graph(bam_file,
                        ncls_dict,
                        ncls_read_dict,
                        ncls_qname_dict,
                        ncls_qname_idx_dict,
                        mean_read_length,
                        edge_weight_cutoff = 0.201,
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
    - The resulting graph is used for subsequent haplotype identification through clique finding (Bron-Kerbosch algorithm).

    See Also:
    ---------
    migrate_bam_to_ncls : For creating the NCLS data structures.
    clique_generator_per_component : For identifying haplotypes in the constructed graph.
    '''

    logger.info(f"There are totally {len(ncls_read_dict)} pair of reads, mean read length is {mean_read_length}. with adequate mapping or base quality which can be used to build the graph")
    # Create an empty graph
    g = gt.Graph(directed = False)
    g.set_fast_edge_removal(fast = True)

    # Create a property map to store the query names for each node
    qname_prop = g.new_vertex_property("string")
    weight = g.new_edge_property("float")

    # When we try to identify a mismatch between two reads are due to a true variant or due to sequencing errors,
    # we need to know the allele depth at the mismatch site. This allows us to estimate the sequencing artifact likelihood by observing the allele's depth fraction at this site.
    # We would like an efficient way to extract AD info for all covered sites across this BAM. So we chose bcftools mpileup + bcftools query.
    # The output file is stored in the same directory as the input BAM file
    nested_ad_dict = stat_ad_dict(bam_file, logger = logger)
    if nested_ad_dict is None:
        return None, None, None, None

    '''
    Now we start to build an undirected graph of read pairs,
    The graph is constructed as follows:
    - Each read pair is represented by a vertex
    - Iterate through all the read pairs in the ncls_read_dict
       - For each read pair, extract its overlapping read pairs
       - Iterate through all the overlapping read pairs, skip the read pairs that have been inspected (undirected)

    '''

    # Create a dictionary to store tuples of vertex indices (v1, v2) as keys, and the value is a boolean indicating whether the pair of read pairs has been inspected
    # This is a temporary data hub to avoid repeated edge inspection for the same pair of reads
    qname_check_dict = {}

    # Initialize the weight matrix (adjacency matrix, which is diagonal symmetric) with np.eye
    total_qname_num = len(ncls_read_dict)
    weight_matrix = np.eye(total_qname_num, dtype=np.float32)

    # Create a dictionary to store the haplotype vectors corresponding to each read (read_ids -> haplotype_vectors), please refer read_ids to function get_read_id
    # This is a temporary data hub to avoid repeated hapvector extraction from the same read
    # The variable is primarily aggregated in the iterative calling of function determine_same_haplotype
    read_hap_vectors = {}

    # Create a dictionary to store the start and end positions for each read pair
    # This is a temporary data hub to avoid repeated read pair start and end position extraction from the same read
    # The variable is primarily aggregated in the iterative calling of function determine_same_haplotype
    read_ref_pos_dict = {}

    # Create a dictionary to map query names to their corresponding nodes (qnames -> graph_vertices_indices)
    # This is a temporary data hub to record whether a read pair has been added to the graph as a vertex(node)
    qname_to_node = {}

    # Create an IntervalTree to record the inspected overlapping intervals
    inspected_overlaps = IntervalTree()

    for qname_idx, paired_reads in ncls_read_dict.items():
        qname = ncls_qname_dict[qname_idx]
        assert len(dict.fromkeys(r.query_name for r in paired_reads)) == 1, f"The query names of the reads in the paired_reads are not the same: {paired_reads}"

        # Check if the query name already exists in the graph
        if qname not in qname_to_node:
            # Add a new node to the graph in case that the read pair has not been added to the graph
            qv = g.add_vertex()
            qname_prop[qv] = qname
            qname_to_node[qname] = int(qv)
        else:
            # Directly retrieve the existing vertex (graph-tool vertex object) from the graph
            qv = g.vertex(qname_to_node[qname])

        # Extract the overlapping span of the current iterating read pair
        # And get a generator of qname indices of the read pairs overlapping with the current read pair
        # Note that the span also include the inner gap between a pair of reads
        chrom = paired_reads[0].reference_name
        start = min(r.reference_start for r in paired_reads)
        end = max(r.reference_end for r in paired_reads)
        qidx_iter = overlapping_qname_idx_generator(ncls_dict, chrom, start, end)

        # Iterate through the overlapping read pairs
        for qidx in qidx_iter:
            other_reads = ncls_read_dict[qidx]
            # Get the query name
            other_qname = other_reads[0].query_name
            # When extracting overlapping reads from NCLS, the result qnames might contain the iterating qname itself
            if qname == other_qname:
                continue

            '''
            # Below is a code block for debugging, this allows you to check whether two read pairs overlap with each other during the inspection of them
            assert len(dict.fromkeys(r.query_name for r in other_reads)) == 1, f"The query names of the reads in the other_reads are not the same: {other_reads}"
            vis_qnames = ["HISEQ1:21:H9V1VADXX:2:1112:21127:38421:PC0",
                            "HISEQ1:26:HA2RRADXX:1:1113:11601:32503:PC0"]
            if qname in vis_qnames and other_qname in vis_qnames:
                logger.info(f"Found the qname {qname} and other_qname {other_qname} might overlap with each other.")
            '''

            # Check if the query name already exists in the graph
            if not other_qname in qname_to_node:
                # Add a new node to the graph
                oqv = g.add_vertex()
                qname_prop[oqv] = other_qname
                qname_to_node[other_qname] = int(oqv)
                # logger.debug(f"Added a new node {v} for qname {qname} to the graph. The current vertices are {g.get_vertices()}")
            else:
                oqv = g.vertex(qname_to_node[other_qname])

            # Check if the pair of read pairs has been inspected, if yes, skip the current iteration
            inspect_res = check_edge(int(qv), int(oqv), qname_check_dict)
            if inspect_res:
                continue

            # Proceeding to this step, we mark the pair of read pairs as inspected
            qname_check_dict[(int(qv), int(oqv))] = True

            # Inspect the overlap between the two pair of reads
            # There might be multiple overlapping intervals between two pair of reads
            overlap_intervals = get_overlap_intervals(paired_reads, other_reads)
            # logger.info(f"Found these overlap intervals for the read pair {qname} and {other_qname}: {overlap_intervals}")

            # If there is no overlap between two pair of reads, skip the current iteration
            # This might be some edge case when a read in pair A completely enclosed by the inner gap of pair B, while the other read in pair A exceeds the coverage span of pair B
            if len(overlap_intervals) == 0:
                continue

            '''
            There are four comparison results (C 2 2) to be recorded between two pair of reads
            - read1 in pair A and read1 in pair B comparison result
            - read1 in pair A and read2 in pair B comparison result
            - read2 in pair A and read1 in pair B comparison result
            - read2 in pair A and read2 in pair B comparison result

            There are also three possible results for the comparison of two reads
             1: meaning we found two reads from distinct pairs overlap with each other with enough size, (or share variants), so they should be in the same haplotype
             0: meaning we found two reads from distinct pairs do not overlap with each other or the overlapping span is too small (e.g < 100bp). So we do not have enough evidence to determine if they are in the same haplotype or not
            -1: meaning we found two reads from distinct pairs overlap with each other and we found clear evidence that they are not in the same haplotype (e.g. different variants)
            '''

            # Initialize the same_hap_bools with 4 zeros
            same_hap_bools = np.zeros(4, dtype=np.int32)

            n = 0
            pair_weight = None

            # Remember the overlap intervals we extracted for two pair of reads
            # Now we iterate through all the overlap intervals and use an IntervalTree to record the inspected intervals
            inspected_overlaps = IntervalTree()
            for (overlap_start, overlap_end), (read1, read2) in overlap_intervals.items():
                '''
                First we identify whehter the iterating overlap interval has some part of it that has been inspected
                If yes, we crop out the inspected part and only keep the uncovered overlapping span for downstream analysis
                Note that after cropping out the inspected part, the remaining span might not be continuous, leaving multiple uncovered intervals
                For example, the inspected intervals is [100, 200], [300, 400], [500, 600]
                If the current iterating overlap interval is [250, 550], then the remaining uncovered intervals are [250, 400] and [500, 550]
                The function IntervalTree.addi will return the following uncovered intervals: [250, 400], [500, 550]
                '''
                uncovered_overlaps = find_uncovered_regions(inspected_overlaps, (overlap_start, overlap_end))
                # logger.info(f"Found the uncovered regions {uncovered_overlaps} for the reads {read1.query_name} and {read2.query_name} in the overlap region ({overlap_start}-{overlap_end})")
                for uncovered_start, uncovered_end in uncovered_overlaps:
                    # Iterate over one uncovered interval at a time, and determine whether two read from distinct pairs overlapping at the current interval show evidence of being in the same haplotype or not
                    bool_res, read_ref_pos_dict, read_hap_vectors, read_weight = determine_same_haplotype(read1, read2,
                                                                                                          uncovered_start, uncovered_end,
                                                                                                          read_hap_vectors = read_hap_vectors,
                                                                                                          nested_ad_dict = nested_ad_dict,
                                                                                                          read_ref_pos_dict = read_ref_pos_dict,
                                                                                                          mean_read_length = mean_read_length,
                                                                                                          logger = logger)
                    # If the bool_res is NaN, we assign 0 to record the inspection result within this overlapping interval
                    # Otherwise, we assign 1 or -1 based on the bool_res
                    if np.isnan(bool_res):
                        same_hap_bools[n] = 0
                    elif bool_res:
                        same_hap_bools[n] = 1
                    else:
                        same_hap_bools[n] = -1

                    if read_weight is not None:
                        # Assign 1 (an aritificial very small number) to read_weight unless it a positive value
                        read_weight = read_weight if read_weight > 0 else 1
                        # Normalize the read_weight to adjust the scale of the weight, why denominate by mean_read_length * 10?
                        # mean_read_length * 10 is just a big enough number to make sure ur normal weight is at the scale between 0 and 1, u can understand it as a maximum normalization factor
                        norm_weight = read_weight/(mean_read_length * 10)
                        if pair_weight is None:
                            # If pair_weight is None, assign norm_weight to it
                            pair_weight = norm_weight
                        else:
                            # If pair_weight is not None, add norm_weight to it
                            pair_weight += norm_weight
                    else:
                        norm_weight = None

                # After inspecting the current uncovered interval, we add it to the inspected_overlaps
                # And adding 1 to n for next overlapping interval.
                inspected_overlaps.addi(overlap_start, overlap_end)
                n += 1

            # Now we iterated through all the overlapping intervals, we can trim the same_hap_bools to the actual number of overlapping intervals
            same_hap_bools = same_hap_bools[:n]

            # Because 1 in the adjacency matrix has a special meaning. The diagonal values are all 1 (representing self-self comparison)
            # So we add a very small number to the pair_weight to distinguish it from 1
            pair_weight = pair_weight if pair_weight != 1 else 1 + 1e-4

            # Below is just a debug code block
            # if qname in vis_qnames and other_qname in vis_qnames:
            #     logger.info(f"Found the qname {qname} (qv is {int(qv)}) and other_qname {other_qname} (oqv is {int(oqv)}) overlap with each other with pair_weight {pair_weight}. The same_hap_bools are {same_hap_bools}")

            # We need to consider all the overlapping intervals between two read pairs to decide the value we record inside the adjacency matrix
            if custom_all_numba(same_hap_bools):
                # logger.info(f"The same_hap_bools are {same_hap_bools}")
                if pair_weight > edge_weight_cutoff:
                    e = g.add_edge(qv, oqv)
                    weight[e] = pair_weight
                weight_matrix[int(qv), int(oqv)] = pair_weight
                weight_matrix[int(oqv), int(qv)] = pair_weight
            elif any_false_numba(same_hap_bools):
                # We found clear evidence that two read pairs are not in the same haplotype
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1

    # Set the query name property for the graph
    g.vertex_properties["qname"] = qname_prop
    g.edge_properties["weight"] = weight

    logger.info(f"Now we finished building up the edges in the graph. There are currently {g.num_vertices()} vertices and {g.num_edges()} edges in the graph")
    return g, weight_matrix, qname_to_node, read_hap_vectors