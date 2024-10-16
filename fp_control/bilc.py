import numba
import numpy as np
import pandas as pd
from numba import prange, types
from highspy import Highs
from collections import defaultdict



def overlapping_reads_generator(ncls_dict, read_dict, qname_dict, chrom, start, end):
    """
    Generator function to lazily yield overlapping read objects.

    Parameters:
    -----------
    ncls_dict : dict
        Dictionary of NCLS objects, keyed by chromosome names.
    read_dict : dict
        Dictionary of read objects, keyed by query name indices.
    qname_dict : dict
        Dictionary mapping query name indices to query names.
    chrom : str
        Chromosome name.
    start : int
        Start position of the interval (0-based, inclusive).
    end : int
        End position of the interval (0-based, exclusive).

    Yields:
    -------
    pysam.AlignedSegment
        Overlapping read object.

    Notes:
    ------
    This function uses NCLS for efficient overlap queries and yields individual reads.
    It filters reads to ensure they actually overlap with the given interval.
    """

    # Perform an overlap query on the NCLS tree
    if chrom not in ncls_dict:
        yield from ()
    elif ncls_dict[chrom] is None:
        yield from ()
    else:
        _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
        # overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
        for qname_idx in list(dict.fromkeys(overlapping_read_qnames)):
            for read in read_dict[qname_idx]:
                if read.reference_start < end and read.reference_end > start:
                    yield read



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




def identify_misalignment_per_region(region,
                                    bam_ncls,
                                    intrin_bam_ncls,
                                    phasing_graph,
                                    qname_hap_info,
                                    hap_qname_info,
                                    weight_matrix,
                                    qname_to_node,
                                    lowqual_qnames,
                                    clique_sep_component_idx,
                                    var_df,
                                    viewed_regions = {},
                                    total_hapvectors = {},
                                    total_errvectors = {},
                                    total_genomic_haps = {},
                                    mean_read_length = 148,
                                    logger = logger ):
    bam_ncls, read_pair_dict, qname_dict, qname_idx_dict = bam_ncls
    '''
    read pair dict use internal qname_idx (created when executing function migrate bam to ncls) as keys instead of qnames
    qname_dict is a dictionary that maps qname_idx to qname
    '''

    # qname_idx_dict = {qname: qname_idx for qname_idx, qname in qname_dict.items()}
    overlap_reads_iter = overlapping_reads_generator(bam_ncls,
                                                    read_pair_dict,
                                                    qname_dict,
                                                    region[0],
                                                    region[1],
                                                    region[2])

    vertices = defaultdict(list)
    vert_inds = set()
    for read in overlap_reads_iter:
        qname = read.query_name
        if qname in lowqual_qnames or qname not in qname_to_node:
            continue
        vert_idx = qname_to_node[qname]
        vertices[(vert_idx, qname)].append(read)
        vert_inds.add(vert_idx)

    subgraphs = group_by_dict_optimized(qname_hap_info, vertices)

    if len(subgraphs) == 0:
        logger.error(f"No haplotype clusters are found for region {region}. These are the qnames found overlapping this region: {vertices}. Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx

    # For each subgraph, we need to generate the consensus sequence for the reads
    # Each subgraph is a connected component in the phasing graph and it represents a haplotype
    region_haplotype_info = {}
    # overlap_haplotype_info = {}
    inspect_results = []
    for subgraph_id, (subgraph_vertices, qnames, read_pair_lists) in subgraphs.items():
        reads = []
        for read_pair_list in read_pair_lists:
            reads += read_pair_list
        # Sort the reads by start positions
        continuous_covered_regions = extract_continuous_regions_dict(reads)

        for span, sreads in continuous_covered_regions.items():
            if span[0] <= region[1] and span[1] >= region[2]:
                region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), subgraph_id, qnames)
            else:
                inspect_results.append((subgraph_id, qnames, span, sreads, (min(span[1], region[2]) - max(span[0], region[1])) / (region[2] - region[1])))

    if len(region_haplotype_info) <= 2:
        logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region. Try recover some regions.\n")
        inspect_results = sorted(inspect_results, key = lambda x: x[4], reverse = True)
        if len(inspect_results) == 0:
            return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
        elif inspect_results[0][4] < 0.8:
            return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
        else:
            recover_results = [t for t in inspect_results if t[4] >= 0.8]
            recover_span = max(recover_results, key = lambda x: x[2][0])[2][0], min(recover_results, key = lambda x: x[2][1])[2][1]
            overlapping_span = (max(recover_span[0], region[1]), min(recover_span[1], region[2]))
            for t in recover_results:
                subgraph_id, qnames, span, sreads, _ = t
                region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), subgraph_id, qnames)
            if len(region_haplotype_info) <= 2:
                logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region even after the recover trial.\n")
                return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
    else:
        overlapping_span = (region[1], region[2])
    region_str = f"{read.reference_name}:{overlapping_span[0]}-{overlapping_span[1]}"

    final_clusters = {}
    for span, (reads, vert_inds, haplotype_idx, qnames) in region_haplotype_info.items():
        hap_vectors = []
        err_vectors = []
        assert all([qname_hap_info[vid] == haplotype_idx for vid in vert_inds]), f"The vertices {vert_inds} are not in the same connected component."
        if len(reads) == 0:
            logger.warning(f"No reads are found for the haplotype across {span}. Skip this haplotype cluster.")
            continue

        read_spans = np.empty((len(reads), 2), dtype=np.int32)
        hap_vectors = np.full((len(reads), 500), -10, dtype=np.int32)
        err_vectors = np.full((len(reads), 500), -10, dtype=np.float32)
        for i, r in enumerate(reads):
            read_spans[i, 0] = r.reference_start
            read_spans[i, 1] = r.reference_end
            rid = get_read_id(r)
            if rid in total_hapvectors:
                hap_vector = total_hapvectors[rid]
            else:
                hap_vector = get_hapvector_from_cigar(r.cigartuples, r.query_sequence, logger = logger)
                total_hapvectors[rid] = hap_vector

            hap_vectors[i, :hap_vector.size] = hap_vector

            if rid in total_errvectors:
                err_vector = total_errvectors[rid]
            else:
                err_vector = get_errorvector_from_cigar(r, r.cigartuples, logger = logger)
                total_errvectors[rid] = err_vector

            err_vectors[i, :err_vector.size] = err_vector
        # logger.info(f"Haplotype across {span} contains {len(reads)} reads.")
        consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)
        # logger.info(f"The consensus sequence for haplotype {haplotype_idx} across {span} composed of {len(reads)} and the qnames are {[r.query_name for r in reads]}. The consensus sequence is \n{consensus_sequence}\n")
        overlapping_con_seq = consensus_sequence[overlapping_span[0] - span[0]:overlapping_span[1] - span[0] + 1]
        final_clusters[haplotype_idx] = (overlapping_con_seq, reads, overlapping_span, qnames)

    # logger.info(f"The Final haplotype dict looks like this \n{final_clusters}\n")
    if len(final_clusters) <= 2:
        logger.info(f"Only {len(final_clusters)} haplotype clusters are found for region {region_str}. Do not need to choose 2 haplotypes, Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx

    record_2d_arr = record_haplotype_rank( final_clusters, mean_read_length = mean_read_length )
    record_df = pd.DataFrame(record_2d_arr, columns = ["start", "end", "total_depth", "hap_id", "hap_depth", "var_count", "indel_count"])
    record_df["chrom"] = region[0]

    logger.info(f"Found {len(final_clusters)} haplotype clusters for region {region_str}. The dataframe recording the haplotypes in this region looks like :\n{record_df.to_string(index=False)}\n")
    return record_df, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx




def inspect_by_haplotypes(input_bam,
                          hap_qname_info,
                          qname_hap_info,
                          bam_ncls,
                          intrin_bam_ncls,
                          phased_graph,
                          weight_matrix,
                          qname_to_node,
                          total_lowqual_qnames,
                          total_hap_vectors,
                          total_err_vectors,
                          total_genomic_haps,
                          var_df,
                          compare_haplotype_meta_tab = "",
                          mean_read_length = 150,
                          logger = logger):
    ncls_dict, read_dict, qname_dict, qname_idx_dict = bam_ncls
    record_dfs = []
    clique_sep_component_idx = 0
    ref_genome_similarities = defaultdict(dict)
    hid_extreme_vard = defaultdict(bool)
    hid_var_count = defaultdict(int)
    scatter_hid_dict = defaultdict(bool)
    hsid_hid_dict = defaultdict(int)
    hsid_seq_dict = {}
    phased_graph_qname_vprop = phased_graph.vertex_properties['qname']
    logger.info("All the haplotype IDs are :\n{}\n".format(list(hap_qname_info.keys())))

    for hid, qnames in hap_qname_info.items():
        logger.info(f"The haplotype {hid} contains {len(qnames)} read pairs in total.")

        qname_indices = [qname_idx_dict[qname] for qname in qnames]
        reads = [read for qname_idx in qname_indices for read in read_dict[qname_idx]]
        conregion_dict = extract_continuous_regions_dict(reads)

        if len(qnames) < 2:
            scatter_hid_dict[hid] = True

        high_vard = False
        mid_vard_dict = {}
        # Generate consensus sequence for each continuously covered region
        for span, reads in conregion_dict.items():
            # span end is exclusive
            logger.info(f"The haplotype {hid} contains {len(reads)} reads in the continuous region {span}.")
            chrom = reads[0].reference_name
            srids = set([])
            read_spans = np.empty((len(reads), 2), dtype=np.int32)
            hap_vectors = np.full((len(reads), 500), -10, dtype = np.int32)
            err_vectors = np.full((len(reads), 500), -10, dtype = np.float32)
            for i, r in enumerate(reads):
                read_spans[i, 0] = r.reference_start
                read_spans[i, 1] = r.reference_end
                rid = get_read_id(r)
                srids.add(rid)
                if rid in total_hap_vectors:
                    hap_vector = total_hap_vectors[rid]
                else:
                    hap_vector = get_hapvector_from_cigar(r.cigartuples, r.query_sequence, logger = logger)
                    total_hap_vectors[rid] = hap_vector

                hap_vectors[i, :hap_vector.size] = hap_vector

                if rid in total_err_vectors:
                    err_vector = total_err_vectors[rid]
                else:
                    err_vector = get_errorvector_from_cigar(r, r.cigartuples, logger = logger)
                    total_err_vectors[rid] = err_vector

                err_vectors[i, :err_vector.size] = err_vector

            consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)
            extreme_vard = judge_misalignment_by_extreme_vardensity(consensus_sequence)
            var_count = count_var(consensus_sequence)
            hid_var_count[hid] += var_count
            logger.info(f"For haplotype {hid} at continuous region {span}, the consensus sequence is {consensus_sequence.tolist()}. The extreme variant density is {extreme_vard}.")
            hid_extreme_vard[hid] = extreme_vard or hid_extreme_vard[hid]
            # We need to extract the overlapping genomic haplotype vectors
            # Get the genomic haplotype vectors for the overlapping region
            max_similarity = 0
            for hread in overlapping_reads_generator(*intrin_bam_ncls[:-1], chrom, span[0], span[1]):
                hread_start = hread.reference_start
                hread_end = hread.reference_end
                hread_qname = hread.query_name
                overlap_span = (max(hread_start, span[0]), min(hread_end, span[1]))
                # logger.info(f"The genomic sequence spanning {hread_start} to {hread_end}. And it overlaps with the haplotype at {overlap_span}")
                if overlap_span[1] - overlap_span[0] < 140:
                    continue
                hregion_str = f"{chrom}:{hread_start}-{hread_end}"
                hread_id = get_read_id(hread) + ":" + hregion_str
                # hread_cigars = tuple(hread.cigartuples)
                if hread_id in total_genomic_haps:
                    hread_hap_vector = total_genomic_haps[hread_id]
                else:
                    hread_hap_vector = get_hapvector_from_cigar(hread.cigartuples, logger = logger)
                    total_genomic_haps[hread_id] = hread_hap_vector

                try:
                    interval_genomic_hap = hread_hap_vector[overlap_span[0] - hread_start:overlap_span[1] - hread_start]
                except IndexError:
                    continue
                interval_con_seq = consensus_sequence[overlap_span[0] - span[0]:overlap_span[1] - span[0]]
                # logger.info(f"Inspect the conesensus sequence (length {len(consensus_sequence)})\n{consensus_sequence.tolist()}\nwith the genomic haplotype (length {len(interval_genomic_hap)})\n{interval_genomic_hap.tolist()}\n")
                varcount, alt_varcount = ref_genome_similarity(interval_con_seq, interval_genomic_hap)
                if hread_qname in ref_genome_similarities[hid]:
                    ref_genome_similarities[hid][hread_qname].append((varcount, alt_varcount))
                else:
                    ref_genome_similarities[hid][hread_qname] = [(varcount, alt_varcount)]

    ref_genome_str = '\n'.join([f"{k}: {v}" for k,v in ref_genome_similarities.items()])
    logger.info(f"ref genome similarities: {ref_genome_str}")
    logger.info(f"extreme variant density haplotypes: {hid_extreme_vard}")
    logger.info(f"scatter haplotypes: {scatter_hid_dict}")
    final_ref_seq_similarities = {}
    ref_seq_delta = {}
    for hid, gdict in ref_genome_similarities.items():
        max_similarity = 0
        for genomic_seq_qname, pairs in gdict.items():
            total_varcount = np.sum([t[0] for t in pairs])
            total_alt_varcount = np.sum([t[1] for t in pairs])
            similarity = total_varcount / total_alt_varcount if total_alt_varcount > 0 else total_varcount
            delta = total_varcount - total_alt_varcount
            # logger.info(f"For haplotype {hid}, comparing to mapping against the genomic sequence {genomic_seq_qname}. Variant count ratio is {similarity} and variant count delta is {delta}.")
            similarity_score = delta + similarity
            max_similarity = max(max_similarity, similarity_score)

        final_ref_seq_similarities[hid] = max_similarity if max_similarity > 0 else 0

    ref_genome_similarities = final_ref_seq_similarities
    logger.info(f"Final ref sequence similarities for all the haplotypes are {final_ref_seq_similarities}")

    sweep_region_bed = input_bam.replace(".bam", ".sweep.bed")
    sweep_regions = sweep_region_inspection(input_bam, sweep_region_bed, logger = logger)
    logger.info(f"Now the sweep regions are saved to {sweep_region_bed}.")
    for interval in sweep_regions:
        region = (interval.chrom, interval.start, interval.end)
        logger.info(f"Inspecting the region {region}.")
        record_df, total_hap_vectors, total_err_vectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx = identify_misalignment_per_region(region,
                                                                                                                                                        bam_ncls,
                                                                                                                                                        intrin_bam_ncls,
                                                                                                                                                        phased_graph,
                                                                                                                                                        qname_hap_info,
                                                                                                                                                        hap_qname_info,
                                                                                                                                                        weight_matrix,
                                                                                                                                                        qname_to_node,
                                                                                                                                                        total_lowqual_qnames,
                                                                                                                                                        clique_sep_component_idx,
                                                                                                                                                        var_df,
                                                                                                                                                        total_hapvectors = total_hap_vectors,
                                                                                                                                                        total_errvectors = total_err_vectors,
                                                                                                                                                        total_genomic_haps = total_genomic_haps,
                                                                                                                                                        mean_read_length = mean_read_length,
                                                                                                                                                        logger = logger )
        if record_df is not None:
            record_dfs.append(record_df)

    failed_lp = False
    if len(record_dfs) > 0:
        total_record_df = pd.concat(record_dfs)
        record_hapids = set(total_record_df["hap_id"].unique())
        logger.info(f"The haplotypes that have been inspected are {record_hapids}.")
        # Is there any other haplotype that has not been inspected?
        remained_hapids = set(hap_qname_info.keys()) - record_hapids
        if len(remained_hapids) > 0:
            remained_hapid_dict = { k:v for k,v in hap_qname_info.items() if k in remained_hapids }
            logger.info(f"The haplotypes that have not been inspected are {remained_hapids}. Their qnames are: \n{remained_hapid_dict}\n")

        total_record_df["extreme_vard"] = total_record_df["hap_id"].map(hid_extreme_vard)
        total_record_df["scatter_hap"] = total_record_df["hap_id"].map(scatter_hid_dict).fillna(False)
        total_record_df["hap_var_count"] = total_record_df["hap_id"].map(hid_var_count)
        total_record_df["ref_genome_similarities"] = total_record_df["hap_id"].map(ref_genome_similarities).fillna(0)
        # total_record_df.loc[:, "coefficient"] = total_record_df["coefficient"] * 100 + total_record_df.loc[:, "var_count"]
        total_record_df.to_csv(compare_haplotype_meta_tab.replace(".tsv", ".raw.tsv"), sep = "\t", index = False)
        logger.info(f"Successfully saved the raw haplotype comparison meta table to {compare_haplotype_meta_tab.replace('.tsv', '.raw.tsv')}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        total_record_df = total_record_df.loc[np.logical_not(total_record_df["extreme_vard"]) & \
                                              np.logical_not(total_record_df["scatter_hap"]) & \
                                              (total_record_df["total_depth"] >= 5) & \
                                              (total_record_df["ref_genome_similarities"] <= 5), :]
        # total_record_df.loc[:, "coefficient"] = np.clip(total_record_df["coefficient"] + total_record_df["ref_genome_similarities"], 10e-3, None)
        if total_record_df.shape[0] == 0:
            failed_lp = True

        if not failed_lp:
            total_record_df.loc[total_record_df["hap_var_count"] == 0, "coefficient"] = -1
            total_record_df = total_record_df.groupby(["chrom", "start", "end"]).filter(lambda x: len(x) > 5)

            if total_record_df.shape[0] == 0:
                failed_lp = True
    else:
        failed_lp = True

    if not failed_lp:
        by_region = total_record_df.groupby(["chrom", "start", "end"], group_keys=False)
        total_record_df = by_region.apply(calulate_coefficient_per_group, logger = logger).reset_index(drop=True)
        total_record_df.to_csv(compare_haplotype_meta_tab, sep = "\t", index = False)
        logger.info(f"Successfully saved the haplotype comparison meta table to {compare_haplotype_meta_tab}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        correct_map_hids, mismap_hids, model_status = lp_solve_remained_haplotypes(total_record_df, logger = logger)
        if model_status == "Infeasible":
            failed_lp = True

    if not failed_lp:
        mismap_hids.update(set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]]))
        mismap_hids.update(set([hid for hid in ref_genome_similarities if ref_genome_similarities[hid] > 5]))
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] > 1]))
        correct_map_hids = correct_map_hids - mismap_hids

        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        correct_map_qnames = set([qname for hid in correct_map_hids for qname in hap_qname_info[hid]])
        logger.info(f"Identify {len(mismap_hids)} mismapped haplotypes, {len(mismap_qnames)} mismapped qnames by solving the BINARY INTEGER LINEAR PROGRAMMING.\n")
        logger.info(f"Here is the qnames seemed mismapped: \n{mismap_qnames}\nAnd the mismap haplotype IDs: \n{mismap_hids}\n")
        logger.info(f"Here is the qnames seemed correctly mapped: \n{correct_map_qnames}\nAnd the correct map haplotype IDs: \n{correct_map_hids}\n")
        return correct_map_qnames, mismap_qnames

    if failed_lp:
        mismap_hids = set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]])
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] > 1]))
        mismap_hids.update(set([hid for hid in ref_genome_similarities if ref_genome_similarities[hid] > 5]))
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        total_qnames = set([qname for qnames in hap_qname_info.values() for qname in qnames])
        correct_map_qnames = total_qnames - mismap_qnames
        return correct_map_qnames, mismap_qnames