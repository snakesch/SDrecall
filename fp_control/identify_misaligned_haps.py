import numba
import numpy as np
import pandas as pd
import pybedtools as pb

from collections import defaultdict
from numba import types, prange

from fp_control.bam_ncls import overlapping_reads_iterator
from fp_control.bilc import lp_solve_remained_haplotypes
from fp_control.numba_operators import numba_sum
from src.log import logger
from src.utils import executeCmd, prepare_tmp_file
from fp_control.pairwise_read_inspection import get_hapvector_from_cigar, \
                                                get_errorvector_from_cigar, \
                                                get_read_id, \
                                                count_var, \
                                                count_continuous_blocks, \
                                                count_continuous_indel_blocks, \
                                                extract_read_qseqs



@numba.njit(types.Tuple((types.int32, types.int32))(types.int16[:], types.int16[:]), fastmath=True)
def ref_genome_similarity(query_read_vector,
                          genomic_hap_vector):
    '''
    This function is used to identify the haplotypes that the query read vector belongs to.
    The genomic_hap_vectors is a 2D numpy array
    The query_read_vector is a 1D numpy array
    The max_vector_distance is a float reflecting the max manhattan distance between a read and a genomic haplotype
    '''
    # assert len(genomic_hap_vectors.shape) == 2, f"The input genomic haplotype vectors should be a 2D numpy array, but the actual input is {genomic_hap_vectors}"
    # Remove the reference haplotype first
    if numba_sum(genomic_hap_vector != 1) == 0:
        return 0, 0

    # ref_read_vector = np.ones(query_read_vector.size, dtype=np.int32)
    # Calculate the Manhattan distance between the query read vector and each genomic haplotype
    alt_var_size = count_continuous_blocks(genomic_hap_vector != query_read_vector)
    var_size = count_var(query_read_vector)

    return var_size, alt_var_size



@numba.njit
def rank_unique_values(arr):
    ## Extract unique values and sort them
    unique_values = np.unique(arr)

    ## Rank unique values by sorted order
    ranks = np.empty(arr.shape, dtype=np.int32)
    for i, it in enumerate(arr):
        for j, jt in enumerate(unique_values):
            if it == jt:
                ranks[i] = j + 1
                break
    
    return ranks


@numba.njit(types.int16[:](types.int16[:,:], types.float32[:,:], types.int32[:,:]), fastmath=True)
def assemble_consensus(seq_arrays, qual_arrays, read_spans):
    # Every read in the reads can have different length
    # Find the start and end position of the consensus sequence
    start_pos = read_spans[:, 0].min()
    end_pos = read_spans[:, 1].max()

    # logger.info(", ".join([f"{r.reference_start}-{r.reference_end}" for r in reads]))

    # Initialize the consensus sequence and quality scores with zeros
    consensus_seq = np.ones(end_pos - start_pos, dtype=np.int16)
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
        read_span = read_spans[i]
        start = read_span[0]
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




@numba.njit(types.float32[:](types.int16[:], types.int32), fastmath=True)
def count_window_var_density(array, padding_size = 25):
    is_var = array != 1
    is_var_size = numba_sum(is_var)
    if is_var_size == 0:
        return np.zeros(array.size, dtype=np.float32)

    density_arr = np.empty(array.size, dtype = np.float32)
    for i in range(array.size):
        start = max(i-padding_size, 0)
        end = min(i+padding_size, array.size)
        iter_arr = array[start:end]
        var_count = count_var(iter_arr)
        density = var_count/(padding_size *2 + 1)
        density_arr[i] = density

    return density_arr


@numba.njit(types.int32[:,:](types.boolean[:]), fastmath=True)
def extract_true_stretches(bool_array):
    n = len(bool_array)

    # Pre-allocate maximum possible space (worst case: every element is True)
    stretches = np.zeros((n, 2), dtype=np.int32)
    count = 0

    if n == 0:
        return stretches[:count]

    in_stretch = False
    start_idx = 0

    for i in range(n):
        if bool_array[i]:
            if not in_stretch:
                in_stretch = True
                start_idx = i
        else:
            if in_stretch:
                stretches[count, 0] = start_idx
                stretches[count, 1] = i - 1
                count += 1
                in_stretch = False

    # Handle the case where the array ends with a True stretch
    if in_stretch:
        stretches[count, 0] = start_idx
        stretches[count, 1] = n - 1
        count += 1

    return stretches[:count]



@numba.njit(types.boolean(types.int16[:]), fastmath=True)
def judge_misalignment_by_extreme_vardensity(seq):
    five_vard = count_window_var_density(seq, padding_size = 42)
    six_vard = count_window_var_density(seq, padding_size = 65)
    read_vard = count_window_var_density(seq, padding_size = 74)
    
    if numba_sum(five_vard >= 7/85) > 0:
        select_bool = five_vard >= 7/85
        padding = 42
    elif numba_sum(six_vard >= 9/131) > 0:
        select_bool = six_vard >= 9/131
        padding = 65
    elif numba_sum(read_vard > 11/148) > 0:
        return True
    else:
        return False

    true_segments = extract_true_stretches(select_bool)
    max_indel_count = 0
    for i in range(true_segments.shape[0]):
        start = true_segments[i, 0] - padding
        end = true_segments[i, 1] + padding
        five_seq = seq[start:end]
        indel_count = count_continuous_indel_blocks(five_seq)
        max_indel_count = max(max_indel_count, indel_count)
    if max_indel_count > 1:
        return True
    return False



@numba.njit(types.float32[:](types.int32[:, :]), fastmath=True)
def calculate_coefficient(arr2d):
    # Span size is the weight for later mean value calculation
    rank = arr2d[:, 7]
    span = (arr2d[:, 1] - arr2d[:, 0])
    depth_frac = (1 - arr2d[:, 4]/arr2d[:, 2]).astype(types.float32)

    return rank * span * np.sqrt(depth_frac)



def sweep_region_inspection(hap_cov_beds,
                            output_bed = None,
                            logger = logger):

    if output_bed is None:
        output_bed = prepare_tmp_file(suffix = ".bed").name

    hap_ids, hap_id_beds = zip(*hap_cov_beds.items())
    x = pb.BedTool()
    sweep_regions = x.multi_intersect(i=hap_id_beds).to_dataframe(disable_auto_names = True, names = ["chrom", "start", "end", "hap_no", "hap_ids"] + [f"unknown_{i}" for i in range(len(hap_id_beds))]).iloc[:, :4]
    multi_haps_regions = sweep_regions.loc[sweep_regions["hap_no"].astype(int) >= 3, :]
    if len(multi_haps_regions) == 0:
        logger.warning(f"No sweep regions are found. Return None. The original sweep regions are:\n{sweep_regions.to_string(index=False)}")
        return None
    multi_haps_regions.to_csv(output_bed, sep = "\t", index = False, header = False)

    op_bed_obj = pb.BedTool(output_bed)
    return op_bed_obj




def calculate_coefficient_per_group(record_df, logger=logger):
    # Generate a rank for haplotypes
    record_df["rank"] = rank_unique_values(record_df["var_count"].to_numpy() * 50 + record_df["indel_count"].to_numpy())
    # logger.info(f"Before calculating the coefficient for this region, the dataframe looks like :\n{record_df[:10].to_string(index=False)}\n")
    record_df["coefficient"] = calculate_coefficient(record_df.loc[:, ["start",
                                                                       "end",
                                                                       "total_depth",
                                                                       "hap_id",
                                                                       "hap_depth",
                                                                       "var_count",
                                                                       "indel_count",
                                                                       "rank"]].to_numpy(dtype=np.int32))
    # logger.info(f"After calculating the coefficient for this region, the dataframe looks like :\n{record_df[:10].to_string(index=False)}\n")
    return record_df



def extract_continuous_regions_dict(reads):
    # Sort reads by reference start position
    reads = sorted(reads, key=lambda read: read.reference_start)

    continuous_regions = {}
    current_region_reads = []
    current_start = None
    current_end = None

    for read in reads:
        read_start = read.reference_start # 0 indexed
        read_end = read.reference_end  # 0 indexed and refers to the one past the last aligned base

        if current_start is None:
            # Initialize the first region
            current_start = read_start
            current_end = read_end
            current_region_reads.append(read)
        else:
            if read_start <= current_end and read_end >= current_start:
                # Extend the current region if the iterating read overlaps with it.
                current_start = min(current_start, read_start)
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



def group_by_dict_optimized(vprop, vertices):
    grouped_keys = {}
    for v, rs in vertices.items():
        v_idx, q = v
        label = vprop[v_idx]
        if label not in grouped_keys:
            grouped_keys[label] = [[], [], []]
        grouped_keys[label][0].append(v_idx)
        grouped_keys[label][1].append(q)
        grouped_keys[label][2].append(rs)
    return grouped_keys



def record_haplotype_rank(haplotype_dict, mean_read_length = 150):
    starts = np.empty(len(haplotype_dict), dtype=np.int32)
    ends = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_depths = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_ids = np.empty(len(haplotype_dict), dtype=np.int32)
    var_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    indel_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    total_depth_col = np.empty(len(haplotype_dict), dtype=np.int32)

    i = 0
    total_depth = 0
    for hid in haplotype_dict:
        seq, reads, span, _ = haplotype_dict[hid]
        starts[i] = span[0]
        ends[i] = span[1]
        depth = len(reads) * mean_read_length/len(seq)
        hap_depths[i] = depth
        hap_ids[i] = hid
        var_counts[i] = count_var(seq)
        indel_counts[i] = count_continuous_indel_blocks(seq)
        total_depth += depth
        i += 1

    total_depth_col.fill(total_depth)
    return np.column_stack((starts, ends, total_depth_col, hap_ids, hap_depths, var_counts, indel_counts))



def summarize_enclosing_haps(hap_subgraphs, qname_to_node, region, logger = logger):
    '''
    For each haplotype cluster, we need to extract the continuous regions that are covered by the reads in the haplotype.
    Each haplotype can contain multiple continuous regions
    '''
    _, start, end = region
    region_haplotype_info = {} # {(start, end): (reads, vert_inds, hap_id, qnames)}
    inspect_results = [] # [(hap_id, qnames, span, sreads, overlap_coef)]
    for hap_id, (_, qnames, read_pair_lists) in hap_subgraphs.items():
        reads = []
        for read_pair_list in read_pair_lists:
            reads += read_pair_list

        # Sort the reads by start positions and return a dictionary of the form {(start, end): [reads]}
        continuous_covered_regions = extract_continuous_regions_dict(reads)

        for span, sreads in continuous_covered_regions.items():
            # Find enclosing haplotype covered regions
            if span[0] <= start and span[1] >= end:
                vert_inds = set([qname_to_node[r.query_name] for r in sreads])
                region_haplotype_info[span] = (sreads, vert_inds, hap_id, qnames)
            else:
                # The continuous region is not fully enclosing the current window
                # We can try to rescue these regions
                overlap_coef = (min(span[1], end) - max(span[0], start)) / (end - start)
                inspect_results.append((hap_id, qnames, span, sreads, overlap_coef))

    if len(region_haplotype_info) > 2:
        overlapping_span = (start, end)
        return region_haplotype_info, overlapping_span

    '''
    Now we dont have enough haplotypes enclosing the inspection window.
    So we try to re-adjust the window size (shrink a little bit) and see if we can get more haplotypes enclosing the inspection window.
    '''

    logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region. Try recover some regions.\n")
    inspect_results = sorted(inspect_results, key = lambda x: x[4], reverse = True)
    if len(inspect_results) == 0:
        # No else haplotypes are found
        return None, None
    elif inspect_results[0][4] < 0.8:
        # Even the largest overlapping haplotype only cover 80% of the inspection window
        # Give up the window
        return None, None
    else:
        # Shrink the window to the largest overlapping haplotype
        recover_results = [t for t in inspect_results if t[4] >= 0.8]
        recover_span = max(recover_results, key = lambda x: x[2][0])[2][0], min(recover_results, key = lambda x: x[2][1])[2][1]
        overlapping_span = (max(recover_span[0], start), min(recover_span[1], end))
        for t in recover_results:
            hap_id, qnames, span, sreads, _ = t
            region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), hap_id, qnames)
        if len(region_haplotype_info) <= 2:
            logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region even after the recover trial.\n")
            return None, None

    return region_haplotype_info, overlapping_span




def identify_misalignment_per_region(region,
                                     bam_ncls,
                                     qname_hap_info,
                                     qname_to_node,
                                     lowqual_qnames,
                                     clique_sep_component_idx,
                                     total_hapvectors = {},
                                     total_errvectors = {},
                                     total_genomic_haps = {},
                                     read_ref_pos_dict = {},
                                     mean_read_length = 148,
                                     logger = logger ):
    bam_ncls, read_pair_dict, *_ = bam_ncls
    '''
    Rank haplotypes enclosed by region by sequence similarity to paralogous reference sequences

    Inputs:
    region: a tuple of (chrom, start, end)
    bam_ncls: a tuple of ({chrom: ncls_object}, {interval_idx: [read1, read2]}, {interval_idx: qname}, {qname: interval_idx})
    qname_hap_info: a dictionary {qname: haplotype_id}
    qname_to_node: a dictionary {qname: vertex_idx}
    lowqual_qnames: a set of qnames that are low quality and ignored
    clique_sep_component_idx: a dictionary {haplotype_id: clique_separated_component_idx}
    total_hapvectors: a dictionary {read_id: haplotype_vector}
    total_errvectors: a dictionary {read_id: error_vector}
    total_genomic_haps: a dictionary {reference_sequence_id: haplotype_vector}
    '''
    ## Identify reads overlapping region
    chrom, start, end = region
    overlap_reads_iter = overlapping_reads_iterator( bam_ncls,
                                                     read_pair_dict,
                                                     chrom,
                                                     start,
                                                     end )

    # Identify a map that
    # (vertex_idx, qname) --> read object (pysam.AlignedSegment)
    vertices = defaultdict(list)
    # a set of vertex indices (qname idx) corresponding to the overlapping reads
    vert_inds = set()
    for read in overlap_reads_iter:
        qname = read.query_name
        if qname in lowqual_qnames or qname not in qname_to_node:
            continue
        vert_idx = qname_to_node[qname]
        vertices[(vert_idx, qname)].append(read)
        vert_inds.add(vert_idx)

    # qname_hap_info is a dictionary {qname: haplotype_id}
    # Establish a dictionary of the form {haplotype_id: [vertex_idices, qnames, read_pair_lists]} based on the overlapping reads identified above
    hap_subgraphs = group_by_dict_optimized(qname_hap_info, vertices)

    if len(hap_subgraphs) == 0:
        logger.error(f"No haplotype clusters are found for region {region}. These are the qnames found overlapping this region: {vertices}. Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict

    # For each hap_subgraph, we need to generate the consensus sequence for the reads
    # Each hap_subgraph is a connected component in the phasing graph and it represents a haplotype
    region_haplotype_info, overlapping_span = summarize_enclosing_haps(hap_subgraphs, qname_to_node, region, logger = logger)
    if region_haplotype_info is None:
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict

    # Prepare region_str for logging
    region_str = f"{read.reference_name}:{overlapping_span[0]}-{overlapping_span[1]}"

    final_clusters = {}
    for span, (reads, vert_inds, haplotype_idx, qnames) in region_haplotype_info.items():
        hap_vectors = []
        err_vectors = []
        assert all([qname_hap_info[vid] == haplotype_idx for vid in vert_inds]), f"The vertices {vert_inds} are not in the same connected component."
        if len(reads) == 0:
            logger.warning(f"No reads are found for the haplotype across {span}. Skip this haplotype cluster.")
            continue

        read_spans, hap_vectors, err_vectors, total_hapvectors, total_errvectors, read_ref_pos_dict = record_hap_err_vectors_per_region(reads,
                                                                                                                                        total_hapvectors,
                                                                                                                                        total_errvectors,
                                                                                                                                        read_ref_pos_dict = read_ref_pos_dict,
                                                                                                                                        logger = logger)

        # logger.info(f"Haplotype across {span} contains {len(reads)} reads.")
        consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)
        # logger.info(f"The consensus sequence for haplotype {haplotype_idx} across {span} composed of {len(reads)} and the qnames are {[r.query_name for r in reads]}. The consensus sequence is \n{consensus_sequence}\n")
        overlapping_con_seq = consensus_sequence[overlapping_span[0] - span[0]:overlapping_span[1] - span[0] + 1]
        final_clusters[haplotype_idx] = (overlapping_con_seq, reads, overlapping_span, qnames)

    # logger.info(f"The Final haplotype dict looks like this \n{final_clusters}\n")
    if len(final_clusters) <= 2:
        logger.warning(f"Only {len(final_clusters)} haplotype clusters are found for region {region_str}. Do not need to choose 2 haplotypes, Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict

    record_2d_arr = record_haplotype_rank( final_clusters, mean_read_length = mean_read_length )
    record_df = pd.DataFrame(record_2d_arr, columns = ["start", "end", "total_depth", "hap_id", "hap_depth", "var_count", "indel_count"])
    record_df["chrom"] = chrom

    logger.info(f"Found {len(final_clusters)} haplotype clusters for region {region_str}. The dataframe recording the haplotypes in this region looks like :\n{record_df.to_string(index=False)}\n")
    return record_df, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict



def record_hap_err_vectors_per_region(reads, 
                                      total_hap_vectors, 
                                      total_err_vectors, 
                                      read_ref_pos_dict = {},
                                      logger = logger):
    # Initialize the read spans (2d array), haplotype vectors (2d array) and error vectors (2d array)
    read_spans = np.empty((len(reads), 2), dtype=np.int32)  # 2d array to store the start and end positions of the reads
    hap_vectors = np.full((len(reads), 500), -10, dtype = np.int16)  # 2d array to store the haplotype vectors (every row stores a haplotype vector, usually haplotype vector is shorter than 500, the remained positions are filled with -10)
    err_vectors = np.full((len(reads), 500), -10, dtype = np.float32)  # 2d array to store the error vectors (every row stores a error vector, usually error vector is shorter than 500, the remained positions are filled with -10)

    # Initialize a set to store the read IDs of the reads overlapping with the iterating continuous region
    srids = set()
    # Iterate over all the reads overlapping with the iterating continuous region
    for i, r in enumerate(reads):
        read_spans[i, 0] = r.reference_start
        read_spans[i, 1] = r.reference_end
        rid = get_read_id(r)
        srids.add(rid)

        # Check if the hap vector of the read has been computed and stored in the total_hap_vectors dictionary, a short-circuit evaluation to avoid redundant computation of function get_hapvector_from_cigar
        if rid in total_hap_vectors:
            hap_vector = total_hap_vectors[rid]
        else:
            _, _, query_sequence_encoded, _, read_ref_pos_dict = extract_read_qseqs(r, read_ref_pos_dict)
            cigar_arr = np.array(r.cigartuples, dtype = np.int32)
            hap_vector = get_hapvector_from_cigar(cigar_arr, query_sequence_encoded)
            total_hap_vectors[rid] = hap_vector

        # Store the haplotype vector of the read into the hap_vectors 2d array
        hap_vectors[i, :hap_vector.size] = hap_vector

        # Check if the error vector of the read has been computed and stored in the total_err_vectors dictionary, a short-circuit evaluation to avoid redundant computation of function get_errorvector_from_cigar
        if rid in total_err_vectors:
            err_vector = total_err_vectors[rid]
        else:
            err_vector = get_errorvector_from_cigar(r, r.cigartuples, logger = logger)
            total_err_vectors[rid] = err_vector

        # Store the error vector of the read into the err_vectors 2d array
        err_vectors[i, :err_vector.size] = err_vector

    return read_spans, hap_vectors, err_vectors, total_hap_vectors, total_err_vectors, read_ref_pos_dict



def stat_refseq_similarity(intrin_bam_ncls,
                           chrom,
                           span,
                           hid,
                           consensus_sequence,
                           total_genomic_haps,
                           read_ref_pos_dict,
                           varcounts_among_refseqs,
                           logger = logger):
    '''
    For each haplotype covered region (continuous region),
    We have a consensus sequence of the haplotype assembled by reads realigned here from some the homologous regions.
    We also have reference sequences aligned to the query region.

    If we want to measure if the haplotype is correctly aligned here, we could compare the below 2 metrics:
      1. the similarity (reversely indicated by variant count) between the haplotype and the reference sequence
      2. the similarity (reversely indicated by variant count) between the haplotype and the homologous counterparts.

    if 2 is much larger than 1, it definitely increase the likelihood that the haplotype is misaligned.
    '''

    intrin_ncls_dict, intrin_read_dict, _, _ = intrin_bam_ncls
    for homo_refseq in overlapping_reads_iterator(intrin_ncls_dict, intrin_read_dict, chrom, span[0], span[1]):
        homo_refseq_start = homo_refseq.reference_start
        homo_refseq_end = homo_refseq.reference_end
        homo_refseq_qname = homo_refseq.query_name
        overlap_span = (max(homo_refseq_start, span[0]), min(homo_refseq_end, span[1]))
        # logger.debug(f"The genomic sequence spanning {homo_refseq_start} to {homo_refseq_end}. And it overlaps with the haplotype at {overlap_span}")

        hregion_str = f"{chrom}:{homo_refseq_start}-{homo_refseq_end}"
        # Some reference sequences can be mapped to multiple places, to have a unique ID, we need to append the current region coordinate string
        homo_refseq_id = get_read_id(homo_refseq) + ":" + hregion_str

        # Avoid redundant computation of function get_hapvector_from_cigar
        if homo_refseq_id in total_genomic_haps:
            homo_refseq_hap_vector = total_genomic_haps[homo_refseq_id]
        else:
            _, _, query_sequence_encoded, _, read_ref_pos_dict = extract_read_qseqs(homo_refseq, read_ref_pos_dict)
            cigar_arr = np.array(homo_refseq.cigartuples, dtype = np.int32)
            homo_refseq_hap_vector = get_hapvector_from_cigar(cigar_arr, query_sequence_encoded)
            total_genomic_haps[homo_refseq_id] = homo_refseq_hap_vector

        # Extract the genomic haplotype vector overlapping with the iterating continuous region
        try:
            interval_genomic_hap = homo_refseq_hap_vector[overlap_span[0] - homo_refseq_start:overlap_span[1] - homo_refseq_start]
        except IndexError:
            continue

        # Extract the consensus sequence overlapping with the iterating continuous region
        interval_con_seq = consensus_sequence[overlap_span[0] - span[0]:overlap_span[1] - span[0]]
        # logger.info(f"Inspect the conesensus sequence (length {len(consensus_sequence)})\n{consensus_sequence.tolist()}\nwith the genomic haplotype (length {len(interval_genomic_hap)})\n{interval_genomic_hap.tolist()}\n")

        # Identify paralogous sequence variants (PSVs) originated from paralogous reference sequences
        varcount, alt_varcount = ref_genome_similarity(interval_con_seq, interval_genomic_hap)
        if homo_refseq_qname in varcounts_among_refseqs[hid]:
            varcounts_among_refseqs[hid][homo_refseq_qname].append((varcount, alt_varcount))
        else:
            varcounts_among_refseqs[hid][homo_refseq_qname] = [(varcount, alt_varcount)]

    return varcounts_among_refseqs, read_ref_pos_dict



def cal_similarity_score(varcounts_among_refseqs):
    '''
    Input is a dictionary with two layers of keys:
    haplotype_id --> reference_sequence_id --> (variant_count, alt_variant_count) where:
    variant_count is the variant count when the haplotype is aligned to the reference sequence
    alt_variant_count is the variant count when the haplotype is aligned to the homologous counterparts of current region

    For each haplotype, we compare its consensus sequence with multiple reference sequences from a group of homologous sequences,
    And calculate a metric calles similarity score for each comparison.

    Similarity Score = Ratio + Delta
    Ratio = variant_count / alt_variant_count if alt_variant_count > 0 else variant_count, the bigger it is, the more likely the haplotype is misaligned
    Delta = variant_count - alt_variant_count, the bigger it is, the more likely the haplotype is misaligned

    For each haplotype, pick the maximum similarity score across all the reference sequences from homologous counterparts.
    Return a dictionary with only one layer of keys:
    haplotype_id --> max_similarity_score
    '''
    # Initialize a dictionary to store the maximum similarity between the haplotype and the reference sequence from homologous counterparts
    hap_max_sim_scores = {}
    # For each haplotype, find the maximum similarity delta.

    # Iterate over all the reference genome similarities for all the haplotypes
    for hid, gdict in varcounts_among_refseqs.items():
        max_varcount_delta = 0
        for _, pairs in gdict.items():
            total_varcount = np.sum([t[0] for t in pairs])
            total_alt_varcount = np.sum([t[1] for t in pairs])
            similarity = total_varcount / total_alt_varcount if total_alt_varcount > 0 else total_varcount
            delta = total_varcount - total_alt_varcount
            # logger.info(f"For haplotype {hid}, comparing to mapping against the genomic sequence {genomic_seq_qname}. Variant count ratio is {similarity} and variant count delta is {delta}.")
            similarity_score = delta + similarity
            max_varcount_delta = max(max_varcount_delta, similarity_score)

        hap_max_sim_scores[hid] = max_varcount_delta if max_varcount_delta > 0 else 0

    return hap_max_sim_scores



def inspect_by_haplotypes(input_bam,
                          hap_qname_info,
                          qname_hap_info,
                          bam_ncls,
                          intrin_bam_ncls,
                          qname_to_node,
                          total_lowqual_qnames,
                          total_hap_vectors,
                          total_err_vectors,
                          total_genomic_haps,
                          read_ref_pos_dict,
                          compare_haplotype_meta_tab = "",
                          mean_read_length = 150,
                          tmp_dir = "/tmp",
                          logger = logger):
    _, read_dict, _, qname_idx_dict, _ = bam_ncls
    record_dfs = []
    clique_sep_component_idx = 0
    hid_extreme_vard = defaultdict(bool) # Initialize a dictionary to store if the haplotype has extreme variant density
    hid_var_count = defaultdict(int) # Initialize a dictionary to store the variant count of the haplotype
    scatter_hid_dict = defaultdict(bool) # Initialize a dictionary to store if the haplotype is scattered
    logger.info("All the haplotype IDs are :\n{}\n".format(list(hap_qname_info.keys())))
    varcounts_among_refseqs = defaultdict(dict)
    hid_cov_beds = {}

    # Iterate over all the haplotypes
    for hid, qnames in hap_qname_info.items():
        logger.debug(f"The haplotype {hid} contains {len(qnames)} read pairs in total.")
        hid_cov_bed = prepare_tmp_file(suffix = f".haplotype_{hid}.cov.bed", tmp_dir = tmp_dir).name

        # Extract all the qname indices and all the reads belonged to the iterating haplotype
        qname_indices = [qname_idx_dict[qname] for qname in qnames]
        reads = [read for qname_idx in qname_indices for read in read_dict[qname_idx]]

        # Extract all the continuous regions covered by the iterating haplotype
        conregion_dict = extract_continuous_regions_dict(reads)
        # Write the continuous regions to the temporary file
        region_str = ""
        for span, reads in conregion_dict.items():
            region_str += f"{reads[0].reference_name}\t{span[0]}\t{span[1]}\n"

        pb.BedTool(region_str, from_string = True).sort().merge().saveas(hid_cov_bed)
        hid_cov_beds[hid] = hid_cov_bed

        # If only one read pair is in the iterating haplotype, it is a scattered haplotype, unless the total involving read pairs are small
        if len(qnames) < 2 and len(read_dict) > 20:
            scatter_hid_dict[hid] = True

        # Iterate over all the continuous regions covered by the iterating haplotype
        for span, reads in conregion_dict.items():
            # span end is exclusive
            logger.debug(f"The haplotype {hid} contains {len(reads)} reads in the continuous region {span}.")
            chrom = reads[0].reference_name

            # Iterate over all the reads overlapping with the iterating continuous region covered by the iterating haplotype
            # Record the haplotype vectors and error vectors for all the reads overlapping with the iterating continuous region to the prepared 2d arrays
            read_spans, hap_vectors, err_vectors, total_hap_vectors, total_err_vectors, read_ref_pos_dict = record_hap_err_vectors_per_region(reads,
                                                                                                                                              total_hap_vectors,
                                                                                                                                              total_err_vectors,
                                                                                                                                              read_ref_pos_dict = read_ref_pos_dict,
                                                                                                                                              logger = logger)

            # Assemble the consensus sequence for the iterating continuous region
            consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)

            # Judge if the consensus sequence of the haplotype within the iterating continuous region contains extremely high variant density
            extreme_vard = judge_misalignment_by_extreme_vardensity(consensus_sequence)

            # Count the variant count of the consensus sequence
            var_count = count_var(consensus_sequence)

            # Update the variant count of the iterating haplotype
            hid_var_count[hid] += var_count

            # Log the extreme variant density and the variant count of the iterating haplotype
            logger.debug(f"For haplotype {hid} at continuous region {span}, the consensus sequence is {consensus_sequence.tolist()}. The extreme variant density is {extreme_vard}.")

            # Update the extreme variant density of the iterating haplotype
            hid_extreme_vard[hid] = extreme_vard or hid_extreme_vard[hid]

            # Within the iterating continuous region, extract the overlapping genomic haplotype vectors (reference sequences elsewhere mapped to current place)
            # The returned varcounts_among_refseqs is a dictionary, it has two layers of keys:
            # haplotype_id --> reference_sequence_id --> (variant_count, alt_variant_count)
            # Where variant_count is the variant count when the haplotype is aligned to the reference sequence
            # And alt_variant_count is the variant count when the haplotype is aligned to the homologous counterparts of current region
            varcounts_among_refseqs, read_ref_pos_dict = stat_refseq_similarity(intrin_bam_ncls,
                                                                                chrom,
                                                                                span,
                                                                                hid,
                                                                                consensus_sequence,
                                                                                total_genomic_haps,
                                                                                read_ref_pos_dict,
                                                                                varcounts_among_refseqs,
                                                                                logger = logger)
    # Log the extreme variant density haplotypes
    logger.info(f"extreme variant density haplotypes: {hid_extreme_vard}")
    # Log the scattered haplotypes
    logger.info(f"scatter haplotypes: {scatter_hid_dict}")

    # Calculate the maximum similarity score for all the haplotypes
    # Detailed explanation of the similarity score is in the docstring of the function cal_similarity_score
    hap_max_sim_scores = cal_similarity_score(varcounts_among_refseqs)
    # Log the final reference sequence similarities for all the haplotypes
    logger.info(f"Final ref sequence similarities for all the haplotypes are {hap_max_sim_scores}")

    sweep_region_bed = input_bam.replace(".bam", ".sweep.bed")
    # Remove the scattered haplotypes from the sweep regions
    hid_cov_beds = {hid:bed for hid, bed in hid_cov_beds.items() if hid not in scatter_hid_dict}
    sweep_regions = sweep_region_inspection(hid_cov_beds, sweep_region_bed, logger = logger)
    logger.info(f"Now the sweep regions are saved to {sweep_region_bed}.")

    '''
    The below iteration is to inspect haplotypes by each window decided by the sweep_region_inspection function.
    Within each window, we rank the enclosing haplotypes according to their similarity with the reference genomic sequence.
    '''
    if sweep_regions is None:
        mismap_hids = set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]])
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        total_qnames = set([qname for qnames in hap_qname_info.values() for qname in qnames])
        correct_map_qnames = total_qnames - mismap_qnames
        logger.warning(f"No regions are found encompassed by more than 2 haplotypes. So we only filter out the haplotypes with extreme variant density. Filtered out {len(mismap_qnames)} read pairs and left {len(correct_map_qnames)} read pairs.")
        logger.info(f"Here is the qnames seemed mismapped: \n{mismap_qnames}\nAnd the mismap haplotype IDs: \n{mismap_hids}\n")
        return correct_map_qnames, mismap_qnames
    
    # Iterate over all the windows
    for interval in sweep_regions:
        region = (interval.chrom, interval.start, interval.end)
        logger.info(f"Inspecting the region {region}.")
        record_df, total_hap_vectors, total_err_vectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict = identify_misalignment_per_region(region,
                                                                                                                                                                            bam_ncls,
                                                                                                                                                                            qname_hap_info,
                                                                                                                                                                            qname_to_node,
                                                                                                                                                                            total_lowqual_qnames,
                                                                                                                                                                            clique_sep_component_idx,
                                                                                                                                                                            total_errvectors = total_err_vectors,
                                                                                                                                                                            total_genomic_haps = total_genomic_haps,
                                                                                                                                                                            read_ref_pos_dict = read_ref_pos_dict,
                                                                                                                                                                            mean_read_length = mean_read_length,
                                                                                                                                                                            logger = logger )
        if record_df is not None:
            record_dfs.append(record_df)

    # Below we need to consider the cases BILC failed to converge.
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

        # Convert the dictionary map to the columns in the dataframe
        total_record_df["extreme_vard"] = total_record_df["hap_id"].map(hid_extreme_vard)
        total_record_df["scatter_hap"] = total_record_df["hap_id"].map(scatter_hid_dict).fillna(False)
        total_record_df["hap_var_count"] = total_record_df["hap_id"].map(hid_var_count)
        total_record_df["hap_max_sim_scores"] = total_record_df["hap_id"].map(hap_max_sim_scores).fillna(0)
        # total_record_df.loc[:, "coefficient"] = total_record_df["coefficient"] * 100 + total_record_df.loc[:, "var_count"]
        total_record_df.to_csv(compare_haplotype_meta_tab.replace(".tsv", ".raw.tsv"), sep = "\t", index = False)
        logger.info(f"Successfully saved the raw haplotype comparison meta table to {compare_haplotype_meta_tab.replace('.tsv', '.raw.tsv')}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        total_record_df = total_record_df.loc[np.logical_not(total_record_df["extreme_vard"]) & \
                                              np.logical_not(total_record_df["scatter_hap"]) & \
                                              (total_record_df["total_depth"] >= 5) & \
                                              (total_record_df["hap_max_sim_scores"] <= 10), :]
        # total_record_df.loc[:, "coefficient"] = np.clip(total_record_df["coefficient"] + total_record_df["hap_max_sim_scores"], 10e-3, None)
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
        total_record_df = by_region.apply(calculate_coefficient_per_group, logger = logger).reset_index(drop=True)
        total_record_df.to_csv(compare_haplotype_meta_tab, sep = "\t", index = False)
        logger.info(f"Successfully saved the haplotype comparison meta table to {compare_haplotype_meta_tab}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        correct_map_hids, mismap_hids, model_status = lp_solve_remained_haplotypes(total_record_df, logger = logger)
        if model_status == "Infeasible":
            failed_lp = True

    if not failed_lp:
        mismap_hids.update(set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]]))
        mismap_hids.update(set([hid for hid in hap_max_sim_scores if hap_max_sim_scores[hid] > 10]))
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
        mismap_hids.update(set([hid for hid in hap_max_sim_scores if hap_max_sim_scores[hid] > 10]))
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        total_qnames = set([qname for qnames in hap_qname_info.values() for qname in qnames])
        correct_map_qnames = total_qnames - mismap_qnames
        return correct_map_qnames, mismap_qnames