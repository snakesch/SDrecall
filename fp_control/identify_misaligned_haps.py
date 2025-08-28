import numba
import re
import numpy as np
import pandas as pd
import pybedtools as pb

from collections import defaultdict
from numba import types, prange

from fp_control.bam_ncls import overlapping_reads_iterator
from fp_control.bilc import lp_solve_remained_haplotypes
from fp_control.numba_operators import numba_sum, extract_true_stretches
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





@numba.njit(types.Tuple((types.int32[:], types.int32[:], types.int32[:]))(types.int16[:], types.int16[:], types.int32), fastmath=True)
def numba_shared_variant_positions(vec1,
                                   vec2,
                                   overlap_start):
    """
    Return absolute genomic positions for shared SNVs and shared indels between two
    haplotype vectors aligned over the same overlap span.

    Inputs:
      - vec1, vec2: int16 haplotype vectors (same length)
        encodings: SNV=-4; deletion=-6 over length; insertion>1 at single index; match=1; NA<=-8
      - overlap_start: absolute genomic start (int32) of the vectors' first index

    Outputs:
      - shared_snv_pos_abs: int32 array of absolute positions where both have SNV (-4)
      - shared_indel_pos_abs: int32 array of absolute positions where both indicate indel
        (either deletion -6 or insertion > 1)
    """
    n = vec1.size
    # Preallocate with max possible and truncate
    snv_out = np.empty(n, dtype=np.int32)
    ins_out = np.empty(n, dtype=np.int32)
    del_out = np.empty(n, dtype=np.int32)
    snv_count = 0
    ins_count = 0
    del_count = 0
    for i in range(n):
        v1 = vec1[i]
        v2 = vec2[i]
        # Shared SNV
        if v1 == -4 and v2 == -4:
            snv_out[snv_count] = overlap_start + i
            snv_count += 1
        # Shared indel: deletion (-6) or insertion marker (>1)
        # Note insertion is marked at a single index with positive value > 1
        is_ins1 = (v1 > 1)
        is_ins2 = (v2 > 1)
        if is_ins1 and is_ins2 and v1 == v2:
            ins_out[ins_count] = overlap_start + i
            ins_count += 1

    # We need to find the identical block of deletion span to determine that they share an deletion event, not 1 single base deviated in deletion span
    vec1_del_spans = extract_true_stretches(vec1 == -6)
    vec2_del_spans = extract_true_stretches(vec2 == -6)
    for vec1_del_span in vec1_del_spans:
        for vec2_del_span in vec2_del_spans:
            if (vec1_del_span == vec2_del_span).all():
                for i in range(vec1_del_span[0], vec1_del_span[1] + 1):
                    del_out[del_count] = overlap_start + i
                    del_count += 1
    return snv_out[:snv_count], ins_out[:ins_count], del_out[:del_count]


@numba.njit(types.void(types.int32[:], types.int32, types.int32[:], types.int8[:], types.int32[:, :]), fastmath=True)
def update_tally_for_read(shared_pos_abs,
                          read_start,
                          ref_qseq_positions,
                          qseq_encoded,
                          tally):
    """
    In-place update of tally for one read.
    - shared_pos_abs: absolute genomic positions (int32[:]) to tally
    - read_start: read.reference_start (int32)
    - ref_qseq_positions: mapping ref_offset -> query_index (int32[:])
    - qseq_encoded: encoded query sequence (int8[:]); A=0,T=1,C=2,G=3,N=4
    - tally: (num_pos, 5) int32 matrix; rows align with shared_pos_abs order
    """
    num_pos = shared_pos_abs.size
    for i in range(num_pos):
        pos = shared_pos_abs[i]
        idx = pos - read_start
        if idx < 0 or idx >= ref_qseq_positions.size:
            continue
        qidx = ref_qseq_positions[idx]
        if qidx < 0:
            continue
        base = qseq_encoded[qidx]
        if base >= 0 and base < 4:  # exclude N (4)
            tally[i, base] += 1


@numba.njit(types.int8[:](types.int32[:], types.int32, types.int32[:], types.int8[:]), fastmath=True)
def map_positions_to_bases(shared_pos_abs,
                           read_start,
                           ref_qseq_positions,
                           qseq_encoded):
    """
    Map absolute genomic positions to the read's encoded base (0..3; N=4) via ref->query index map.
    Returns int8 array with -1 for positions not covered by the read.
    """
    n = shared_pos_abs.size
    out = np.full(n, np.int8(-1))
    for i in range(n):
        idx = shared_pos_abs[i] - read_start
        if idx < 0 or idx >= ref_qseq_positions.size:
            continue
        qidx = ref_qseq_positions[idx]
        if qidx >= 0:
            out[i] = qseq_encoded[qidx]
    return out


@numba.njit(types.int32(types.int8[:], types.int8[:]), fastmath=True)
def count_equal_alt(cons_codes,
                    hom_bases):
    """
    Count positions where both consensus ALT and homolog base are defined (>=0) and equal.
    """
    cnt = 0
    for i in range(cons_codes.size):
        if cons_codes[i] >= 0 and hom_bases[i] >= 0 and cons_codes[i] == hom_bases[i]:
            cnt += 1
    return cnt


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
    consensus_qual = np.full(end_pos - start_pos, 0.031622777, dtype=np.float32)

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
        mask = mask & (qual < 0.031622777)

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






@numba.njit(types.boolean(types.int16[:]), fastmath=True)
def judge_misalignment_by_extreme_vardensity(seq):
    five_vard = count_window_var_density(seq, padding_size = 42)
    six_vard = count_window_var_density(seq, padding_size = 65)
    read_vard = count_window_var_density(seq, padding_size = 74)
    # logger.debug(f"For the sequence {seq.tolist()}, the five_vard is {five_vard.tolist()}, the six_vard is {six_vard.tolist()}, the read_vard is {read_vard.tolist()}")

    if numba_sum(five_vard >= 5/85) > 0:
        select_bool = five_vard >= 5/85
        padding = 42
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

    if numba_sum(six_vard >= 6/131) > 0:
        select_bool = six_vard >= 6/131
        padding = 65
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
    
    if numba_sum(read_vard >= 12/148) > 0:
        return True
    
    return False



@numba.njit(types.float32[:](types.int32[:, :]), fastmath=True)
def calculate_coefficient(arr2d):
    # Span size is the weight for later mean value calculation
    rank_f = arr2d[:, 8].astype(np.float32)
    span = (arr2d[:, 1].astype(np.float32) - arr2d[:, 0].astype(np.float32)) / np.float32(100.0)
    depth_frac = (np.float32(1.0) - (arr2d[:, 4].astype(np.float32) / arr2d[:, 2].astype(np.float32)))

    res = (rank_f * span) * depth_frac
    return res



def sweep_region_inspection(input_bam, 
                            output_bed = None,
                            depth_cutoff = 5,
                            window_size = 120,
                            step_size = 30,
                            logger = logger):

    if output_bed is None:
        output_bed = prepare_tmp_file(suffix = ".bed").name

    cmd = f"samtools depth {input_bam} | \
            mawk -F '\\t' '$3 >= {depth_cutoff}{{printf \"%s\\t%s\\t%s\\n\", $1, $2-1, $2;}}' | \
            bedtools merge -i stdin | \
            bedtools makewindows -b stdin -w {window_size} -s {step_size} > {output_bed}"

    executeCmd(cmd, logger = logger)

    op_bed_obj = pb.BedTool(output_bed)
    return op_bed_obj


def calculate_coefficient_per_group(record_df, logger=logger):
    # Generate a rank for haplotypes
    record_df["rank"] = rank_unique_values(record_df["hap_max_sim_scores"].to_numpy(dtype=np.float32))
    # logger.info(f"Before calculating the coefficient for this region, the dataframe looks like :\n{record_df[:10].to_string(index=False)}\n")
    record_df["coefficient"] = calculate_coefficient(record_df.loc[:, ["start",
                                                                       "end",
                                                                       "total_depth",
                                                                       "hap_id",
                                                                       "hap_depth",
                                                                       "var_count",
                                                                       "indel_count",
                                                                       "psv_count",
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
        v_idx, qname = v
        label = vprop[v_idx]
        if label not in grouped_keys:
            grouped_keys[label] = [[], [], []]
        grouped_keys[label][0].append(v_idx)
        grouped_keys[label][1].append(qname)
        grouped_keys[label][2].append(rs)
    return grouped_keys



def record_haplotype_rank(haplotype_dict, mean_read_length = 150, hap_max_psv_pos = {}):
    starts = np.empty(len(haplotype_dict), dtype=np.int32)
    ends = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_depths = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_ids = np.empty(len(haplotype_dict), dtype=np.int32)
    var_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    indel_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    psv_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    total_depth_col = np.empty(len(haplotype_dict), dtype=np.int32)

    i = 0
    total_depth = 0
    for hid in haplotype_dict:
        seq, reads, span, _ = haplotype_dict[hid]
        starts[i] = span[0]
        ends[i] = span[1]
        # Calculate average depth over this span as:
        # total number of reference-aligned bases contributed by all reads within the span
        # divided by the span length (using the consensus sequence length for consistency).
        region_start = span[0]
        region_len = len(seq)
        region_end_excl = region_start + region_len
        cov_bases = 0
        for r in reads:
            s = max(r.reference_start, region_start)
            e = min(r.reference_end, region_end_excl)
            if e > s:
                cov_bases += (e - s)
        depth = int(cov_bases / region_len) if region_len > 0 else 0
        hap_depths[i] = depth
        hap_ids[i] = hid
        var_counts[i] = count_var(seq)
        indel_counts[i] = count_continuous_indel_blocks(seq)
        # Count the number of PSV positions for this haplotype that falls within this region
        psv_counts[i] = len(np.intersect1d(hap_max_psv_pos[hid], np.arange(region_start, region_end_excl)))
        total_depth += depth
        i += 1

    total_depth_col.fill(total_depth)
    return np.column_stack((starts, ends, total_depth_col, hap_ids, hap_depths, var_counts, indel_counts, psv_counts))



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
        logger.info(f"Shrink the region {region} to include more enclosing haplotypes: {inspect_results}")
        recover_results = [t for t in inspect_results if t[4] >= 0.8]
        recover_span = max(recover_results, key = lambda x: x[2][0])[2][0], min(recover_results, key = lambda x: x[2][1])[2][1]
        overlapping_span = (max(recover_span[0], start), min(recover_span[1], end))
        logger.info(f"Now the region is shrinked to {overlapping_span} and it contains the following haplotypes: {recover_results}")
        for t in recover_results:
            hap_id, qnames, span, sreads, _ = t
            region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), hap_id, qnames)
        if len(region_haplotype_info) <= 1:
            logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region even after the recover trial.\n")
            return None, None

    return region_haplotype_info, overlapping_span




def identify_misalignment_per_region(region,
                                     bam_ncls,
                                     qname_hap_info,
                                     qname_to_node,
                                     lowqual_qnames,
                                     clique_sep_component_idx,
                                     hap_max_psv_pos = {},
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
    qname_hap_info: a dictionary {vertex_idx: haplotype_id}
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
    debug_vertices = defaultdict(list)
    # a set of vertex indices (qname idx) corresponding to the overlapping reads
    vert_inds = set()
    for read in overlap_reads_iter:
        qname = read.query_name
        if qname in lowqual_qnames or qname not in qname_to_node:
            continue
        vert_idx = qname_to_node[qname]
        vertices[(vert_idx, qname)].append(read)
        debug_vertices[(vert_idx, qname)].append(get_read_id(read))
        vert_inds.add(vert_idx)
    logger.debug(f"Building the table for BILP, within region {region}, the vertices are:\n{'\n'.join([f'vertex_idx: {t[0]}, qname: {t[1]}, haplotype_id: {qname_hap_info[qname_to_node[t[1]]]}, reads: {read_ids}' for t, read_ids in debug_vertices.items()])}")

    # qname_hap_info is a dictionary {qname: haplotype_id}
    # Establish a dictionary of the form {haplotype_id: [vertex_idices, qnames, read_pair_lists]} based on the overlapping reads identified above
    hap_subgraphs = group_by_dict_optimized(qname_hap_info, vertices)

    if len(hap_subgraphs) == 0:
        logger.warning(f"No haplotype clusters are found for region {region}. These are the qnames found overlapping this region: {vertices}. Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict

    # For each hap_subgraph, we need to generate the consensus sequence for the reads
    # Each hap_subgraph is a connected component in the phasing graph and it represents a haplotype
    region_haplotype_info, overlapping_span = summarize_enclosing_haps(hap_subgraphs, qname_to_node, region, logger = logger)
    if region_haplotype_info is None:
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict
    logger.debug(f"The region {chrom}:{overlapping_span} contains the following haplotypes:\n{'\n'.join([f'{span}: {haplotype_idx} {qnames}' for span, (reads, vert_inds, haplotype_idx, qnames) in region_haplotype_info.items()])}")
    # Prepare region_str for logging
    region_str = f"{chrom}:{overlapping_span[0]}-{overlapping_span[1]}"

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

    logger.debug(f"The Final haplotype dict looks like this \n{'\n'.join([f'{haplotype_idx}: encoded_consensus_sequence: {con_seq.tolist()}, reads: {reads}, span: {span}, qnames: {qnames}' for haplotype_idx, (con_seq, reads, span, qnames) in final_clusters.items()])}")
    if len(final_clusters) <= 0:
        logger.warning(f"Only {len(final_clusters)} haplotype clusters are found for region {region_str}. Do not need to choose 2 haplotypes, Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict

    record_2d_arr = record_haplotype_rank( final_clusters, mean_read_length = mean_read_length, hap_max_psv_pos = hap_max_psv_pos )
    record_df = pd.DataFrame(record_2d_arr, columns = ["start", "end", "total_depth", "hap_id", "hap_depth", "var_count", "indel_count", "psv_count"])
    record_df["chrom"] = chrom

    logger.debug(f"Found {len(final_clusters)} haplotype clusters for region {region_str}. The dataframe recording the haplotypes in this region looks like :\n{record_df.to_string(index=False)}\n")
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
                           reads,
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
        if alt_varcount == 0:
            # Parse the homo_refseq_qname to get the original interval of the homologous sequence
            matched_groups = re.search(r"([a-zA-Z0-9]+):(\d+)-(\d+)[:RG0-9]*", homo_refseq_qname)
            if matched_groups is not None:
                homo_refseq_origin_chrom, homo_refseq_origin_start, homo_refseq_origin_end = matched_groups.groups()
                if homo_refseq_origin_chrom == chrom and int(homo_refseq_origin_start) <= overlap_span[0] and int(homo_refseq_origin_end) >= overlap_span[1]:
                    logger.warning(f"The homologous genomic sequence {homo_refseq_qname} overlapping interval {overlap_span[0]}-{overlap_span[1]} is aligned to its origin. So ignore this homologous sequence.")
                    continue

        # Compute positions with shared SNV and shared indel statuses in both vectors
        shared_snv_pos_abs, shared_ins_pos_abs, shared_del_pos_abs = numba_shared_variant_positions(interval_con_seq.astype(np.int16),
                                                                                  interval_genomic_hap.astype(np.int16),
                                                                                  np.int32(overlap_span[0]))
        # Default: empty verified SNV list; indels are directly verified as positions where both indicate indel
        verified_shared_snv_pos_abs = np.empty(0, dtype=np.int32)

        if shared_snv_pos_abs.size > 0:
            # Tally consensus ALT codes at shared positions from member reads
            num_pos = shared_snv_pos_abs.size
            tally = np.zeros((num_pos, 5), dtype=np.int32)  # bases 0..4

            # Use Numba helper to update tally per read
            for r in reads:
                ref_qseq_positions, qseq_ref_positions, qseq_encoded, _, read_ref_pos_dict = extract_read_qseqs(r, read_ref_pos_dict)
                update_tally_for_read(shared_snv_pos_abs.astype(np.int32),
                                      np.int32(r.reference_start),
                                      ref_qseq_positions,
                                      qseq_encoded,
                                      tally)

            # Majority base among A/T/C/G (codes 0..3)
            atcg = tally[:, :4]
            coverage = atcg.sum(axis=1)
            cons_codes = np.where(coverage > 0, atcg.argmax(axis=1).astype(np.int8), np.int8(-1))

            # ALT codes from the homologous read using Numba-mapped bases
            h_ref_positions, h_qseq_ref_positions, h_qseq_encoded, _, read_ref_pos_dict = extract_read_qseqs(homo_refseq, read_ref_pos_dict)
            h_base = map_positions_to_bases(shared_snv_pos_abs.astype(np.int32),
                                            np.int32(homo_refseq_start),
                                            h_ref_positions,
                                            h_qseq_encoded)

            # Verified shared PSV SNV positions: both defined and equal ALT
            eq_mask = (cons_codes >= 0) & (h_base >= 0) & (cons_codes == h_base)
            verified_shared_snv_pos_abs = shared_snv_pos_abs[eq_mask]

        # Shared PSV count for scoring: verified SNVs + shared indels
        shared_psv_ins_effect_size = 0
        for ins_pos in shared_ins_pos_abs:
            # Retrieve the encoded alignment status at the position of indels
            ins_encode_event = interval_genomic_hap[ins_pos - overlap_span[0]]
            if ins_encode_event > 0:
                logger.info(f"The homologous genomic sequence aligned at interval {chrom}:{overlap_span[0]}-{overlap_span[1]} shared an insertion at position {ins_pos}. The encoded event is {interval_genomic_hap[ins_pos - overlap_span[0]]}, the alignment status on consensus sequence is {interval_con_seq[ins_pos - overlap_span[0]]}")
                shared_psv_ins_effect_size += np.ceil(ins_encode_event/4)
            elif ins_encode_event < 0:
                logger.warning(f"The homologous genomic sequence aligned at interval {chrom}:{overlap_span[0]}-{overlap_span[1]} shared an insertion at position {ins_pos}. The encoded event is {interval_genomic_hap[ins_pos - overlap_span[0]]}, the alignment status on consensus sequence is {interval_con_seq[ins_pos - overlap_span[0]]}")

        shared_psv_del_effect_size = 0
        for del_pos in shared_del_pos_abs:
            # Retrieve the encoded alignment status at the position of indels
            del_encode_event = interval_genomic_hap[del_pos - overlap_span[0]]
            if del_encode_event < 0:
                logger.info(f"The homologous genomic sequence aligned at interval {chrom}:{overlap_span[0]}-{overlap_span[1]} shared a deletion at position {del_pos}. The encoded event is {interval_genomic_hap[del_pos - overlap_span[0]]}, the alignment status on consensus sequence is {interval_con_seq[del_pos - overlap_span[0]]}")
                shared_psv_del_effect_size += 1
            elif del_encode_event >= 0:
                logger.warning(f"The homologous genomic sequence aligned at interval {chrom}:{overlap_span[0]}-{overlap_span[1]} shared a deletion at position {del_pos}. The encoded event is {interval_genomic_hap[del_pos - overlap_span[0]]}, the alignment status on consensus sequence is {interval_con_seq[del_pos - overlap_span[0]]}")

        shared_psv = verified_shared_snv_pos_abs.size + shared_psv_ins_effect_size + shared_psv_del_effect_size
        verified_shared_indel_pos_abs = np.unique(np.concatenate([shared_ins_pos_abs, shared_del_pos_abs]))

        logger.info(f"For haplotype {hid}, within region {span}, comparing to the genomic sequence {homo_refseq_qname}. Variant count is {varcount}, alt-var count vs homologs is {alt_varcount}, shared PSV (ALT-consistent) count is {shared_psv}.")
        if homo_refseq_qname in varcounts_among_refseqs[hid]:
            # Each tuple records the stats across one continuous region of the haplotype hid
            varcounts_among_refseqs[hid][homo_refseq_qname].append((alt_varcount,
                                                                    shared_psv,
                                                                    verified_shared_snv_pos_abs,
                                                                    verified_shared_indel_pos_abs))
        else:
            varcounts_among_refseqs[hid][homo_refseq_qname] = [(alt_varcount,
                                                                shared_psv,
                                                                verified_shared_snv_pos_abs,
                                                                verified_shared_indel_pos_abs)]

    return varcounts_among_refseqs, read_ref_pos_dict



def cal_similarity_score(varcounts_among_refseqs, hid_var_count, logger = logger):
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
    hap_max_psvs = {}
    hap_max_psv_pos = {}
    # For each haplotype, find the maximum similarity delta.

    # Iterate over all the reference genome similarities for all the haplotypes
    for hid, gdict in varcounts_among_refseqs.items():
        max_psv = -100
        max_psv_c = 0
        max_psv_pos = np.empty(0, dtype=np.int32)
        for homo_refseq_qname, pairs in gdict.items():
            # pairs now contain tuples of (alt_varcount, shared_psv, verified_shared_snv_pos_abs, verified_shared_indel_pos_abs)
            total_shared_psv = np.sum([t[1] for t in pairs])
            alt_varcount = np.sum([t[0] for t in pairs])
            total_varcount = hid_var_count[hid]
            # Concat the verified shared SNV positions and indel positions
            verified_shared_snv_pos_abs = np.concatenate([t[2] for t in pairs])
            verified_shared_indel_pos_abs = np.concatenate([t[3] for t in pairs])
            # Unique the positions and sort them and conat snv and indel to a single array of abs genomic positions, sort in ascending order
            unique_psv_pos_abs = np.unique(np.concatenate([verified_shared_snv_pos_abs, verified_shared_indel_pos_abs]))
            psv_pos_abs = np.sort(unique_psv_pos_abs)
            # Count the number of unique positions
            shared_psv_ratio = total_shared_psv / total_varcount if total_varcount > 0 else min(1, total_shared_psv)
            mixed_psv_metric = 5 * shared_psv_ratio + total_shared_psv - (alt_varcount - total_shared_psv) + 0.1 * total_varcount
            logger.info(f"For haplotype {hid}, comparing to the reference sequence {homo_refseq_qname}, the similarity score is 4*{shared_psv_ratio} + {total_shared_psv} - ({alt_varcount} - {total_shared_psv}) + 0.1*{total_varcount} = {mixed_psv_metric}")

            if mixed_psv_metric > max_psv:
                max_psv_c = total_shared_psv
                max_psv = mixed_psv_metric
                max_psv_pos = psv_pos_abs

        hap_max_sim_scores[hid] = max_psv
        hap_max_psvs[hid] = max_psv_c
        hap_max_psv_pos[hid] = max_psv_pos  
        logger.info(f"For haplotype {hid}, the maximum similarity score is {max_psv}, the maximum PSV count is {max_psv_c}, the maximum PSV positions are {max_psv_pos}")

    return hap_max_sim_scores, hap_max_psvs, hap_max_psv_pos



def select_regions_with_min_haplotypes_from_hapbeds(hid_cov_beds,
                                                    output_bed=None,
                                                    min_haplotypes=3,
                                                    use_cli=True,
                                                    logger=logger):
    """
    Select intervals covered by >= min_haplotypes distinct haplotypes.

    hid_cov_beds: dict {hap_id: bed_path} for per-haplotype merged coverage.
    """
    if output_bed is None:
        output_bed = prepare_tmp_file(suffix=".hap_overlap.bed").name

    # Keep determinism
    bed_list = [hid_cov_beds[hid] for hid in sorted(hid_cov_beds)]
    if len(bed_list) < min_haplotypes:
        logger.warning(f"Only {len(bed_list)} haplotypes have coverage; need at least {min_haplotypes}.")
        return None

    if use_cli:
        # Build a small mapping file once
        map_file = prepare_tmp_file(suffix=".hid_beds.tsv").name
        with open(map_file, "w") as f:
            for hid in sorted(hid_cov_beds):
                f.write(f"{hid}\t{hid_cov_beds[hid]}\n")

        # Use BEDOPS to avoid argv limits and reproduce multiinter semantics:
        # 1) Tag each BED with haplotype ID in col4, then sort with sort-bed
        # 2) Partition the genome at all interval breakpoints
        # 3) Map partitioned intervals to tagged inputs, count distinct IDs, threshold by k
        # 4) Sort and merge
        tagged_bed = prepare_tmp_file(suffix=".tagged.bed").name
        parts_bed = prepare_tmp_file(suffix=".parts.bed").name

        cmd_tag_sort = (
            f"while IFS=$'\\t' read -r hid bed; do "
            f"  mawk -v id=\"$hid\" 'BEGIN{{FS=OFS=\"\\t\"}}{{print $1,$2,$3,id}}' \"$bed\"; "
            f"done < {map_file} | sort-bed - > "
            f"{tagged_bed}"
        )
        executeCmd(cmd_tag_sort, logger=logger)

        executeCmd(f"bedops --partition {tagged_bed} > {parts_bed}", logger=logger)

        cmd_map_filter = (
            "bedmap --ec --delim $'\\t' --echo --echo-map-id-uniq "
            f"{parts_bed} {tagged_bed} | "
            f"mawk -F '\\t' -v k={min_haplotypes} "
            "'{ ids=$NF; n=(ids==\"N/A\")?0:split(ids,a,\";\"); if (n>=k) print $1\"\\t\"$2\"\\t\"$3 }' | "
            f"sort-bed - > {output_bed}"
        )
        executeCmd(cmd_map_filter, logger=logger)

        bt = pb.BedTool(output_bed)
        return bt if bt.count() > 0 else None
    else:
        multi = pb.BedTool().multi_intersect(i=bed_list)
        # Keep only chrom,start,end; count in field 4
        filt = multi.filter(lambda iv: int(iv[3]) >= min_haplotypes).cut([0, 1, 2]).sort().merge()
        filt.saveas(output_bed)
        return pb.BedTool(output_bed) if filt.count() > 0 else None


def inspect_by_haplotypes(input_bam,
                          bam_ncls,
                          hap_qname_info,
                          qname_hap_info,
                          read_dict,
                          node_read_ids,
                          intrin_bam_ncls,
                          qname_to_node,
                          total_lowqual_qnames,
                          total_hap_vectors,
                          total_err_vectors,
                          total_genomic_haps,
                          read_ref_pos_dict,
                          compare_haplotype_meta_tab = "",
                          mean_read_length = 148,
                          logger = logger):
    record_dfs = []
    clique_sep_component_idx = 0
    hid_extreme_vard = defaultdict(bool) # Initialize a dictionary to store if the haplotype has extreme variant density
    hid_var_count = defaultdict(int) # Initialize a dictionary to store the variant count of the haplotype
    scatter_hid_dict = defaultdict(bool) # Initialize a dictionary to store if the haplotype is scattered
    logger.info("All the haplotype IDs are :\n{}\n".format(list(hap_qname_info.keys())))
    varcounts_among_refseqs = defaultdict(dict)
    total_qnames = set()
    hid_cov_beds = {}

    # Iterate over all the haplotypes
    for hid, qnames in hap_qname_info.items():
        # If only one read pair is in the iterating haplotype, it is a scattered haplotype, it should not be considered since the poor coverage
        qnames = [ qn for qn in qnames if qn not in total_lowqual_qnames ]
        total_qnames.update(set(qnames))
        if len(qnames) < 2:
            scatter_hid_dict[hid] = True

        logger.debug(f"The haplotype {hid} contains {len(qnames)} read pairs in total. And the qnames are: \n{qnames}\n")
        # hid_cov_bed = prepare_tmp_file(suffix = f".haplotype_{hid}.cov.bed", tmp_dir = tmp_dir).name

        # Extract all the qname indices and all the reads belonged to the iterating haplotype
        qname_indices = [qname_to_node[qname] for qname in qnames]
        read_ids = [rid for qname_idx in qname_indices for rid in node_read_ids[qname_idx]]
        reads = [read_dict[rid] for rid in read_ids]

        # Extract all the continuous regions covered by the iterating haplotype
        conregion_dict = extract_continuous_regions_dict(reads)
        # Write the continuous regions to a temporary bed file for this haplotype
        region_str = ""
        for span, reads in conregion_dict.items():
            region_str += f"{reads[0].reference_name}\t{span[0]}\t{span[1]}\n"

        hid_cov_bed = prepare_tmp_file(suffix = f".haplotype_{hid}.cov.bed").name
        pb.BedTool(region_str, from_string=True).sort().merge().saveas(hid_cov_bed)
        hid_cov_beds[hid] = hid_cov_bed

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
            logger.info(f"For haplotype {hid}, within region {chrom}:{span[0]}-{span[1]}, the consensus sequence is {consensus_sequence.tolist()}. The extreme variant density is {extreme_vard}. The varcount is {var_count}, till now the total varcount is {hid_var_count[hid]}")

            # Log the extreme variant density and the variant count of the iterating haplotype
            # logger.debug(f"For haplotype {hid} at continuous region {span}, the consensus sequence is {consensus_sequence.tolist()}. The extreme variant density is {extreme_vard}.")

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
                                                                                reads,
                                                                                total_genomic_haps,
                                                                                read_ref_pos_dict,
                                                                                varcounts_among_refseqs,
                                                                                logger = logger)
    # Log the extreme variant density haplotypes
    logger.info(f"extreme variant density haplotypes: {hid_extreme_vard}")
    # Log the scattered haplotypes
    logger.info(f"scatter haplotypes: {scatter_hid_dict}")
    # Log qnames that does not enter the ILP, and being directly categorized as mismap qnames
    scatter_qnames = set([qname for hid in scatter_hid_dict if scatter_hid_dict[hid] for qname in hap_qname_info[hid]])
    logger.info(f"{len(scatter_qnames)} qnames that does not enter the ILP for belonging to very small haplotypes, and being directly categorized as mismap qnames: {scatter_qnames}")

    # Calculate the maximum similarity score for all the haplotypes
    # Detailed explanation of the similarity score is in the docstring of the function cal_similarity_score
    hap_max_sim_scores, hap_max_psvs, hap_max_psv_pos = cal_similarity_score(varcounts_among_refseqs, hid_var_count, logger = logger)
    for hid in hap_qname_info.keys():
        if hid not in hap_max_sim_scores:
            logger.warning(f"The haplotype {hid} does not have any reference sequence similarity score. So we set the similarity score to 0.")
            hap_max_psv_pos[hid] = np.empty(0, dtype=np.int32)
    # Log the final reference sequence similarities for all the haplotypes
    logger.info(f"Final ref sequence similarities for all the haplotypes are {hap_max_sim_scores}")
    logger.info(f"Final ref sequence PSVs for all the haplotypes are {hap_max_psvs}")

    sweep_region_bed = input_bam.replace(".bam", ".sweep.bed")
    # Select regions overlapped by at least 3 haplotypes from cached per-haplotype coverage
    # sweep_regions = sweep_region_inspection(  input_bam,
    #                                           output_bed=sweep_region_bed,
    #                                           depth_cutoff=5,
    #                                           window_size=120,
    #                                           step_size=10,
    #                                           logger=logger )
    sweep_regions = select_regions_with_min_haplotypes_from_hapbeds(
        hid_cov_beds,
        output_bed=sweep_region_bed,
        min_haplotypes=1,
        use_cli=True,
        logger=logger
    )
    logger.info(f"Now the sweep regions are saved to {sweep_region_bed}.")

    '''
    The below iteration is to inspect haplotypes by each window decided by the sweep_region_inspection function.
    Within each window, we rank the enclosing haplotypes according to their similarity with the reference genomic sequence.
    '''
    if sweep_regions is None:
        mismap_hids = set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]])
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        mismap_qnames.update(scatter_qnames)
        correct_map_qnames = total_qnames - mismap_qnames
        logger.warning(f"No regions are found encompassed by more than 2 haplotypes. So we only filter out the haplotypes with extreme variant density. Filtered out {len(mismap_qnames)} read pairs and left {len(correct_map_qnames)} read pairs.")
        logger.info(f"Here is the qnames seemed mismapped: \n{mismap_qnames}\nAnd the mismap haplotype IDs: \n{mismap_hids}\n")
        return correct_map_qnames, mismap_qnames
    
    # Iterate over all the windows
    for interval in sweep_regions:
        region = (interval.chrom, interval.start, interval.end)
        # if interval.end - interval.start < 20:
        #     logger.debug(f"The region {region} is too small (spanning only {interval.end - interval.start} bp) for calculating the coefficients of haplotype items in the binary integer linear programming model.")
        #     continue
        # logger.debug(f"Inspecting the region {region} to compose the binary integer linear programming model.")
        record_df, total_hap_vectors, total_err_vectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx, read_ref_pos_dict = identify_misalignment_per_region(region,
                                                                                                                                                                            bam_ncls,
                                                                                                                                                                            qname_hap_info,
                                                                                                                                                                            qname_to_node,
                                                                                                                                                                            total_lowqual_qnames,
                                                                                                                                                                            clique_sep_component_idx,
                                                                                                                                                                            hap_max_psv_pos = hap_max_psv_pos,
                                                                                                                                                                            total_errvectors = total_err_vectors,
                                                                                                                                                                            total_genomic_haps = total_genomic_haps,
                                                                                                                                                                            read_ref_pos_dict = read_ref_pos_dict,
                                                                                                                                                                            mean_read_length = mean_read_length,
                                                                                                                                                                            logger = logger )
        if record_df is not None:
            record_dfs.append(record_df)
        else:
            logger.warning(f"No valid haplotypes covered by region {region} for downstream BILC")

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
        total_record_df["hap_max_psvs"] = total_record_df["hap_id"].map(hap_max_psvs).fillna(0)
        # total_record_df.loc[:, "coefficient"] = total_record_df["coefficient"] * 100 + total_record_df.loc[:, "var_count"]
        total_record_df.to_csv(compare_haplotype_meta_tab.replace(".tsv", ".raw.tsv"), sep = "\t", index = False)
        logger.info(f"Successfully saved the raw haplotype comparison meta table to {compare_haplotype_meta_tab.replace('.tsv', '.raw.tsv')}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        remove_hids = total_record_df.loc[(total_record_df["scatter_hap"]) | \
                                          (total_record_df["hap_max_sim_scores"] > 10) | \
                                          (total_record_df["extreme_vard"]), "hap_id"].unique()

        kept_scatter_hids = total_record_df.loc[(total_record_df["hap_var_count"] > 1) & (total_record_df["scatter_hap"]), "hap_id"].unique()
        remove_hids = set(remove_hids) - set(kept_scatter_hids)

        total_record_df = total_record_df.loc[np.logical_not(total_record_df["hap_id"].isin(remove_hids)), :]
        total_record_df = total_record_df.loc[total_record_df["total_depth"] > 5, :]
        if total_record_df.shape[0] == 0:
            logger.warning(f"No regions are left for BILC after filtering out the scattered haplotypes and the haplotypes with extreme variant density.")
            failed_lp = True
    else:
        logger.warning(f"No regions are left for BILC after filtering out the scattered haplotypes and the haplotypes with extreme variant density.")
        failed_lp = True

    if not failed_lp:
        by_region = total_record_df.groupby(["chrom", "start", "end"], group_keys=False)
        total_record_df = by_region.apply(calculate_coefficient_per_group, logger = logger).reset_index(drop=True)
        total_record_df.to_csv(compare_haplotype_meta_tab, sep = "\t", index = False)
        logger.info(f"Successfully saved the haplotype comparison meta table to {compare_haplotype_meta_tab}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        correct_map_hids, mismap_hids, model_status = lp_solve_remained_haplotypes(total_record_df, logger = logger)
        if model_status == "Infeasible":
            logger.warning(f"The BINARY INTEGER LINEAR PROGRAMMING failed to converge. So we need to filter out the haplotypes with extreme variant density and the haplotypes with high similarity with the reference genomic sequence.")
            failed_lp = True

    if not failed_lp:
        mismap_hids.update(set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]]))
        mismap_hids.update(set([hid for hid in hap_max_sim_scores if hap_max_sim_scores[hid] > 10]))
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] <= 1]))
        correct_map_hids = correct_map_hids - mismap_hids

        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        correct_map_qnames = total_qnames - mismap_qnames
        logger.info(f"Identify {len(mismap_hids)} mismapped haplotypes, {len(mismap_qnames)} mismapped qnames by solving the BINARY INTEGER LINEAR PROGRAMMING.\n")
        logger.info(f"Here is the qnames seemed mismapped: \n{mismap_qnames}\nAnd the mismap haplotype IDs: \n{mismap_hids}\n")
        logger.info(f"Here is the qnames seemed correctly mapped: \n{correct_map_qnames}\nAnd the correct map haplotype IDs: \n{correct_map_hids}\n")
        return correct_map_qnames, mismap_qnames

    if failed_lp:
        mismap_hids = set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]])
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] <= 1]))
        mismap_hids.update(set([hid for hid in hap_max_sim_scores if hap_max_sim_scores[hid] >= 8]))
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        total_qnames = set([qname for qnames in hap_qname_info.values() for qname in qnames])
        correct_map_qnames = total_qnames - mismap_qnames
        return correct_map_qnames, mismap_qnames


