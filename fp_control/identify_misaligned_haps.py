import numba
import numpy as np
import pandas as pd
from collections import defaultdict
from numba import types


from bilc import lp_solve_remained_haplotypes
from numba_operators import numba_sum
from fp_control.pairwise_read_inspection import get_hapvector_from_cigar, \
                                   get_errorvector_from_cigar, \
                                   get_read_id, \
                                   count_var, \
                                   count_continuous_blocks


@numba.njit(types.Tuple((types.int32, types.int32))(types.int32[:], types.int32[:]), fastmath=True)
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