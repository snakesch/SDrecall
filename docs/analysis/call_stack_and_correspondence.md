# Call Stack & Python-Rust Function Correspondence

## Call Stack with Data Flow

Legend: `->` = produced, `-->` = passed/consumed, `x` = erased/dropped after this step
Rust type abbreviations: `HM` = HashMap, `HS` = HashSet, `FxHM` = FxHashMap, `Arr2` = Array2

```
inspect_haplotypes_rust()                          [python_bindings.rs]   <- Python entry point
 |
 |  --> Inputs from Python (via PyO3):
 |      input_bam: &str                            (path to input BAM)
 |      intrinsic_bam: &str                        (path to intrinsic/homologous BAM)
 |      hap_qname_info: HM<i32, Vec<String>>       (hap_id -> [qnames])
 |      qname_hap_info: HM<i32, i32>               (vertex_idx -> hap_id)
 |      qname_to_node: HM<String, u32>              (qname -> vertex_idx)
 |      total_lowqual_qnames: HS<String>            (low-quality qnames to skip)
 |      compare_haplotype_meta_tab: &str            (output TSV path)
 |      mean_read_length: f64                       (default 148.0)
 |
 +-- build_lapper_from_bam(input_bam)              [bam_lappers.rs]
 |    |  -> bam_lapper: (HM<String, Lapper<u32,u32>>,    lapper_dict per chrom
 |    |                  FxHM<u32, Vec<Record>>,          read_dict: qname_idx -> [BAM records]
 |    |                  FxHM<u32, String>,               qname_dict: qname_idx -> qname
 |    |                  HM<String, u32>)                 qname_idx_dict: qname -> qname_idx
 |    +-- samtools_available()                       [bam_lappers.rs]
 |    +-- spawn_collate_pipe() [paired, Unix]        [bam_lappers.rs]
 |    +-- collate_bam_file() [paired, fallback]      [bam_lappers.rs]
 |    +-- process_qname_group()                      [bam_lappers.rs]
 |    |    +-- is_read_noisy()                       [bam_lappers.rs]
 |    |         +-- fast_median()                    [bam_lappers.rs]
 |    +-- get_reference_name()                       [bam_lappers.rs]
 |    +-- get_next_reference_name()                  [bam_lappers.rs]
 |
 +-- build_lapper_from_bam(intrinsic_bam)          [bam_lappers.rs]
 |    |  -> intrin_bam_lapper: (same 4-tuple type as above)
 |
 +-- inspect_haplotypes()                          [identify_misaligned_haps.rs]  (~475 lines, 10 phases)
 |    |
 |    |  --> Inputs:
 |    |      bam_lapper, intrin_bam_lapper          (from above)
 |    |      hap_qname_info, qname_hap_info, qname_to_node, total_lowqual_qnames
 |    |
 |    |  -> Accumulators initialized at top of function:
 |    |      hid_extreme_vard:    HM<i32, bool>        hap_id -> has extreme variant density?
 |    |      hid_var_count:       HM<i32, i32>         hap_id -> total variant count across regions
 |    |      scatter_hid_dict:    HM<i32, bool>        hap_id -> is scattered? (< 3 qnames)
 |    |      hid_max_local_density: HM<i32, f64>       hap_id -> max local density across regions
 |    |      varcounts_among_refseqs: HM<i32, HM<String, Vec<(i32,i32,i32,i32,Vec<i32>,Vec<i32>)>>>
 |    |                                                hap_id -> homo_refseq_id -> [(alt_varcount,
 |    |                                                   alt_snv_count, alt_indel_count,
 |    |                                                   shared_psv, verified_snv_pos, verified_indel_pos)]
 |    |      total_qnames:        HS<String>            all qnames seen (excl. lowqual)
 |    |      hid_cov_beds:        HM<i32, Vec<(String,i64,i64)>>  hap_id -> [(chrom,start,end)] merged
 |    |
 |    |  -> Caches (persist across all haplotypes and regions, grow monotonically):
 |    |      hap_cache:           HM<String, Array1<i16>>    read_id -> hap vector (cached)
 |    |      err_cache:           HM<String, Array1<f32>>    read_id -> err vector (cached)
 |    |      total_genomic_haps:  HM<String, Array1<i16>>    homo_refseq_id -> hap vector (cached)
 |    |      qseq_cache:          HM<String, ReadQseqData>   read_id -> ref/query pos maps (cached)
 |    |
 |    +---- [FOR EACH (hid, qnames) in hap_qname_info] --------
 |    |    |
 |    |    |  -> qnames (filtered): Vec<String>      (remove lowqual_qnames)
 |    |    |  -> total_qnames += qnames
 |    |    |  -> scatter_hid_dict[hid] = true if len(qnames) < 3
 |    |    |  -> reads: Vec<&Record>                 (from qname_to_records[qname], Option A lookup)
 |    |    |
 |    |    +-- extract_continuous_regions_dict()    [identify_misaligned_haps.rs]
 |    |    |    -> conregion_dict: BTreeMap<(i64,i64), Vec<&Record>>
 |    |    |
 |    |    |  -> hid_cov_beds[hid] = conregion_dict keys as bed intervals
 |    |    |
 |    |    +---- [FOR EACH (span, reads) in conregion_dict] --------
 |    |    |    |
 |    |    |    +-- record_hap_err_vectors_per_region()  [identify_misaligned_haps.rs]
 |    |    |    |    |  -> read_spans: Array2<i32>, hap_vectors: Array2<i16>, err_vectors: Array2<f32>
 |    |    |    |    +-- read_id()                       [pairwise_read_inspection.rs]
 |    |    |    |    +-- extract_hap_err_vectors()       [pairwise_read_inspection.rs]
 |    |    |    |    |    +-- extract_hap_vector()       [pairwise_read_inspection.rs]
 |    |    |    |    |    +-- extract_error_vector()     [pairwise_read_inspection.rs]
 |    |    |    |    |         +-- phred_to_prob()       [pairwise_read_inspection.rs]
 |    |    |    |    +-- extract_read_qseqs()            [pairwise_read_inspection.rs]
 |    |    |    |         +-- encode_base()              [pairwise_read_inspection.rs]
 |    |    |    |
 |    |    |    +-- assemble_consensus()                 [identify_misaligned_haps.rs]
 |    |    |    |    --> hap_vectors, err_vectors, read_spans
 |    |    |    |    -> consensus_sequence: Array1<i16>
 |    |    |    |
 |    |    |    +-- judge_misalignment_by_extreme_vardensity() [identify_misaligned_haps.rs]
 |    |    |    |    --> consensus_sequence
 |    |    |    |    -> extreme_vard: bool, region_max_density: f32
 |    |    |    |    +-- count_window_var_density()            [identify_misaligned_haps.rs]
 |    |    |    |    +-- count_snv()                           [pairwise_read_inspection.rs]
 |    |    |    |    +-- count_continuous_indel_blocks()       [pairwise_read_inspection.rs]
 |    |    |    |    +-- extract_true_stretches()              [identify_misaligned_haps.rs]
 |    |    |    |
 |    |    |    |  -> hid_max_local_density[hid] = max(old, region_max_density)
 |    |    |    |
 |    |    |    +-- count_var()                               [pairwise_read_inspection.rs]
 |    |    |    |    -> var_count: i32
 |    |    |    |
 |    |    |    |  -> hid_var_count[hid] += var_count
 |    |    |    |  -> hid_extreme_vard[hid] |= extreme_vard
 |    |    |    |
 |    |    |    +-- stat_refseq_similarity()                  [identify_misaligned_haps.rs]
 |    |    |         --> intrin_bam_lapper, chrom, span, hid, consensus_sequence, reads
 |    |    |         -> varcounts_among_refseqs[hid][homo_refseq_id] updated
 |    |    |         +-- query_overlapping_reads()            [bam_lappers.rs]
 |    |    |         +-- read_id()                            [pairwise_read_inspection.rs]
 |    |    |         +-- extract_hap_vector()                 [pairwise_read_inspection.rs]
 |    |    |         +-- extract_read_qseqs()                 [pairwise_read_inspection.rs]
 |    |    |         +-- ref_genome_similarity()              [identify_misaligned_haps.rs]
 |    |    |         +-- parse_origin_region()                [identify_misaligned_haps.rs]
 |    |    |         +-- numba_shared_variant_positions()     [identify_misaligned_haps.rs]
 |    |    |         +-- update_tally_for_read()              [identify_misaligned_haps.rs]
 |    |    |         +-- map_positions_to_bases()             [identify_misaligned_haps.rs]
 |    |    |         +-- verify_shared_snv_positions()        [identify_misaligned_haps.rs]
 |    |    |         +-- merge_unique_sorted()                [identify_misaligned_haps.rs]
 |    |    |
 |    |    |  x read_spans, hap_vectors, err_vectors, consensus_sequence dropped
 |    |    |
 |    |    +---- [END continuous region loop] --------
 |    |
 |    |  x conregion_dict, reads dropped
 |    |
 |    +---- [END haplotype loop] --------
 |    |
 |    |  Accumulators fully populated at this point.
 |    |
 |    +-- cal_similarity_score()                   [identify_misaligned_haps.rs]
 |    |    --> varcounts_among_refseqs, hid_var_count, hid_max_local_density
 |    |    -> hap_max_sim_scores: HM<i32, f64>, hap_max_psvs: HM<i32, i32>,
 |    |       hap_max_psv_pos: HM<i32, Vec<i32>>
 |    |
 |    |  x varcounts_among_refseqs dropped
 |    |
 |    +-- select_regions_with_min_haplotypes()     [identify_misaligned_haps.rs]
 |    |    --> hid_cov_beds
 |    |    -> sweep_regions: Option<Vec<(String, i64, i64)>>
 |    |
 |    |  -- Early return if sweep_regions is None --
 |    |  |  mismap_hids = {hid | hid_extreme_vard[hid]}
 |    |  |  return (total_qnames - mismap_qnames, mismap_qnames)
 |    |  --------
 |    |
 |    +---- [FOR EACH (chrom, start, end) in sweep_regions] --------
 |    |    |
 |    |    +-- identify_misalignment_per_region()  [identify_misaligned_haps.rs]
 |    |         --> region, bam_lapper, qname_hap_info, qname_to_node, etc.
 |    |         -> Option<IdentifyMisalignmentResult>:
 |    |             record_2d_arr: Arr2<i32> (N_haps x 8 cols)
 |    |             chrom: String
 |    |         +-- query_overlapping_reads()      [bam_lappers.rs]
 |    |         +-- group_by_dict_optimized()      [identify_misaligned_haps.rs]
 |    |         +-- summarize_enclosing_haps()     [identify_misaligned_haps.rs]
 |    |         +-- record_hap_err_vectors_per_region() [identify_misaligned_haps.rs]
 |    |         +-- assemble_consensus()           [identify_misaligned_haps.rs]
 |    |         +-- record_haplotype_rank()        [identify_misaligned_haps.rs]
 |    |
 |    +---- [END sweep region loop] --------
 |    |
 |    |  x caches dropped
 |    |
 |    |  === BILC PRE-PROCESSING (Python lines 1331-1408) ===
 |    |
 |    |  (a) Flatten record_results -> total_records: Vec<BilcRecord>
 |    |  (b) Enrich from accumulators: extreme_vard, scatter_hap, hap_var_count,
 |    |      hap_max_sim_scores (rounded 1dp), hap_max_psvs
 |    |  (c) Filter remove_hids (scatter | sim>10 | psvs>=12 | extreme_vard)
 |    |      EXCEPT kept_scatter_hids (var_count>=1 & scatter & psvs<12 & sim<=10 & !extreme)
 |    |      Then filter total_depth <= 5
 |    |  (d) Per-region: rank_unique_values() -> varc_rank, then calculate_coefficient()
 |    |      = psv_metric * sqrt(span/100) * sqrt(1 - hap_depth/total_depth)
 |    |  (e) Span-weighted average coefficient per hap_id
 |    |  (f) coefficient += hap_max_sim_scores
 |    |  (g) Rank coefficient per region -> varc_rank for ILP
 |    |  (h) Deduplicate total_records
 |    |
 |    |  === BILC SOLVER (bilc.py:62-148) [bilc_solver.rs] ===
 |    |
 |    +-- lp_solve_remained_haplotypes()           [bilc_solver.rs]
 |    |    Variables: x_i in {0,1} per unique hap_id
 |    |    Objective: minimize sum(-coefficient_i * x_i)
 |    |    Constraints per region: 0 <= sum(x_j) <= upper_bound
 |    |      upper_bound: n_haps>4 AND sum(var_count)>=1 -> n_haps-2
 |    |                   n_haps<=1 -> n_haps
 |    |                   else -> n_haps - count(varc_rank<=1)
 |    |    -> select_hap_ids, drop_hap_ids, status, objective
 |    |
 |    |  === POST-ILP AUGMENTATION (Python lines 1415-1431) ===
 |    |
 |    |  mismap_hids = drop_hap_ids
 |    |    U {hid | extreme_vard}
 |    |    U {hid | sim_score > 15}
 |    |    U {hid | scatter & var_count <= 3}
 |    |    U remove_hids
 |    |
 |    |  -- Fallback if failed_lp --
 |    |  |  mismap_hids = {extreme_vard} U {scatter} U {sim>15}
 |    |  --------
 |    |
 |    +-- return (correct_qnames, mismap_qnames)
 |
 +-- [Return to Python]
```

---

## Python-Rust Function Correspondence (46 functions)

| # | Rust Function | Rust File | Python Function | Python File | Notes |
|---|--------------|-----------|-----------------|-------------|-------|
| **BAM I/O & Interval Queries** | | | | | |
| 1 | `build_lapper_from_bam` | bam_lappers.rs | `migrate_bam_to_ncls` | bam_ncls.py:288 | |
| 1a | `samtools_available` | bam_lappers.rs | *(implicit `shutil.which`)* | bam_ncls.py:199 | |
| 1b | `spawn_collate_pipe` | bam_lappers.rs | `_collate_bam_file` | bam_ncls.py:196 | pipe via `/dev/fd/N` |
| 1c | `collate_bam_file` | bam_lappers.rs | `_collate_bam_file` | bam_ncls.py:196 | temp file fallback |
| 2 | `process_qname_group` | bam_lappers.rs | `_process_qname_group` | bam_ncls.py:219 | |
| 3 | `is_read_noisy` | bam_lappers.rs | `is_read_noisy` | bam_ncls.py:115 | |
| 4 | `fast_median` | bam_lappers.rs | *(inline)* | bam_ncls.py:133 | |
| 5 | `query_overlapping_reads` | bam_lappers.rs | `overlapping_reads_iterator` | bam_ncls.py:13 | |
| 6 | `query_overlapping_qname_indices` | bam_lappers.rs | `overlap_qname_idx_iterator` | bam_ncls.py:56 | |
| 7 | `query_overlapping_read_pairs` | bam_lappers.rs | *(implicit)* | -- | |
| **Pairwise Read Inspection** | | | | | |
| 8 | `read_id` | pairwise_read_inspection.rs | `get_read_id` | pairwise_read_inspection.py:91 | |
| 9 | `extract_hap_vector` | pairwise_read_inspection.rs | `get_hapvector_from_cigar` | pairwise_read_inspection.py:165 | |
| 10 | `extract_error_vector` | pairwise_read_inspection.rs | `get_errorvector_from_cigar` | pairwise_read_inspection.py:269 | |
| 11 | `extract_hap_err_vectors` | pairwise_read_inspection.rs | *(combined #9+#10)* | -- | |
| 12 | `phred_to_prob` | pairwise_read_inspection.rs | *(inline)* | pairwise_read_inspection.py:283 | |
| 13 | `record_hap_err_vectors_per_region` | identify_misaligned_haps.rs | `record_hap_err_vectors_per_region` | identify_misaligned_haps.py:742 | |
| 14 | `assemble_consensus` | identify_misaligned_haps.rs | `assemble_consensus` | identify_misaligned_haps.py:289 | |
| 15 | `encode_base` | pairwise_read_inspection.rs | *(inline base_dict)* | pairwise_read_inspection.py:446 | |
| 16 | `extract_read_qseqs` | pairwise_read_inspection.rs | `extract_read_qseqs` + `prepare_ref_query_idx_map` | pairwise_read_inspection.py:446,99 | Fused into single CIGAR walk |
| 17 | `extract_read_qseqs_cached` | pairwise_read_inspection.rs | *(implicit dict caching)* | -- | |
| **Similarity & Variant Analysis** | | | | | |
| 18 | `count_snv` | pairwise_read_inspection.rs | `count_snv` | pairwise_read_inspection.py:45 | |
| 19 | `count_continuous_indel_blocks` | pairwise_read_inspection.rs | `count_continuous_indel_blocks` | pairwise_read_inspection.py:61 | |
| 20 | `count_var` | pairwise_read_inspection.rs | `count_var` | pairwise_read_inspection.py:71 | |
| 21 | `count_continuous_blocks` | pairwise_read_inspection.rs | `count_continuous_blocks` | pairwise_read_inspection.py:25 | |
| 22 | `count_window_var_density` | identify_misaligned_haps.rs | `count_window_var_density` | identify_misaligned_haps.py:347 | |
| 23 | `extract_true_stretches` | identify_misaligned_haps.rs | *(inline)* | identify_misaligned_haps.py:395 | |
| 24 | `judge_misalignment_by_extreme_vardensity` | identify_misaligned_haps.rs | `judge_misalignment_by_extreme_vardensity` | identify_misaligned_haps.py:370 | |
| 25 | `edit_distance` | pairwise_read_inspection.rs | *(not in Python)* | -- | Added for future use |
| 26 | `merge_unique_sorted` | identify_misaligned_haps.rs | `merge_unique_sorted` | identify_misaligned_haps.py:244 | |
| 27 | `ref_genome_similarity` | identify_misaligned_haps.rs | `ref_genome_similarity` | identify_misaligned_haps.py:32 | |
| 28 | `numba_shared_variant_positions` | identify_misaligned_haps.rs | `numba_shared_variant_positions` | identify_misaligned_haps.py:58 | |
| 29 | `map_positions_to_bases` | identify_misaligned_haps.rs | `map_positions_to_bases` | identify_misaligned_haps.py:174 | |
| 30 | `update_tally_for_read` | identify_misaligned_haps.rs | `update_tally_for_read` | identify_misaligned_haps.py:111 | |
| 31 | `verify_shared_snv_positions` | identify_misaligned_haps.rs | `verify_shared_snv_positions` | identify_misaligned_haps.py:208 | |
| 32 | `parse_origin_region` | identify_misaligned_haps.rs | `re.search(origin pattern)` | identify_misaligned_haps.py:947 | Regex-free manual parser |
| 33 | `stat_refseq_similarity` | identify_misaligned_haps.rs | `stat_refseq_similarity` | identify_misaligned_haps.py:888 | |
| **Inspection Loop & Scoring** | | | | | |
| 34 | `extract_continuous_regions_dict` | identify_misaligned_haps.rs | `extract_continuous_regions_dict` | identify_misaligned_haps.py:477 | |
| 35 | `cal_similarity_score` | identify_misaligned_haps.rs | `cal_similarity_score` | identify_misaligned_haps.py:1022 | |
| 36 | `inspect_haplotypes` | identify_misaligned_haps.rs | `inspect_by_haplotypes` | identify_misaligned_haps.py:1153 | ~475 lines, 10 phases |
| 37 | `identify_misalignment_per_region` | identify_misaligned_haps.rs | `identify_misalignment_per_region` | identify_misaligned_haps.py:637 | |
| 38 | `group_by_dict_optimized` | identify_misaligned_haps.rs | `group_by_dict_optimized` | identify_misaligned_haps.py:516 | |
| 39 | `summarize_enclosing_haps` | identify_misaligned_haps.rs | `summarize_enclosing_haps` | identify_misaligned_haps.py:573 | |
| 40 | `record_haplotype_rank` | identify_misaligned_haps.rs | `record_haplotype_rank` | identify_misaligned_haps.py:530 | |
| 41 | `select_regions_with_min_haplotypes` | identify_misaligned_haps.rs | `select_regions_with_min_haplotypes_from_hapbeds` | identify_misaligned_haps.py:1089 | In-memory sweep-line |
| **ILP Post-processing** | | | | | |
| 42 | `rank_unique_values` | identify_misaligned_haps.rs | `rank_unique_values` | identify_misaligned_haps.py:270 | Cross-validated vs numba |
| 43 | `calculate_coefficient` | identify_misaligned_haps.rs | `calculate_coefficient` | identify_misaligned_haps.py:425 | Cross-validated vs numba |
| **BILC ILP Solver** | | | | | |
| 44 | `lp_solve_remained_haplotypes` | bilc_solver.rs | `lp_solve_remained_haplotypes` | bilc.py:62 | Cross-validated (1,568 files) |
| 45 | `lp_solve_remained_haplotypes_with_obj` | bilc_solver.rs | *(extended)* | -- | Returns objective value |
| **Python Bindings** | | | | | |
| 46 | `inspect_haplotypes_rust` | python_bindings.rs | *(new entry point)* | -- | |
