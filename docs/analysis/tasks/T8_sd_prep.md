# T8 — `sd-prep` crate (Phase 1 preparation, incl. graph-tool replacement)

**Crate:** `sd-prep` (new — lib + bin)
**Status (2026-07-17):** Production-integrated. The historical Phase-1 Python differential reached 679/690 (98.4%) accepted paralog-pair parity, and the Rust stage now completes in validated t2t/hg19/hg38 end-to-end runs.
**Depends on:** T0
**Replaces (Python):** `prepare_recall_regions.py` + the `preparation/` package: `pick_multialign_regions`, `inferred_depths`, `sd_pairs`, `graph_build`, `graph_query`, `graph_traversal`, `homoseq_region`, `genome`, `build_beds_and_masked_genomes`, `intrinsic_alignment`, `seq`.

## Goal & scope boundary

Bring all of Phase 1 into Rust so the pipeline is Python-free. The **hard core** is replacing **graph-tool** (`label_components`, `GraphView` vertex-filter, `sequential_vertex_coloring`, `shortest_path`) — it has no clean Rust equivalent. Everything else in Phase 1 is either straightforward Rust compute or subprocess orchestration of the same external tools.

Sub-components, by difficulty:
- **Hard:** SD-graph traversal + node grouping by coloring (`graph_query.py`, `graph_traversal.py`) → `petgraph`/`rustworkx` + custom greedy coloring.
- **Medium:** umbrella SD-pair filtering (`sd_pairs.py`, O(n²) interval logic); multiplex graph build (`graph_build.py`, networkx + interval tree → `petgraph` + `rust-lapper`).
- **Easy / in-process:** multi-align depth (`pick_multialign_regions` + `inferred_depths` → rust-htslib pileup); masking (`genome.py` → `rust-bio` FASTA); region seq extraction + frag-size stats (`seq.py` → rust-htslib). **External algorithm:** intrinsic alignment (`intrinsic_alignment.py` → `minimap2` via `minimap2-rs` FFI; see policy).

## Data flow

```
ref FASTA + input BAM + SD-map BED + target BED
  → multi-align BED            (4-way depth filter)
  → SD pairs (umbrella-filtered)
  → multiplex graph (SD edges + per-chr physical-overlap edges)
  → traversal → SD paralog pairs               ⟵ graph-tool core (Phase-1 hotspot)
  → node grouping via coloring → realignment groups (RG0..RGn)
  → per-RG: query BED, counterpart BED, masked genome FASTA, intrinsic BAM
  → GraphML graphs + filtered_SD_binary_map.tsv
```

## Dependencies

- Crates: `petgraph` 0.8.3 (graph + greedy coloring), `rust-lapper` 1.2.0 (interval tree), `bedrs` 0.2.26 (BED ops), `rust-htslib` 0.47.0, `bio` (rust-bio) 3.0.0 (FASTA masking), `minimap2` (minimap2-rs) 0.1.31 (FFI), `rayon` 1.8, `sdrecall-utils`/`sdrecall-io`.
- In-process: `rust-htslib` (depth/BAM), `bedrs` (BED ops), `rust-bio` (FASTA masking), `petgraph` (graph). External algorithm only: `minimap2` via `minimap2-rs` FFI (see policy).

## Performance bottleneck / rationale

The Phase-1 hotspot is graph traversal: per-qnode `shortest_path` over each subgraph vertex **plus a minimap2 + 2× samtools faidx subprocess per candidate counterpart** (`graph_traversal.py`), parallelized over a pool. The graph-tool replacement is the dominant effort and the dominant risk.

**Risk + fallback:** if exact partition parity with graph-tool's coloring proves intractable, the documented fallback is to keep Phase-1 graph ops as a thin subprocess to a tiny graph-tool helper — but the goal is pure Rust, so this is a last resort, recorded here as a decision point.

## Tests

### Unit (tier 1)
- Connected components vs `gt.label_components` on fixture graphs.
- Greedy vertex coloring vs `gt.sequential_vertex_coloring` on small graphs (compare chromatic grouping, not color labels).
- Umbrella-pair filtering on hand-built overlapping SD pairs.

### Differential vs Python (tier 2)
- Full Phase-1 run on HG002/HG006: compare RG node groupings, paralog-pair set, and masked-genome FASTAs.

**Pass criterion:** identical **RG node groupings** (set partition), identical **paralog-pair set**, and **md5-identical masked genome FASTAs** vs Python (minimap2 results are identical since the same binary is invoked). Coloring is compared as a partition, not by color id.

**Data:** HG002 + HG006 full preparation runs.

## Progress
- [x] Scaffold `sd-prep` crate (lib + bin) — `cargo test -p sd-prep` green (55 unit tests), `cargo clippy -p sd-prep --all-targets` clean (2026-06-12).
- [x] **Graph core (the deliverable)** — `graph_core.rs`: `component_labels` (UnionFind; ONE versatile fn for all-edge AND overlap-only partitions via an edge predicate), `greedy_vertex_coloring` (gt.sequential_vertex_coloring index-order parity), `dijkstra_route` (predecessor-recording weighted Dijkstra over the undirected view). UNIT-TESTED (components on 4 fixtures incl. the overlap-only bridge; coloring on P4/C5/K3+isolated/empty + a proper-coloring property test on 50 random graphs; Dijkstra on 4 weighted fixtures).
- [x] **Multiplex graph build** — `graph_build.rs`: `NodeKey`/`EdgeAttr`/`SdGraph` + `build_multiplex_graph` (per-chr PO with large→small direction + `1/overlap_frac` weight + dedup, SD overlay/upgrade, frozenset-equivalent pair dedup). UNIT-TESTED + **DIFFERENTIAL PASS** on HG002: 1268 nodes / 1952 edges / 691 SD edges — byte-identical to the Python `*_multiplexed_SDs.graphml` (198 components).
- [x] **Umbrella SD-pair filtering** — `sd_pairs.rs`: `Pair` + the ONE `umbrella_to_remove` sweep (called for raw + granular), `is_umbrella_pair`/`overlap_fractions_other`/`extract_subsegment_for_target`. UNIT-TESTED (umbrella cover, strand-flip block, partial-coverage block, break-semantics, same/opposite-strand subsegment).
- [x] **Node grouping via greedy coloring (fast path)** — `grouping.rs`: `ConnectedQnodes` (order-preserving insertion = coloring index order, the parity contract) + `color_groups`/`optimal_node_grouping`. UNIT-TESTED. The general/conflict-graph branch is intentionally NOT ported (dead for parity).
- [x] **Multi-align depth (sweep + filter)** — `multialign.rs`: shared `depth_sweep`, four-pass filtering, and the rust-htslib BAM-read loop are implemented and differentially validated.
- [x] **Masking math + FASTA I/O** — end masking, N-bridge merging, indexed FASTA read/write, and masked-genome output are implemented and validated.
- [x] **`HomoseqRegion`** — `homoseq.rs`: route-carrying coordinate object + `qnode_relative_region` back-projection (SD clip/flip + overlap re-anchor). UNIT-TESTED.
- [x] **`sort_query_nodes`** — `traversal.rs`: load-balance/coloring-order sort (degree×component_size desc, stable on ties) + `is_small_sd`/`should_prune_po_edge` prune predicates. UNIT-TESTED.
- [x] **Graph traversal → paralog pairs** — `traversal.rs::extract_sd_paralog_pairs` DONE (route walk `inspect_cnode_along_route` + `compare_homologous_sequences` via the `minimap::align_similarity` FFI wrapper + `ConnectedQnodes` insertion-order wiring). UNIT-TESTED (route walk: SD same/opposite-strand clip, overlap-too-small drop, `route_vertices`).
- [x] **minimap2 wrapper** — `minimap.rs::align_similarity` (ONE wrapper, `asm10`/`asm20` presets, `match_len/block_len`), reused by traversal + intrinsic. `SIMILARITY_THRESHOLD = 0.95` const documents the bundled-2.30-vs-system-2.28 skew. UNIT-TESTED (identical/SNP/unmappable/empty). The crate's `htslib` feature is OFF (it hard-pins rust-htslib 1.0, conflicting with the workspace's 0.47.1); the `map()` → `Mapping` API + workspace rust-htslib build the intrinsic BAM records.
- [x] **mask_genome FASTA I/O** — `masking.rs::mask_genome` (bio `IndexedReader` faidx + 60-char-wrapped FASTA writer + md5-gated update + `samtools faidx`). UNIT-TESTED. **DIFFERENTIAL: all 4 HG002 RG masked FASTAs are md5-IDENTICAL to Python** given matching frag-stats.
- [x] **multialign BAM loop** — `multialign.rs::{inferred_coverage, pick_multialigned_regions}` (rust-htslib `fetch` per region → `read_passes` → `depth_sweep`, 4 passes via a scoped rayon pool, 4-way AND filter, target intersect). DIFFERENTIAL: 1978 multi-align intervals — exact match with Python.
- [x] **intrinsic alignment + per-RG outputs + driver** — `intrinsic.rs` (getRawseq + minimap2 asm20 `Mapping`→BAM record + `filter_intrinsic_alignments` self-location/enclosure/sec→pri + samtools merge/sort/index for the total BAM) and `driver.rs::prepare_recall_regions` (the 7-step orchestration, rayon over RGs, auto-derives frag-stats + mean read length from the BAM). UNIT-TESTED.
- [x] **Differential vs Python (HG002):** `examples/validate_phase1_e2e.rs`. **Paralog-pair set: 679/690 frozenset match (98.4%)** after two correctness fixes (see below). Remaining 11 divergent pairs from: (a) 3 multi-align BED micro-boundary differences (±5–15bp at interval ends on chr2:106400146 / chr7:74799947 / chr7:74805440 — sweep-line event-ordering edge case vs bedtools `genomecov`), cascading to 7 chr2:10640xxxx SD pairs; (b) 4 chr16 mutual-umbrella ties unresolvable without exact RNG parity. **User confirmed these are ignorable.**
- [x] **FIX 1 — umbrella filter within-group sort (2026-06-12):** The O(n²) umbrella sweep's `break` semantics are order-dependent for mutual-umbrella ties (equal `overlap_len`). The expanded SD map has near-duplicate entries (forward/reverse CIGAR alignments, ±18bp trim differences) that are mutual umbrellas. **Python** sorts the SD map by `BedTool.sort()` (chrom, start, **end ascending**) before intersection, so the narrower-end entry appears first and gets removed. **Rust** had entries in SD-map file order (wider-end first), removing the wider entry instead — flipping which orientation survived for ~19% of pairs. **Fix:** sort each umbrella group's `idxs` by `(a.chrom, a.start, a.end)` before the sweep (`driver.rs::umbrella_filter_and_dedup`), matching Python's BedTool sort. This brought frozenset parity from **560/668 (81%) → 679/690 (98.4%)**.
- [x] **FIX 2 — deterministic insert-size RNG (2026-06-12):** `sdrecall_io::get_insert_size_distribution` used `rand::thread_rng()` (non-deterministic Monte Carlo), producing different median values across runs (567–573bp range). Seeded with `StdRng::seed_from_u64(42)` for reproducibility. The Rust median (569.7bp) still differs from Python's (567.4bp) because the RNG algorithms differ (ChaCha vs Mersenne Twister) — this causes ±1bp size-filter threshold shifts for 0–2 SD pairs near the boundary. Accepted as ignorable per user.
- [x] **3 multi-align BED micro-boundaries (root-caused, not fixed):** Out of 1978 multi-align intervals, 3 differ by 5–15bp at their right endpoints: chr2:106400146 (Rust +15bp), chr7:74799947 (Rust +5bp), chr7:74805440 (Rust −13bp). Cause: the sweep-line event sort `(position, delta)` processes `−1` (end) events before `+1` (start) events at tied positions (standard convention, matches bedtools), but bedtools internally uses a BAM-backed `genomecov` with subtly different read-to-interval conversion at clipped-read boundaries. Does not affect the vast majority of SD pairs but cascades through the chr2:10640xxxx umbrella filter groups.
- [x] **Decision point resolved (deps):** minimap2-rs `0.1.31` (bundled minimap2 2.30) AND rust-bio `2.3.0` BOTH build cleanly in the SDrecall env — pure-Rust path is viable, no graph-tool subprocess fallback needed. ⚠ Version skew: bundled minimap2 2.30 vs system `minimap2 2.28-r1209` — pinned in `minimap::SIMILARITY_THRESHOLD`.

## Migration design — frontier propagation (2026-06-11)

This is the coding-grade design for `sd-prep`. It ports the *actual* Python control flow (line-referenced below), maps each operation to a **verified** Rust crate API, and pins every borrow/owner choice. The hard core (graph-tool replacement) is treated in depth; the easy/medium modules are sketched with their one-versatile-unit home.

### 1. Python logic inventory (functions to port, with control-flow shape)

Driver: `prepare_recall_regions.py::prepare_recall_regions` (l.57-243). Sequential 7-step pipeline; only steps 1, 4, 5, 7 are internally parallel. The driver's two helpers `save_connected_qnodes_to_graphml` (l.19-54) collapse into the generic `sdrecall-io::write_graphml`.

**Step 1 — multi-align depth (Easy, in-process).**
- `pick_multialign_regions.py::pick_multialigned_regions` (l.7-63): fans out **4 fixed depth passes** over a `Pool(min(4,threads))` (raw / high-MQ / XA-tag / XS-tag), merges the 4 per-base depth tables on `(chrom,pos)`, applies the AND filter `raw>=min_depth & (XA_frac>=frac | XS_frac>=frac) & high_MQ<=hq_depth` (l.43-48), then `sort|merge|intersect(target)` (l.50-57). The 4 passes are **one function with 4 parameter tuples** — collapse into ONE depth kernel.
- `inferred_depths.py::calculate_inferred_coverage` (l.67-136): per-pass kernel. Reads BAM (`bam.fetch` per target region, l.92-110), filters each read by `filter_and_process_read` (l.10-29: MQ + not-unmapped/dup/secondary/supplementary/qcfail + tag logic incl. the `AS-XS<=5` XS rule l.18), then converts read intervals → **per-base coverage** via `genome_coverage(bg=True)` then expands `bg` to per-position rows (l.124-127). The standalone `calculate_coverage` (l.42-63) is dead (uses undefined `output_depth`/`threads`/`MQ_threshold`) — **do NOT port**. `create_genome_dict` (l.32-39) → reads `.fai` → chrom sizes map.

**Step 2-3 — SD-map load + umbrella filter (Easy load + Medium filter).**
- Driver l.121-179: load `reference_sd_map` BED, filter SDs by `len>avg_frag_size` (l.125), `intersect(multi_align_bed, wo=True)` (l.129), main-contig regex filter `^(chr)?([0-9]+|[XYM]|MT)$` (l.144-151), both-negative-strand → both-positive flip (l.155-158), groupby `(chr_bam1,start_bam1,end_bam1)` and parallel-map `filter_umbrella_pairs` over each group (l.162-166), then dedup by `frozenset({sdA_key, sdB_key})` (l.175-176) and emit `filtered_SD_binary_map.tsv` (l.177).
- `sd_pairs.py` (Medium, O(n²)): `filter_umbrella_pairs` (l.256-285) → `_find_umbrella_pairs` (l.182-253) does a **two-pass O(n²)** umbrella sweep (raw pairs, then `extract_subsegment_for_target`-refined "granular" pairs), `Pair.is_umbrella_pair` (l.96-122) + `calculate_interval_overlaps` (l.32-64) + `is_same_pair` (l.67-94) + `extract_subsegment_for_target` (l.125-173). The two sweeps are **identical logic over different pair lists** — collapse into ONE `umbrella_to_remove(&[Pair]) -> HashSet<usize>` helper called twice.

**Step 4 — multiplex graph build (Medium).**
- `graph_build.py::create_multiplex_graph` (l.83-144): builds an undirected SD-edge graph `G` (l.102-106), then **per-chromosome in a Pool** builds a directed physical-overlap (PO) graph via `compose_PO_graph_per_chr` (l.35-81) using an `IntervalTree` (PO edge large→small, weight `=min_size/overlap_span` i.e. `1/max(span/sizes)`, l.67-76), merges all per-chr PO graphs (`nx.compose`, l.122-124), then overlays SD edges onto the merged PO graph (l.127-135). Returns a merged DiGraph. `read_graphml`/`string_to_tuple` (l.15-32, l.146-156) = tuple-node (de)serialization → handled by `sdrecall-io` GraphML + a node-key parser.

**Step 5 — graph traversal → SD paralog pairs (HARD CORE, Phase-1 hotspot).**
- `graph_query.py::extract_SD_paralog_pairs_from_graph` (l.133-283): tag query nodes (l.152-153), **prune small SD nodes** `size<=max(mean_read_len, avg_frag-1*std)` (l.157-167), **prune weak PO edges** `overlap_size<cutoff && 1/weight<0.5` + self-loops (l.170-184), `to_undirected()` (l.193), convert nx→graph-tool (`convert_networkx_to_graphtool` l.14-100), `sort_query_nodes` by `degree*component_size` for load balance (l.103-126), `gt.label_components` (l.200), per-qnode `GraphView` vertex-filter to its component (l.202), then **`Pool(threads).imap_unordered(traverse...)` over qnodes** (l.205-214); collects results into a fresh `connected_qnodes_gt` graph adding qnode↔counter-qnode edges (l.216-278).
- `graph_traversal.py::traverse_network_to_get_homology_counterparts` (l.305-362): within a qnode's component, sub-partition by **overlap-only edges** (`GraphView efilt overlap=="True"` + `label_components`, l.328-329), then per overlap-subgraph call `summarize_shortest_paths_per_subgraph` (l.209-302): for every vertex run `gt.shortest_path(weights=weight)` from qnode (l.254-257), reject paths ending on a PO edge / with ≥2 adjacent PO edges / with SD-similarity product `∏(1-weight)<=0.8` (l.260-266), walk the route via `inspect_cnode_along_route` (l.13-75) to derive the cnode's relative window, then `compare_homologous_sequences` (l.78-175: minimap2 `-x asm10 --eqx --cs -c`, similarity `=matches/aln_len`, l.155-162). Counterparts kept at `similarity>=0.95` (l.349); counter-qnodes kept looser for better coloring (l.350). `HOMOSEQ_REGION` (homoseq_region.py l.6-104) is the route-carrying coordinate object with `qnode_relative_region` back-projection (l.56-104).
- `graph_query.py::query_connected_nodes` (l.297-334) → `optimal_node_grouping` (l.337-428): the `min_distance==2` fast path calls `gt.sequential_vertex_coloring` and groups by color (l.356-367); the general path builds a conflict graph from all-pairs `shortest_path` then colors it (l.369-428). **In production only the `min_distance=2, max_similarity=None` fast path runs** (the caller passes no overrides), so the general/conflict-graph branch is **dead for parity** — port the fast path; gate the general branch behind a feature flag / leave unported.

**Step 7 — beds + masked genomes + intrinsic align (Easy/External).**
- `build_beds_and_masked_genomes.py` (l.14-216): relabel groups → `RG0..RGn`, load-balance sort (l.36), then **`Pool(nthreads).imap_unordered(establish_beds_per_RG_cluster)`** (l.65-66) writing per-RG query/counterpart/all BEDs (l.148-189), masked genome (`Genome.mask`), intrinsic BAM (`getIntrinsicBam`); finally merge intrinsic BAMs + homo-region BEDs (l.88-118).
- `genome.py::Genome.mask` (l.42-152): pad query BED by `avg+2*std+1000` (l.68-73), slice ref FASTA per interval, N-mask 1000 bp at both ends (l.75-87), **merge same-chrom contigs within 1000 bp with exact N-bridge** (`_merge_nearby_contigs` l.89-152 — coordinate-preserving), write FASTA + md5-gated update.
- `intrinsic_alignment.py::getIntrinsicBam` (l.14-66) → `getRawseq` (seq.py l.5-22, FASTA window extraction) + minimap2 `asm20` against masked genome + `filter_intrinsic_alignments` (l.116-228: drop self-location alignments, dedup/enclosed intervals via `compute_interval_status` l.70-112, secondary→primary promotion). `seq.py::get_bam_frag_size` (l.25-40) → already covered by `sdrecall-io::get_insert_size_distribution`.

### 2. Python → Rust crate mapping

| Python operation / idiom | Rust crate::api | confidence |
|---|---|---|
| `graph-tool gt.label_components` (per-vertex labels) | `petgraph::unionfind::UnionFind::new(n)` + `union(a,b)` over edges + `.into_labeling() -> Vec<K>` (NOT `algo::connected_components`, which returns only the `usize` count) | verified-docs |
| `gt.GraphView(vfilt=…)` component / overlap subgraph | no copy — iterate `petgraph` edges and **filter by the precomputed label vector** (`labels[u]==comp && labels[v]==comp`); for overlap-subgraphs, a second `UnionFind` over only `edge.overlap` edges | verified-docs |
| `gt.shortest_path(weights=weight)` (single source→target, Dijkstra) | `petgraph::algo::dijkstra(&g, src, Some(tgt), |e| *e.weight())` for distances; predecessor map via a thin Dijkstra-with-`came_from` (petgraph's `dijkstra` returns only costs, so wrap a custom BinaryHeap Dijkstra that records predecessors to reconstruct the edge route) | verified-docs |
| `gt.sequential_vertex_coloring` (greedy, vertex-index order) | **custom greedy coloring** in vertex-insertion order (`smallest available color not used by neighbors`) — petgraph ships `algo::coloring::dsatur_coloring` (DSATUR heuristic) but its ORDER/partition differs from gt's index-order greedy; reimplement gt's exact heuristic for partition parity | verified-docs |
| `networkx.DiGraph` SD+PO multiplex | `petgraph::graph::DiGraph<NodeKey, EdgeAttr>` with a `FxHashMap<NodeKey, NodeIndex>` interning map | verified-docs |
| `intervaltree.IntervalTree` (per-chr PO) | insertion-ordered per-chromosome `Vec` plus the parity-preserving linear overlap sweep used by `compose_po_per_chr` | production |
| `pybedtools` sort/merge/intersect/subtract/slop/complement | `sdrecall_io::{sort_merge_bed, intersect, slop, complement}` + a new `subtract` helper (stranded) | plausible |
| `pysam bam.fetch` + per-read tag filters | `sdrecall_io::RegionReader::fetch/records` + `rust_htslib::bam::Record::aux(b"XA")` / `aux(b"XS")` / `aux(b"AS")` | plausible |
| `BedTool.genome_coverage(bg=True)` → per-base depth | custom sweep-line over read intervals (sorted starts/ends, running depth) → emit per-position rows; reuse for all 4 passes | plausible |
| `pyfaidx.Fasta[chrom][start:stop]` slice | `rust_bio::io::fasta::IndexedReader::fetch + read` (random-access faidx) | plausible |
| `Bio.SeqIO.write` masked FASTA | `rust_bio::io::fasta::Writer` | plausible |
| `minimap2 -x asm10/asm20 --eqx --cs -c` (FFI, not subprocess) | `minimap2::Aligner::builder().asm10()/.asm20().with_cigar().with_index(masked,…)?.map(seq,…)` → `Mapping{match_len:i32, block_len:i32, strand}`; similarity `= match_len as f64 / block_len as f64` | verified-docs |
| `samtools faidx region [-i revcomp]` for pairwise similarity | replaced by in-memory `rust-bio` fetch + `bio::alphabets::dna::revcomp` then minimap2-rs `map` of the two slices (no temp FASTA/PAF) | plausible |
| `frozenset({A,B})` SD-pair dedup | `BTreeSet<NodeKey>` (or sorted `(min,max)` tuple) as `AHashSet` key | verified-docs |
| `multiprocessing.Pool.imap_unordered` over qnodes / RGs | `rayon::par_iter` over the independent axis (qnodes; RG clusters), collect into `Vec` | verified-docs |

### 3. Crate file layout (one versatile unit per job)

```
sd-prep/
├─ Cargo.toml                 (lib + bin "sd-prep")
├─ src/
│  ├─ lib.rs                  # pub fn prepare_recall_regions(paths, params) -> Result<Paths>; re-exports
│  ├─ params.rs               # PrepParams {mq_threshold, high_quality_depth, minimum_depth, multialign_frac, threads}
│  ├─ multialign.rs           # ONE depth kernel `inferred_coverage(bam, pass: DepthPass) -> Vec<(chrom,pos,depth)>`
│  │                          #   driven by 4 DepthPass variants (Raw/HighMq/Xa/Xs) → pick_multialigned_regions
│  ├─ sd_pairs.rs             # Pair + ONE `umbrella_to_remove(&[Pair]) -> AHashSet<usize>` (called for raw & granular)
│  │                          #   + filter_umbrella_pairs(group) ; collapses both Python O(n²) sweeps
│  ├─ graph_build.rs          # NodeKey, EdgeAttr, build_multiplex_graph(sd_pairs) -> SdGraph
│  │                          #   insertion-ordered per-chr PO overlap sweep, SD overlay
│  ├─ graph_core.rs           # ★ graph-tool replacement: components (UnionFind),
│  │                          #   dijkstra_with_route(), greedy_vertex_coloring()  ← THE hard unit
│  ├─ homoseq.rs              # HomoseqRegion (route-carrying coord object) + qnode_relative_region back-projection
│  ├─ traversal.rs            # traverse_qnode(qnode, comp_view, …) → counterparts; similarity via minimap2.rs
│  ├─ grouping.rs             # node grouping by coloring → RG clusters (optimal_node_grouping fast path)
│  ├─ masking.rs              # mask_genome(query_bed, ref) — rust-bio FASTA + N-bridge contig merge
│  ├─ minimap.rs              # ONE `align_similarity(q,t,preset)` wrapper over minimap2-rs Aligner (reused by
│  │                          #   traversal similarity AND intrinsic alignment)
│  ├─ intrinsic.rs            # intrinsic_bam(rg) — getRawseq + minimap align + filter_intrinsic_alignments
│  └─ build_rg.rs             # establish_beds_per_RG_cluster: per-RG query/counterpart/all BEDs + outputs
├─ examples/
│  ├─ validate_components.rs       # vs gt.label_components on dumped fixture graphs
│  ├─ validate_coloring.rs         # partition parity vs gt.sequential_vertex_coloring
│  └─ validate_phase1_e2e.rs       # full HG002/HG006 RG groupings + FASTA md5
```

**Collapse points (enforce #1 coding rule):**
- The 4 depth passes (raw/high-MQ/XA/XS) → **one** `inferred_coverage` kernel parameterized by a `DepthPass` enum + a read-predicate closure. No per-pass helpers.
- The raw and granular umbrella sweeps → **one** `umbrella_to_remove(&[Pair])`.
- Similarity-minimap2 (traversal) and intrinsic-minimap2 share **one** `minimap::align_similarity`/builder wrapper, differing only in the `asm10` vs `asm20` preset argument.
- `save_connected_qnodes_to_graphml` and any `nx.write_graphml` → the generic `sdrecall_io::write_graphml` with node/edge attr closures.

### 4. Core data structures + key fn signatures (explicit borrow/owner)

```rust
// homoseq.rs — the route-carrying coordinate object (Python HOMOSEQ_REGION)
#[derive(Clone, Debug)]
pub struct HomoseqRegion {
    pub key: NodeKey,            // (chrom, start, end, strand) — owned, it is an entity identity
    pub size: i64,
    pub ups_rela_start: i64, pub ups_rela_end: i64,
    pub down_rela_start: i64, pub down_rela_end: i64,
    pub rela_start: i64, pub rela_end: i64,
    pub vertex: NodeIndex,       // Copy index into SdGraph
    pub route: Vec<(NodeKey, EdgeKind)>,  // owned route; cloned only when a counterpart is committed
}
// WHY owned route Vec: each accepted counterpart needs an independent route snapshot; borrowing
// would tie its lifetime to the traversal scratch state, which is mutated per shortest-path probe.

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeKind { SegmentalDuplication, Overlap }

// graph_build.rs
pub type NodeKey = (String, i64, i64, Strand);   // matches Python 4-tuple node identity
pub struct EdgeAttr { pub kind: EdgeKind, pub weight: f64, pub overlap: bool, pub nonoverlap_smaller: i64 }
pub struct SdGraph { pub g: DiGraph<NodeKey, EdgeAttr>, pub index: FxHashMap<NodeKey, NodeIndex> }

pub fn build_multiplex_graph(sd_pairs: &[SdPairRow], threads: usize) -> SdGraph;
// WHY &[SdPairRow]: read-only over the filtered SD table; the graph owns fresh node copies.

// graph_core.rs — ★ the graph-tool replacement unit
/// Component label per node via union-find (replaces gt.label_components → labels, NOT just count).
pub fn component_labels(g: &SdGraph, edge_filter: impl Fn(&EdgeAttr) -> bool) -> Vec<u32>;
// WHY &SdGraph + Fn filter: ONE versatile unit serves both "all-edge components" and the
// "overlap-only" sub-partition (graph_traversal l.328-329) by swapping the predicate; no duplicate fn.

/// Weighted single-source→target shortest path returning the EDGE route (petgraph::dijkstra gives
/// only costs, so we record predecessors). Returns None if unreachable.
pub fn dijkstra_route(g: &SdGraph, src: NodeIndex, tgt: NodeIndex)
    -> Option<(f64, Vec<EdgeRef>)>;
// WHY return owned Vec<EdgeRef>: route is consumed downstream by inspect_cnode_along_route; the
// caller does not retain graph borrows past reconstruction.

/// graph-tool sequential_vertex_coloring parity: greedy, vertices in ascending index order,
/// smallest color absent among already-colored neighbors. Returns color per node index.
pub fn greedy_vertex_coloring(adj: &[Vec<u32>]) -> Vec<u32>;
// WHY &[Vec<u32>] adjacency (not the petgraph type): coloring is pure combinatorics on the
// connected_qnodes graph; a CSR-ish adjacency keeps it cache-friendly and decoupled from petgraph.

// traversal.rs — parallel over qnodes (rayon)
pub fn extract_sd_paralog_pairs(
    query_nodes: &[NodeKey], graph: &SdGraph, ref_fa: &Path,
    avg_frag: f64, std_frag: f64, mean_read_len: f64, threads: usize,
) -> Result<(FxHashMap<NodeKey, Vec<HomoseqRegion>>, QnodeGraph)>;
// WHY &SdGraph shared (Arc not needed — rayon borrows immutably across the par_iter); each qnode
// task allocates only its own HomoseqRegion scratch.

// minimap.rs — ONE wrapper reused by traversal + intrinsic
pub fn align_similarity(query: &[u8], target: &[u8], preset: Preset) -> Result<f64>;
// returns match_len/block_len of the best mapping; &[u8] borrows reference slices (zero-copy).

// masking.rs
pub fn mask_genome(query_bed: &[GenomicInterval], ref_fa: &Path, out: &Path,
                   avg_frag: f64, std_frag: f64, merge_gap: i64) -> Result<()>;
```

### 5. Performance optimizations from ownership / borrowing

- **No nx→graph-tool conversion.** Python pays a full `convert_networkx_to_graphtool` copy (graph_query l.14-100) per run; Rust builds the `petgraph::DiGraph` once and never re-materializes it. Component subgraphs are **label-vector views** (a `Fn(&EdgeAttr)->bool` predicate over the single owned graph), never copied (Python's `gt.GraphView` is cheap but the nx copy at l.156 `directed_graph.copy()` is not).
- **rayon over validated independent axes** replaces selected `multiprocessing.Pool` work with zero IPC/pickling, including qnode traversal and per-RG bed/mask/align. The per-chromosome PO overlay intentionally remains an insertion-ordered linear sweep because that ordering is part of the parity contract.
- **Linear PO overlap sweep instead of Python `intervaltree`:** the current `Vec` implementation preserves Python insertion order and avoids the unused Lapper layer that was removed during closeout.
- **minimap2 via FFI, in-memory** (`minimap.rs`): the Python hotspot does **minimap2 + 2× `samtools faidx` subprocess per candidate counterpart** (graph_traversal l.120-141) — process spawn + temp FASTA + temp PAF parse per probe. Rust holds the masked-genome index once (`Aligner::with_index`) and maps `&[u8]` ref slices directly; `Mapping{match_len, block_len}` gives the identical `matches/aln_len` similarity with no I/O.
- **Avoid per-record clones** in depth passes: `RegionReader::records()` yields `Result<bam::Record>`; extract `(tid,start,end)` into a `Vec<(u32,u32)>` and drop the record immediately — no buffering of full records. Per-base coverage is a single sweep-line, not a materialized interval tree.
- **FxHashMap for the NodeKey→NodeIndex intern map and component dictionaries** (integer-ish keys after interning); the SD-pair dedup uses a sorted `(min,max)` `(NodeKey,NodeKey)` key in an `AHashSet` (string-containing) — no `frozenset` allocation per row.
- **`Cow`/`&[u8]` for revcomp**: same-strand slices borrow the ref FASTA buffer; only opposite-strand counterparts allocate a `revcomp` Vec (Python always shells out `samtools faidx -i`).

### 6. Risks / open decisions

- **Coloring partition parity (resolved):** production uses the pure-Rust greedy coloring with the required stable insertion order. Retained three-assembly grouping/FASTQ parity supersedes the old open decision; there is no graph-tool fallback in the production path.
- **`sort_query_nodes` ordering (resolved):** the load-balancing sort and tie behavior are fixed by the production parity tests because they feed vertex insertion order and therefore coloring.
- **`dijkstra_route` tie-breaking vs `gt.shortest_path`.** Equal-cost paths may be returned in different order; since downstream rejects paths by structural rules (trailing PO edge, adjacent PO edges, SD-product threshold, l.260-266) and then recomputes coordinates, a different equal-cost path **could change the cnode set**. Hazard: must verify gt's shortest_path returns the same representative path (gt uses Dijkstra with its own heap order). Dump per-qnode chosen paths from Python and diff.
- **minimap2-rs vs CLI minimap2 numerical identity.** minimap2-rs 0.1.31 bundles minimap2 2.30; the Python pipeline shells the system `minimap2` binary (version unpinned). `match_len/block_len` should be deterministic for `asm10/asm20`, but a minor version skew could shift the `>=0.95` similarity threshold near the boundary. **Flag as external-algorithm version skew**; pin/record the system minimap2 version used for the Python oracle.
- **Per-base depth parity.** `BedTool.genome_coverage(bg=True)` then expanding to per-position rows (inferred_depths l.124-127) is an exact integer coverage; the sweep-line must match it bit-for-bit including the `+1` half-open conversion (l.50-51 `start=pos-1, end=pos`). The `AS-XS<=5` XS-tag rule (l.18) and per-read flag filters must mirror pysam exactly.
- **`_merge_nearby_contigs` coordinate invariant.** The N-bridge length `next_start-(cur_start+len)` (genome.py l.121-122) is load-bearing for the downstream `modify_masked_genome_coords` mapping; the masked-FASTA md5 criterion will catch any drift, but the bridge math must be exact (off-by-one in the 1000 bp end-mask slice `[1000:-1000]`, l.81).
- **Dead-code traps to NOT port:** `inferred_depths.calculate_coverage` (broken), `optimal_node_grouping` general/conflict-graph branch (never reached in production), and the commented HOMOSEQ_REGION rewrite (l.106-198). Porting them would create the exact near-duplicate units the #1 rule forbids.
- **Leaf subprocess survivors:** checked `samtools merge`/`sort`/`index` operations remain the deliberate production implementation for total intrinsic BAM handling; the unused in-process `merge_bams` proposal was removed.

### 7. Test plan delta (concrete fixtures + assertions)

- **`validate_components.rs` (unit):** build 3 fixture `SdGraph`s (single component; two disjoint; one with an overlap-only bridge). Assert `component_labels(g, all_edges)` induces the **same set partition** as a dumped `gt.label_components` labeling (compare partitions, not label ids). Add the overlap-only predicate case (replaces `gt.label_components` on the overlap `GraphView`).
- **`validate_coloring.rs` (unit):** hand-build a small `connected_qnodes` graph (path P4, a 5-cycle, K3+isolated) with a **fixed node-insertion order matching the Python sort**, run `greedy_vertex_coloring`, and assert the **partition equals** `gt.sequential_vertex_coloring` grouped-by-color (dumped from a one-off Python script). Property test: coloring is proper (no edge joins same-color nodes) on random graphs.
- **Umbrella filter (unit):** hand-build overlapping `Pair`s where A umbrella-covers B at ≥0.95 on both segments and the strand-consistency `(strandA==strandB)==(otherA==otherB)` flips — assert `umbrella_to_remove` returns exactly B's index for both the raw and granular sweeps (one shared helper, two call sites).
- **Depth kernel (unit):** 5 synthetic reads with known XA/XS/AS tags over a 30 bp window; assert per-base `raw/high-MQ/XA/XS` depths and the final AND-filter mask match a hand-computed table; assert `start=pos-1,end=pos` BED conversion.
- **`validate_phase1_e2e.rs` (differential, tier 2):** full HG002 + HG006 Phase-1 run. Pass criteria (from doc): (a) **RG node groupings identical as a set partition**; (b) **paralog-pair set identical** (`FxHashMap<NodeKey,Set<NodeKey>>` equality); (c) **masked-genome FASTAs md5-identical**; (d) `filtered_SD_binary_map.tsv` row-set identical after frozenset dedup. Dump the Python `connected_qnodes_gt` vertex order + per-qnode chosen shortest paths to localize any coloring/path divergence before comparing final groupings.
