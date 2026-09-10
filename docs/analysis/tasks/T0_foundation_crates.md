# T0 — Foundation crates: `sdrecall-utils` + `sdrecall-io`

**Crates:** `sdrecall-utils` (lib), `sdrecall-io` (lib)
**Status (2026-07-17):** Production-integrated. Both crates are in the 11-crate workspace and support the end-to-end orchestrator; `Paths` ultimately remained in `sdrecall`, and selected `samtools`/`bcftools` operations are allowed by the current production policy.
**Depends on:** —
**Replaces (Python):** scattered helpers in `src/utils.py`, path logic in `src/const.py` (`SDrecallPaths`), and ad-hoc BAM/BED/VCF parsing duplicated across modules.

> The detailed sections below preserve the June interface design. Where that design conflicts with the live workspace, the July status above and current source are authoritative.

## Goal & scope boundary

Two foundation libraries every stage crate depends on, so that **no stage re-implements file I/O or shared types** (directly enforces the user's #1 coding rule against overlapping helpers). This is also where the **inter-crate interface contract** is defined: the common in-memory structs that flow between stages in-process, and the common file formats (BAM/BED/VCF/GraphML/TSV) that each stage's CLI reads and writes.

- `sdrecall-utils`: pure types + errors + logging. No I/O, no heavy deps.
- `sdrecall-io`: shared format readers/writers built around `rust-htslib`, interval utilities, GraphML, and TSV support. The production pipeline may also use approved external-tool wrappers.

Out of scope: any domain compute (that lives in stage crates).

## What lives here (Python origins)

Neither crate is a pipeline step — they are the shared toolbox, collecting plumbing currently scattered across `src/` and copy-pasted as pysam/pybedtools idioms in every module.

**`sdrecall-utils`** — pure logic, decisions, labels; no file engine:

| Rust home | Python origin |
|-----------|---------------|
| Where every input/output/temp file for a run lives | `src/const.py::SDrecallPaths` (all `*_path()` methods) |
| Logging (colored console + per-subprocess capture → `(ok, result, log)` tuples used by the `imap_*` workers) | `src/log.py` (`configure_logger`, `log_command`, `error_handling_decorator`, `init_logger`) |
| Parallelism budget (jobs × threads) | `src/utils.py::configure_parallelism` |
| Output-freshness check, NaN detection | `src/utils.py` (`is_file_up_to_date`, `na_value`) |
| Warning silencing | `src/suppress_warning.py` |
| Shared data types: interval, region key, hap id, qname index | currently bare tuples/dicts passed around |

**`sdrecall-io`** — shared file-format utilities centered on in-process Rust libraries. The table records the original replacement targets; current production policy also permits approved external `samtools`/`bcftools` wrappers:

| Rust home (in-process) | Library | Python origin it replaces |
|------------------------|---------|---------------------------|
| BAM/CRAM read per region (fetch/iterate + QC filters) | rust-htslib | pysam idioms in `bam_ncls.py`, `slice_bam_by_cov.py`, `cal_edge_NM_values.py`, `annotate_HP_tag_to_vars.py` (already done in Rust `read_extraction`) |
| BAM merge / sort / index / validity | rust-htslib | `src/utils.py::merge_bams` + `shell_utils.sh` (was samtools) |
| VCF/BCF read (sorted cursor), write / concat / sort | rust-htslib | `identify_common_vars.py`, `merge_variants_with_priority.py`, `src/utils.py::combine_vcfs` (was bcftools) |
| BED read/write + intersect / merge / slop / complement | bedrs | `src/utils.py` (`sortBed_and_merge`, `merge_bed_files`), pybedtools (was bedtools) |
| GraphML write (read via quick-xml only if needed) | petgraph-graphml | graph-tool `.save`/`.load` in `prepare_recall_regions.py`, `graph_query.py` |
| Fragment-size estimation | rust-htslib | `src/insert_size.py::get_insert_size_distribution`; `preparation/seq.py` |
| Temp file; md5-compare-and-replace | std + a hash crate | `src/utils.py` (`prepare_tmp_file`, `update_plain_file_on_md5`) |

> `src/utils.py::executeCmd` (the shell-runner) does **not** survive migration — its callers (samtools/bcftools/bedtools wrappers) become the in-process libraries above. The genuinely external *algorithms* (alignment, variant calling) are handled at the orchestrator level — see the External tools & libraries policy in the plan.

**Why two crates, not one?** `sdrecall-io` needs the heavy, slow-to-build rust-htslib engine (the one requiring `LIBCLANG_PATH` + the OpenSSL env vars in the build notes); `sdrecall-utils` needs nothing heavy. Splitting them lets a math-only stage (e.g. `region-prep`) or a quick test depend on `sdrecall-utils` without compiling htslib. If that isolation ever proves not worth it, they can be merged into one `sdrecall-common` crate later.

## Interface contract (the shared vocabulary)

`sdrecall-utils` types (initial set, grown as stages land):

| Type | Meaning |
|------|---------|
| `GenomicInterval { chrom: String, start: i64, end: i64, strand: Strand }` | half-open BED interval |
| `RegionKey` | hashable (chrom,start,end) for per-region maps |
| `HapId(i32)`, `QnameIdx(u32)` | newtypes for haplotype id / qname index |
| `Paths` | port of `SDrecallPaths` (all derived I/O paths for a run) |
| `SdError` (thiserror) | crate-wide error enum |

`sdrecall-io` API surface (all in-process):

| Function | Format — library |
|----------|------------------|
| `read_bam` / `BamReader` (indexed + streaming) | BAM/CRAM — rust-htslib |
| `write_bam` / sort / index / merge | BAM — rust-htslib (no samtools) |
| `read_bed` / `write_bed` + intersect/merge/slop/complement | BED — bedrs |
| `read_vcf` / `write_vcf` (sorted cursor) + concat / sort | VCF/BCF — rust-htslib (no bcftools) |
| `write_graphml` (read via `quick-xml` only if needed) | GraphML — petgraph-graphml |
| `read_tsv` / `write_tsv` | TSV |

## Data flow

Library only — no pipeline position. Provides the types and codecs that all other tasks consume.

## Dependencies

- `sdrecall-utils` crates: `serde` 1.0, `thiserror` 1.0, `log` 0.4, `rustc-hash` 1.1, `ahash` 0.8 — no heavy deps.
- `sdrecall-io` crates: `rust-htslib` 0.47.0 (BAM/CRAM/BCF/VCF), `bedrs` 0.2.26 (BED + interval ops; its `rust-htslib` feature turns BAM records into intervals), `petgraph` 0.8.3 + `petgraph-graphml` 5.0.0 (GraphML write), `quick-xml` 0.40.1 (GraphML read, only if needed).
- Foundation: root of the dependency DAG; `sdrecall-io` depends on `sdrecall-utils`.
- External tools: selected checked `samtools` and `bcftools` leaf operations remain by the July production policy; core parsing and interval work stays in process. Full pinned list: see the plan's *Dependency versions (pinned)* section.

## Performance bottleneck / rationale

Not hot itself. Its existence removes two recurring costs: (1) duplicated parsing/encoding code, and (2) — once `fp-control` (T4) is fused — the ability to open and parse a BAM **once** and share the reader/records, instead of today's double parse (phasing-graph + haplotype-inspection each open the BAM).

## Tests

### Unit (tier 1)
- Round-trip read→write→read for each format on small fixtures; assert structural equality.
- `GenomicInterval` overlap/containment/merge edge cases.
- `Paths` derives the same strings as the Python `SDrecallPaths` for a known run config.

### Differential vs Python (tier 2)
- Re-emit HG002 fixture BED/VCF through `sdrecall-io` and diff against the Python-written equivalents.

**Pass criterion:** byte-identical (BED/TSV) or record-identical (VCF normalized) re-emit on HG002 fixtures; `Paths` strings match `SDrecallPaths` exactly.

## Progress
- [x] Scaffold both crates in the workspace.
- [x] Implement shared errors, geometry, logging, parallelism, fatal handling, and resource controls in `sdrecall-utils`.
- [x] Implement BAM/BED/VCF readers and writers plus insert-size support in `sdrecall-io`.
- [x] Implement GraphML + TSV codecs.
- [x] Exercise the foundation layer through focused tests and three-assembly end-to-end production parity.
- [x] Port path derivation for the production pipeline; implementation ownership settled in `sdrecall::paths` rather than `sdrecall-utils`.

## Review findings (2026-06-11)

From the migrated-code review — full detail + IDs in [`../REVIEW_FINDINGS.md`](../REVIEW_FINDINGS.md). These confirm T0 targets the right redundancy: the dominant duplication in the existing crates is exactly what these two foundation crates absorb.

**→ `sdrecall-io` absorbs:**
- **DUP-1 (HIGH):** "read BAM → drop noisy reads → build interval index" is implemented **three times** — `build_phasing_graph/src/bam_reading.rs:16-365`, `haplotype_inspection/src/bam_lappers.rs:39-618`, and (lighter) `read_extraction/src/lib.rs:42-103`. Same `samtools collate`, same qname-streaming, same noisy-read filter. Collapse to one reader + one filter here.
- **DUP-3 (MED):** interval-overlap done three ways (hand-rolled `structs.rs:115`, `rust-lapper` `bam_lappers.rs:663`, endpoint sweeps in `graph_builder.rs`). Pick one interval backend.
- **HYG-6 (LOW):** `bam_lappers.rs:103` leaks one file handle per island via `mem::forget`; when this logic moves here, prefer a plain `drop` (htslib already dup'd the fd) — verify with a pipe-read + fd-count test.

**→ `sdrecall-utils` absorbs:**
- **DUP-2 (HIGH):** the per-read vector toolkit (`extract_hap_vector`, `extract_error_vector`, variant counting, the A/T/C/G/N encoder ×3, read-id, caches) is duplicated across `haplotype_determination.rs` ↔ `pairwise_read_inspection.rs`.
- **DUP-4 (MED):** shared-SNV / position→base mapping duplicated (`haplotype_determination.rs:1090` ↔ `identify_misaligned_haps.rs:389`).
- **DUP-5 (LOW):** small interval-merge / sort-dedup helpers overlap.

**Resolved — DUP-2a (HIGH):** the production encoders now reject `M` CIGAR operations through typed `CigarError::UnsupportedMatchOp` propagation, and the golden indel encoding (`INDEL_UNIT=10`, `HAP_PAD=-20`, compound SNV/insertion decoding) is covered by focused oracle tests. The former panic/tolerate divergence is historical, not open work.
