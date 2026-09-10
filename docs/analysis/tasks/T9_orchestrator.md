# T9 — `sdrecall` orchestrator (retire Python + PyO3)

**Crate:** `sdrecall` (new — bin)
**Status (2026-07-17):** Production-integrated. `run`, `prepare`, and `realign` execute the full Rust pipeline; retained t2t/hg19/hg38 runs pass semantic parity. PyO3/Python binding cleanup is complete; formal HG006 Rust-versus-Python validation remains.
**Depends on:** T0–T8
**Replaces (Python):** the `SDrecall` CLI, `realign_and_recall.py`, `misalignment_elimination.py` orchestration, and the thin wrappers absorbed here: `realign_recall/realign_per_RG.py`, `slice_bam_by_cov.py` (island detection), `stat_realign_group_regions.py`, and `annotate_HP_tag_to_vars.py` (if not in T5).

## Goal & scope boundary

The top-level Rust binary threads stage libraries and manages the complete pipeline. The current production contract is **no Python interpreter**; approved `samtools`, `bcftools`, and `minimap2` subprocesses may remain for I/O-heavy and algorithmic operations.

This task also **retires `py-bindings`**: as each stage's orchestration moves in-process, its PyO3 entry point is deleted.

## Data flow

Mirrors the pipeline spine (see [`../RUST_MIGRATION_PLAN.md`](../RUST_MIGRATION_PLAN.md) flowcharts), with all glue in Rust:

```
prepare:  sd-prep
realign:  region-prep → read-extraction → minimap2 (minimap2-rs FFI) + variant call (FFI / leaf bcftools)
          → merge + markdup (rust-htslib; markdup reimpl or leaf) → slice islands (rust-htslib)
          → fp-control (per island, rayon) → filter BAM (rust-htslib) + variant call
          → vcf-ops (priority merge, rust-htslib) → subset to target (rust-htslib)
post:     vcf-ops (inhouse-common) → vcf-ops (merge w/ conventional)
```

Per-island parallelism uses `rayon` (replaces `multiprocessing.Pool` with `spawn`); the thread budget mirrors the current `job_num × threads_per_job = total_threads` rule.

## Dependencies

- All stage libs (`sd-prep`, `region-prep`, `read-extraction`, `fp-control`, `vcf-ops`, `nm-stats`) + `sdrecall-utils`/`sdrecall-io`.
- `clap` 4.6.1 (CLI), `rayon` 1.8 (parallelism).
- External production tools: minimap2/minimap2-rs for alignment, bcftools for calling/normalization, and selected samtools BAM operations. Core orchestration and most parsing remain in-process Rust. See the July policy.

## Performance bottleneck / rationale

Orchestration correctness + parallelism. The win is the end state: one binary, no Python, no PyO3, in-process data flow between stages.

## Tests

### Unit (tier 1)
- Per-stage CLIs already validated in their own tasks; here, test the thread-budget calc and subcommand wiring.

### Differential vs Python (tier 2)
- Full-pipeline run on HG006: final recall VCF vs the Python pipeline's final VCF.

**Pass criterion:** `final_recall` VCF record-identical (normalized) vs Python on HG006. Once green, delete the corresponding Python modules + `py-bindings` entries stage-by-stage.

**Data:** HG006 full run.

## Progress
- [x] Scaffold `sdrecall` bin with `run`/`prepare`/`realign`.
- [x] Port `Paths`/`SDrecallPaths` behavior into `sdrecall::paths`.
- [x] Thread all production stages through the Rust orchestrator.
- [x] Implement per-island parallelism, shared resource leases, and bounded helper threading.
- [x] Absorb the Python orchestration wrappers; approved external bioinformatics commands remain behind Rust wrappers.
- [ ] Full-pipeline HG006 differential
- [x] Remove the production PyO3/Python binding surfaces, `cdylib` targets, build metadata, checked-in wheel, and pyo3/numpy/pyo3-log dependencies.

## Review findings (2026-06-11)

From the migrated-code review — full detail + IDs in [`../REVIEW_FINDINGS.md`](../REVIEW_FINDINGS.md). One theme matters specifically at the single-process orchestrator stage.

- **ROB-1 (fixed 2026-06-11):** `panic="abort"` has been removed from the workspace root profile (and the dead override in `read_extraction/Cargo.toml` deleted), so the panic strategy is now the default `unwind`. A panic in any stage no longer aborts the whole binary — it unwinds and is catchable. **Decision recorded:** unwind, so that under the single-process orchestrator one island's panic can be isolated per-`rayon`-task (catch with `std::panic::catch_unwind` or `rayon`'s propagation) instead of taking down the run. The earlier "keep abort vs switch to unwind" question is settled in favour of unwind.
  - **ROB-2 (fixed):** `read_extraction` returns a normal Rust error when a read name is not valid UTF-8.
  - **ROB-3 (fixed):** `M`-CIGAR input returns typed `CigarError::UnsupportedMatchOp`; per-island panic isolation remains as a last-resort containment boundary, not the normal error path.
