# SDrecall Rust Workspace

This directory is a **Cargo workspace** holding the Rust port of the SDrecall
pipeline. The production path is a single `sdrecall` Rust binary with no Python
interpreter or PyO3 bindings. Core stages call one another as Rust libraries;
`minimap2`, `bcftools`, and `samtools` remain supported external tools for
alignment, variant calling, and selected BAM/VCF operations.

This README is the **developer/operator reference** for the workspace: how the
crates fit together, how the orchestrator threads them, and a copy-pasteable
command + input/output contract for every module.

> For the **user-facing** science tool (installation, biology, the production
> `SDrecall run` CLI) see the top-level [`../README.md`](../README.md).
> For migration status, per-task contracts and validation evidence see the hub
> [`../docs/analysis/RUST_MIGRATION_SUMMARY.md`](../docs/analysis/RUST_MIGRATION_SUMMARY.md).

---

## 1. Design in one picture: a workspace of `lib + bin` crates

Every pipeline stage is its own crate built as **both a library and a thin
binary**:

- the **library** is the in-process API the `sdrecall` orchestrator calls
  directly (function in → values out, no subprocess, no temp files on the hot
  path);
- the **binary** is a thin `files-in → files-out` CLI so the stage can be run,
  tested and differentially validated against the Python output on its own.

Two **foundation crates** hold the shared plumbing so no stage re-implements it:

| Crate | Role | CLI? | Replaces (Python) |
|-------|------|------|-------------------|
| `sdrecall-utils` | shared types, errors, logging, parallelism math — **no file I/O** | — (library) | `utils.py` / `log.py` / `const.py` helpers |
| `sdrecall-io` | in-process BAM/BED/VCF/GraphML/TSV + insert-size I/O | — (library) | `utils.py` I/O, `insert_size.py` |
| `read_extraction` | BAM → FASTQ (region/multi-align filtered) | — (library) | `realign_recall` read extraction |
| `haplotype_inspection` | consensus / similarity / BILC solve | — (library) | `fp_control/identify_misaligned_haps.py` |
| `phasing` | BAM → phasing graph + weight matrix → GCE partition → HP-tagged BAM | bin = phaser + diff harness | `fp_control/graph_build.py` + `phasing.py` + `gce_algorithm.py` |
| `region-prep` | per-RG fc/nfc realignment-region projection | `region-prep` | `prepare_masked_align_region.py` |
| `fp-control` | **fused** Phase-2c: graph → phasing → inspect → BILC | `fp-control` | `realign_filter_per_cov.py` wiring |
| `vcf-ops` | priority VCF merge + inhouse-common annotation | `vcf-ops` | `merge_variants_with_priority.py`, `identify_common_vars.py` |
| `nm-stats` | per-BAM NM (edit-distance) Poisson cutoff | `nm-stats` | `cal_edge_NM_values.py` |
| `sd-prep` | Phase-1 SD graph + region prep (**partial CLI**) | `sd-prep` (graph/mask) | `prepare_recall_regions.py` + `preparation/*` |
| `sdrecall` | **top-level orchestrator** — threads all stages | `sdrecall` | `SDrecall` CLI + `realign_and_recall.py` + `misalignment_elimination.py` |

The obsolete PyO3 bindings, Maturin manifests, extension build scripts, and
checked-in wheel were removed on 2026-07-17. `read_extraction` and
`haplotype_inspection` now expose only the Rust library APIs used by the
orchestrator. The graph builder absorbed from the former `build_phasing_graph`
crate lives in `phasing`.

---

## 2. How the orchestration works

`sdrecall` (the binary in [`sdrecall/`](sdrecall/)) is the **conductor**. It
depends on every stage crate as a path-dependency (see
[`sdrecall/Cargo.toml`](sdrecall/Cargo.toml)) and calls their **library
functions in-process** — there is no CLI-to-CLI piping and no Python on the
control path. Its own glue logic lives in four local modules:

| Module | Job | Ports |
|--------|-----|-------|
| [`sdrecall/src/pipeline.rs`](sdrecall/src/pipeline.rs) | the stage call-graph (the spine below) | `realign_and_recall.py` + `misalignment_elimination.py` |
| [`sdrecall/src/tools.rs`](sdrecall/src/tools.rs) | subprocess wrappers for minimap2/bcftools/samtools | `shell_utils.sh` |
| [`sdrecall/src/island.rs`](sdrecall/src/island.rs) | coverage-island detection + per-island BAM slicing | `slice_bam_by_cov.py` |
| [`sdrecall/src/rg_discovery.rs`](sdrecall/src/rg_discovery.rs) | RG discovery + load-balancing sort | `stat_realign_group_regions.py` |
| [`sdrecall/src/bam_filter.rs`](sdrecall/src/bam_filter.rs) | in-process BAM filter by qname set | `realign_filter_per_cov.py` filter loop |

### Pipeline spine

Each `→` is one **in-process library call** unless marked `[ext]` (leaf
subprocess):

```
prepare:  sd_prep::prepare_recall_regions                                  [Phase 1]

realign:  rg_discovery::stat_all_rg_region_size            (which RGs, biggest first)
          ├─ per RG (rayon) ─ region_prep::prepare_masked_align_region_per_rg
          │                   → read_extraction::bam_to_fastq
          │                   → minimap2 [ext] → bcftools call [ext]
          ├─ samtools merge + markdup [ext]      → pooled deduped raw BAM
          ├─ sdrecall_io::concat_sort_vcfs        → pooled raw VCF
          ├─ island::split_bams_into_islands      (samtools depth [ext] + BED math)
          ├─ nm_stats::nm_distribution_poisson    (NM cutoff)
          ├─ per island (rayon, catch_unwind) ─ fp_control::run_fp_control
          │                   → bam_filter::filter_bam_by_qnames
          │                   → bcftools call [ext]
          ├─ samtools merge + sdrecall_io::concat_sort_vcfs   (clean BAM + VCF)
          └─ vcf_ops::merge_with_priority (raw vs clean) → bcftools view -R [ext]  → SDrecall VCF

post:     vcf_ops::annotate_inhouse_common (cohort)  → vcf_ops::merge_with_priority (conventional)
```

This mirrors the Python `SDrecall` / `realign_and_recall.py` /
`misalignment_elimination.py` call graph 1:1.

### Parallelism

`ThreadBudget` (in `pipeline.rs`) mirrors the Python `configure_parallelism`
rule: it splits the total thread count into `num_jobs × threads_per_job`. There
are two **rayon** fan-out points — **per-RG** (realignment) and **per-island**
(fp-control). Each island runs inside `std::panic::catch_unwind`, so one island
panicking is logged and skipped rather than aborting the whole sample.

---

## 3. Build & environment

All crates build inside the `SDrecall` conda env. The env vars let `rust-htslib`
find `libclang` and force the system (conda) `curl`/OpenSSL instead of a vendored
build:

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

cd /paedyl01/disk1/yangyxt/SDrecall-rust-migration/rust_modules

cargo build --release                 # whole workspace
cargo build --release -p vcf-ops      # one stage binary
cargo test                            # all unit tests
cargo run -p sdrecall -- --help       # the orchestrator CLI
```

Run a built stage binary either via `cargo run -p <crate> --release -- <args>`
or directly as `./target/release/<crate> <args>`. The examples below use
`cargo run` form; **placeholder paths are in angle brackets** — substitute real
files.

---

## 4. Module reference (command + inputs + outputs)

### Foundation crates (library only — no standalone CLI)

#### `sdrecall-utils`
Shared types and pure helpers used by every crate; **no file I/O**.

- **Library surface:** `SdError` / `Result` (the workspace error type),
  `configure_parallelism(total, per_job) -> (num_jobs, threads_per_job)`,
  `init_console_logger(level)`, `GenomicInterval`, `Strand`.
- **Consumed by:** every other crate. There is nothing to run on its own.

#### `sdrecall-io`
Shared in-process file I/O plus the checked `samtools`/`bcftools` leaf wrappers
that remain deliberate production dependencies.

- **Library surface (modules):** `bam`, `bed`, `vcf`, `graphml`, `tsv`,
  `insert_size`. Frequently-used fns: `read_bed`, `sort_merge_bed`, `slop`,
  `concat_sort_vcfs`, `get_insert_size_distribution`.
- **Consumed by:** `region-prep`, `sd-prep`, `island.rs`, the VCF concat steps.
  Nothing to run on its own.

---

### Stage crates

#### `region-prep` — per-RG fc/nfc realignment-region projection
```bash
cargo run -p region-prep --release -- \
  --rg_label RG0 \
  --whole_region_bed <RG0_related_homo_regions.bed> \
  --target_bed       <sample_target.bed> \
  --ref_genome       <hg38.fasta> \
  --fc_target_out    <out/RG0.fc_target.bed> \
  --nfc_out_dir      <out/RG0/>
# optional: --subgroups 0,1,2   (default: every FC: subgroup in the whole-region BED)
```
- **Inputs:** the 7-column `{rg}_related_homo_regions.bed` from Phase 1; the
  sample target BED (first 3 cols); the reference (its `.fasta.fai` supplies
  contig sizes for slop).
- **Outputs:** one shared **FC target BED** (`--fc_target_out`); one
  **`{rg}_{sub}.nfc.bed`** per subgroup in `--nfc_out_dir`; prints one CSV row
  per subgroup to stdout: `rg,sub,fc_bed,nfc_bed,fc_size,nfc_size`.

#### `nm-stats` — NM (edit-distance) Poisson cutoff
```bash
cargo run -p nm-stats --release -- \
  --bam <deduped_raw.bam> \
  --conf_level 0.01 --sample_size 3000000 --threads 4
```
- **Inputs:** one coordinate-sorted BAM.
- **Outputs:** a single line to **stdout**: `cutoff<TAB>mean` (the per-read NM
  edit-distance cutoff and the sampled mean) — same shape as the Python dump.

#### `fp-control` — fused Phase-2c FP control (the headline hotspot)
Runs `phasing` (graph build + GCE) → `haplotype_inspection` on one island BAM
in a single in-process call.
```bash
cargo run -p fp-control --release -- \
  --bam        <island.realigned.bam> \
  --intrinsic  <island.intrinsic.bam> \
  --reference  <hg38.fasta> \
  -o           <island.fpcontrol.tsv> \
  --mapq_cutoff 10 --edge_weight_cutoff 0.301 --threads 4
```
- **Inputs:** the per-island realigned BAM; the matching intrinsic
  (reference-sequence) BAM; the reference FASTA (`.fai` required — used by the
  graph's allele-depth mpileup).
- **Outputs:** a TSV `<output>` with header `qname<TAB>label`, where `label` is
  `correct` or `mismap`. An island skipped for being trivial (≤2 haplotypes)
  yields a present-but-empty TSV.

#### `vcf-ops` — priority merge + inhouse-common annotation
Two subcommands.
```bash
# priority merge (query wins; reference fills gaps) — the raw-vs-clean merge
cargo run -p vcf-ops --release -- merge \
  --query_vcf     <recall.raw.vcf.gz> \
  --reference_vcf <recall.clean.vcf.gz> \
  --output_vcf    <recall.merged.vcf.gz> \
  --ref_genome    <hg38.fasta> \
  --added_filter MISALIGNED --qv_tag RAW --rv_tag CLEAN

# inhouse-common binomial annotation against a control cohort
cargo run -p vcf-ops --release -- inhouse-common \
  --query_vcf  <sdrecall.vcf.gz> \
  --cohort_vcf <control_cohort.vcf.gz> \
  --output_vcf <sdrecall.inhouse_common.vcf.gz> \
  --ref_genome <hg38.fasta> \
  --inhouse_common_cutoff 0.01 --conf_level 0.999
```
- **Inputs:** two bgzipped+indexed VCFs (`merge`: query + reference;
  `inhouse-common`: query + cohort) and the reference FASTA.
- **Outputs:** one bgzipped+indexed VCF at `--output_vcf` (path echoed to
  stdout). `merge` tags each record's source (`--qv_tag`/`--rv_tag`) and adds
  `--added_filter` to query-only records; `inhouse-common` adds the
  `--added_filter` (default `INHOUSE_COMMON`) flag to variants judged common in
  the cohort.

#### `sd-prep` — Phase-1 SD graph & masking (**partial CLI**)
The full Phase-1 driver runs in-process via the library
(`sd_prep::prepare_recall_regions`); the binary currently exposes only the two
**validated units** for differential checking — `traverse`/`prepare`
subcommands are intentionally not wired yet.
```bash
# build the multiplex SD+PO graph and report counts + component partition
cargo run -p sd-prep --release -- graph --sd_map <filtered_SD_binary_map.tsv>

# build a masked genome from a query BED and print its md5
cargo run -p sd-prep --release -- mask \
  --query_bed <RG0.query.bed> --ref_fa <hg38.fasta> --out <RG0.masked.fasta>
```
- **Inputs:** `graph` — the 9-column `filtered_SD_binary_map.tsv`; `mask` — a
  query BED + reference FASTA.
- **Outputs:** `graph` — node/edge/SD-edge/component counts + component-size
  histogram to stdout; `mask` — the masked FASTA at `--out` plus its md5
  (for the byte/md5 differential vs Python `RG<n>.masked.fasta`).

#### `phasing` — graph phasing + GCE (**differential harness**, not a stage CLI)
The production phasing is the **library** (`phasing::phase`,
`phasing::qname_partition`), called inside `fp-control`. The binary is a
validation harness that replays per-island dumps written by
`fp_control/diff_dump.py` and asserts the qname partition matches Python up to
relabeling.
```bash
cargo run -p phasing --release -- --path <dump_root>/         # all islands
cargo run -p phasing --release -- --path <dump_root>/island_7 --single
```
- **Inputs:** a dump directory per island (`weight_matrix.npy`, `edges.json`,
  `node_read_ids.json`, `vertex_qname.json`, `phasing_hap_qname_info.json`, …).
- **Outputs:** a per-island `MATCH`/`MISMATCH` report and an `N/total islands
  match` summary; non-zero exit if any island mismatches.

#### `read_extraction` — BAM → FASTQ (library, no CLI)
- **Library surface:** `bam_to_fastq(input_bam, region_bed, r1, r2,
  multi_aligned, threads)`. The orchestrator links this API directly; there is
  no Python feature or extension-module build.
- **Inputs:** an indexed BAM + a region BED. **Outputs:** paired `r1`/`r2`
  FASTQ files of reads overlapping the regions (`multi_aligned` toggles the
  multi-alignment recruitment filter).

#### `phasing` / `haplotype_inspection` — Phase-2c kernels
The compute core of `fp-control`: BAM → phasing graph + weight matrix → GCE
partition (`phasing`, which absorbed the former `build_phasing_graph` crate) and
consensus / similarity / BILC solve (`haplotype_inspection`). Both are called
in-process by `fp-control`; `phasing` also retains a standalone phaser binary
and differential `examples/` harnesses.

---

### Orchestrator

#### `sdrecall` — top-level binary (`run` / `prepare` / `realign`)
The Rust equivalent of the Python `SDrecall` executable; argument names, short
flags and defaults match 1:1.
```bash
# full pipeline (preparation → realign + recall → post-merge)
cargo run -p sdrecall --release -- run \
  -i <input.bam> -r <hg38.fasta> -m <sd_map.tsv> -b <target.bed> -o <outdir> \
  -s <sample_id> -t 16 --target_tag exome \
  --conventional_vcf <deepvariant.vcf.gz> --caller_name DeepVariant \
  --cohort_vcf <control_cohort.vcf.gz>

cargo run -p sdrecall --release -- prepare -i ... -r ... -m ... -b ... -o ...  # Phase 1 only
cargo run -p sdrecall --release -- realign -i ... -r ... -m ... -b ... -o ...  # Phase 2 only
```
- **Inputs (common):** `-i` indexed BAM, `-r` reference `.fasta` (suffix
  enforced), `-m` reference SD map, `-b` target BED, `-o` output dir; `-s`
  sample id (else derived from the BAM name), `-t` threads (default 10),
  `--target_tag` (default `exome`), `--mq_cutoff` (default 41). `run`/`realign`
  also accept the optional `--conventional_vcf` / `--cohort_vcf` groups.
- **Outputs:** the SDrecall VCF under
  `<outdir>/<sample>_<assembly>_<tag>_SDrecall/recall_results/`, and — when the
  optional VCFs are supplied — the inhouse-common-annotated and
  conventional-merged VCFs. The final path is printed on success.

> **Status:** the Rust orchestrator has completed production end-to-end runs on
> t2t/chm13, hg19, and hg38, with retained-output parity checks. The remaining
> formal migration gates are the HG006 Rust-versus-Python differential and
> wiring the Python pipeline's NM Poisson cutoff into per-island inspection.

---

## 5. Validation harnesses (`examples/`)

Differential/validation harnesses link a finished library and live in each
crate's `examples/`. Build/run with `cargo run -p <crate> --example <name> --
<args>`:

| Crate | Example | What it checks |
|-------|---------|----------------|
| `fp-control` | `diff_vs_hybrid` | fused Rust path vs the Rust/Python hybrid (per-island parity) |
| `vcf-ops` | `diff_vcf_ops` | merge / inhouse-common vs Python on real VCFs |
| `region-prep` | `diff_region_prep` | fc/nfc BEDs byte-identical vs Python |
| `haplotype_inspection` | `test_bam_lapper`, `test_consensus`, `test_phase1`, `validate_bilc_solver`, `validate_select_regions` | submodule oracles |

Per the project test convention, run harnesses with `RUST_LOG=debug` and tee the
output to `/paedyl01/disk1/yangyxt/test_tmp/<name>.log`.

---

## 6. Where to look next

- Workspace members & lints: [`Cargo.toml`](Cargo.toml)
- Orchestration spine: [`sdrecall/src/pipeline.rs`](sdrecall/src/pipeline.rs)
- Migration status, per-task contracts, validation evidence:
  [`../docs/analysis/RUST_MIGRATION_SUMMARY.md`](../docs/analysis/RUST_MIGRATION_SUMMARY.md)
  and [`../docs/analysis/tasks/`](../docs/analysis/tasks/)
- Code-review tracker:
  [`../docs/analysis/REVIEW_FINDINGS.md`](../docs/analysis/REVIEW_FINDINGS.md)
