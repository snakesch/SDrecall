# T5 — `vcf-ops` crate (sorted-VCF merge + inhouse-common + HP annotate)

**Crate:** `vcf-ops` (new — lib + bin)
**Status (2026-07-17):** Production-integrated. Priority merge and inhouse-common logic run in Rust; `bcftools` normalization remains an allowed external dependency. The next accuracy track is canonical-caller merge followed by HG002 truth benchmarking.
**Depends on:** T0
**Replaces (Python):** `src/merge_variants_with_priority.py` (`merge_with_priority`), `identify_common_vars.py` (`annotate_inhouse_common`), and — optional sub-scope — `realign_recall/annotate_HP_tag_to_vars.py`.

## Goal & scope boundary

The two big Python post-processing files share an **identical sorted two-pointer VCF co-iteration skeleton** (`compare_variant_positions` / `deal_with_same_loc_variants`), differing only in the per-locus callback. Migrate the engine **once** (per the #1 coding rule) and parameterize the callback:

- **priority merge** (`merge_with_priority`): merge GT from `AD`/`GQ`/`HPSUP` FORMAT fields, apply source/priority tags.
- **inhouse-common** (`annotate_inhouse_common`): binomial test `binom.cdf(AC, AN, cutoff) > conf_level` against a cohort VCF → add `INHOUSE_COMMON` filter.

Optional third callback: HP-tag annotation from BAM pileup (`annotate_HP_tag_to_vars`) — a BAM→VCF per-variant pileup op; include here or defer to T9.

Out of scope: nothing shells out by default — sort/concat are in-process via rust-htslib. **Indel left-normalization (`bcftools norm`) is the one genuine gap** (needs ref-FASTA + left-align): reimplement on rust-htslib, or keep as an isolated leaf subprocess — see the plan's External tools & libraries policy.

## Data flow

```
query VCF (sorted) ─┐
                    ├─ 2-pointer co-iteration (sequential contig pass)
cohort/ref VCF ─────┘     per matched locus → callback:
                              merge-priority  → merged VCF
                              inhouse-common  → annotated VCF
```

Both Python paths ran per-chromosome under `multiprocessing.Pool`. Rust currently uses a deterministic sequential contig pass because `bcf::Record` carries a non-`Send` header reference; safe contig parallelism remains optional.

## Dependencies

- Crates: `rust-htslib` 0.47.0 (VCF/BCF), `statrs` 0.18.0 (binomial), and `sdrecall-utils`/`sdrecall-io`.
- I/O: `rust-htslib` (VCF read/write/sort/concat, in-process). Only indel left-normalization may need a `bcftools norm` leaf subprocess (see policy).

## Performance bottleneck / rationale

The per-record pysam (de)serialization across the `multiprocessing` pickling boundary dominates today. Rust eliminates the pickling and the interpreter overhead; the binomial/ratio math is trivial.

## Tests

### Unit (tier 1)
- 2-pointer merge cases: query-only, ref-only, same-locus, multiallelic, overlapping indels.
- Binomial cutoff at AC/AN boundaries.
- GT-correction rules from `AD`/`GQ`/`HPSUP` thresholds (force `GT=(1,1)`).

### Differential vs Python (tier 2)
- Run both callbacks on HG006: Rust output VCF vs Python output VCF.

**Pass criterion:** record-identical output VCF (after `bcftools norm` normalization) vs Python on HG006, for both the priority-merge and inhouse-common paths.

**Data:** HG006 query VCF + a cohort VCF.

## Progress
- [x] Scaffold `vcf-ops` crate (lib + bin with `merge` / `inhouse-common` subcommands) — `rust_modules/vcf-ops/`
- [x] Shared 2-pointer co-iteration engine (callback-parameterized) — `coiterate.rs::coiterate_sorted_vcfs` + `LocusOp` trait (the single unit; both Python copies collapsed)
- [x] priority-merge callback + GT correction — `priority_merge.rs` + `gt_rules.rs::should_force_hom` (the **4** Python GT ladders collapsed into one evaluator over `static` rule tables; exact FILTER-tag END-append ordering preserved)
- [x] inhouse-common binomial callback — `inhouse_common.rs::determine_common` (`statrs::Binomial::new(cutoff, AN).cdf(AC) > conf_level`, correct (p,n) arg order, cross-checked vs scipy ~1e-7)
- [~] rayon over chromosomes — **deferred (documented):** `bcf::Record` holds `Rc<HeaderView>` (non-atomic refcount) → sequential per-contig pass. The headline win (eliminating the `multiprocessing` pickling boundary) is already delivered; contig-axis parallelism is a future safe-refactor (per-worker header-cloned reader).
- [x] Differential vs Python — **PASS on real HG002 data** (HG006 still the formal target)
- [ ] (optional) HP-tag-from-pileup callback — deferred to T9 (BAM-pileup op, needs sdrecall-io pileup helper)
- [x] Merge/benchmark workflow implemented in `benchmarks/benchmark_hg002_hg38.sh` and completed on PBS job `72607` (exit 0).

### Benchmark handoff (2026-07-17)

- Rust recall VCF: `/paedyl01/disk1/yangyxt/SDrecall-test/results/sdrecall_rust_external_helper_limits_20260717/hg38/HG002_hg38_exome_avg50x_SDrecall/recall_results/HG002.sdrecall.vcf.gz`
- Canonical caller candidate: `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002.deepvariant.sorted.vcf.gz`
- GIAB truth VCF: `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.vcf.gz`
- GIAB confident BED: `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.bed`
- Reference: `/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta`
- SD target BED: `/paedyl01/disk1/yangyxt/SDrecall-test/results/sdrecall_rust_external_helper_limits_20260717/hg38/HG002_hg38_exome_avg50x_SDrecall/realign_groups/all_target_recall_SD_regions.bed`

`benchmarks/benchmark_hg002_hg38.sh` now evaluates DeepVariant alone and DeepVariant+accepted SDrecall calls over the same SD-target ∩ GIAB-confident BED. Its bundled comparator is an exact normalized allele/dosage regression screen, not haplotype-aware GIAB reconciliation; publication-grade claims still require hap.py or RTG vcfeval. Start with hg38 because the retained t2t run and available canonical/truth resources use incompatible CHM13 versions/naming.

### Exact-screen result (2026-07-17)

PBS job `72607` completed in 1m31s with exit status 0 and about 3.9 GiB reported memory. The callable scope is 3,941 intervals / 579,383 bases. `MISALIGNED` SDrecall records were excluded; the normalized merged set contains every one of the 416 DeepVariant calls plus 1,137 accepted SDrecall additions.

| Callset | Genotype TP/FP/FN | Genotype P/R/F1 | Exact-allele TP/FP/FN | Exact-allele P/R/F1 |
|---|---:|---:|---:|---:|
| DeepVariant | 376 / 40 / 356 | 0.903846 / 0.513661 / 0.655052 | 403 / 13 / 329 | 0.968750 / 0.550546 / 0.702091 |
| DeepVariant + SDrecall | 607 / 946 / 125 | 0.390856 / 0.829235 / 0.531291 | 647 / 906 / 85 | 0.416613 / 0.883880 / 0.566302 |

The accepted SDrecall additions recover 244 exact truth alleles (231 with matching alternate dosage) and add 893 exact-allele false positives (906 genotype-level false positives). Thus recall rises by 0.333333 at the allele level and 0.315574 at the genotype level, while exact-screen F1 falls by 0.135789 and 0.123761 respectively. This is evidence of a sensitivity/precision trade-off, not a definitive GIAB score: haplotype-aware reconciliation may rescue some complex-representation mismatches, but the false-positive burden still requires stratification and review.

Artifacts: `/paedyl01/disk1/yangyxt/SDrecall-test/benchmarks/hg002_hg38_sdrecall_deepvariant_20260717/metrics/HG002.callset_comparison.tsv` and `/paedyl01/disk1/yangyxt/SDrecall-test/benchmarks/hg002_hg38_sdrecall_deepvariant_20260717/logs/benchmark.log`.

### Result (2026-06-12) — 32 unit tests + real-data differential PASS

- **Unit:** 32 tests pass, clippy clean. Logs `/paedyl01/disk1/yangyxt/test_tmp/vcf_ops_test_final_20260612.log`, `…/vcf_ops_clippy_final_20260612.log`.
- **Differential vs Python on real HG002:**
  - **priority-merge:** query `HG002.sdrecall.raw.vcf.gz` vs ref `HG002.sdrecall.clean.vcf.gz` (ucsc.hg38) → **rust 1814 == py 1814, diffs = 0** (record-identical on CHROM/POS/REF/ALT/FILTER/GT/AD/GQ after `bcftools norm`). Log `…/diff_vcf_ops_merge.log`.
  - **inhouse-common:** real data rust 415 == py 415 (no `SDrecall` FILTER present → tagging branch not exercised); a real-ref-validated 3-record **synthetic** fixture exercises positive tagging → Rust and Python **byte-identical** incl. `SDrecall;INHOUSE_COMMON` ordering.
- **bcftools-norm:** isolated leaf subprocess (`norm.rs`, `set -o pipefail`, hard error on non-zero), per the decided external-algorithm policy.
- **#1 rule honored:** one co-iteration engine + one GT evaluator; the two contig-filter constants (`is_main_contig` / `is_inhouse_contig`) kept separate as the Python requires.

## Migration design — frontier propagation (2026-06-11)

> Coding-grade design for the `vcf-ops` crate. Builds on the data-flow + test plan above. The
> central decision is realized: **one** generic sorted two-pointer co-iteration engine, parameterized
> by a per-locus callback + per-set finalizers (the #1 coding rule). All VCF I/O is in-process via
> `rust-htslib` 0.47.0; the only deliberately-external piece is indel left-normalization (see Risks).

### 1. Python logic inventory (functions to port, with file:line)

Both production files carry a **verbatim-duplicated** co-iteration skeleton. The engine is ported
once; the differences are isolated to (a) the per-locus match callback and (b) what happens to the
three result sets (matched / query-only / ref-only).

**`src/merge_variants_with_priority.py`** (the priority-merge path):
- `compare_variant_positions(v1, v2)` — L148-165. Total order on `(contig, start, stop)`; returns
  `-1` (different contig → error inside a region), `0` (same start AND same stop), `1` (v1 before v2),
  `2` (v1 after v2). Note: ordering is by **0-based `start`/`stop`**, NOT by alleles.
- `VariantRecordWrapper.__eq__` — L34-39. Two records are EQUAL only if `chrom == pos == ref == alts`
  (allele-level), `__hash__` L31-32 on `(chrom, pos, tuple(alleles))`. So "same location" (order==0)
  is a coarser bucket than "same variant" (eq). This is the crux of `deal_with_same_loc_variants`.
- `deal_with_same_loc_variants(...)` — L168-238. When two records share a *location*, drains both
  iterators of all co-located records into `current_loc_1recs` / `current_loc_2recs` (pushing the
  first downstream record into `buffer_var{1,2}`), then does an O(n·m) allele-equality match between
  the two co-located lists; matched pairs → callback → `merged_records`; unmatched v1 → `non_overlap_var1`;
  unmatched v2 → `non_overlap_var2`. Single shared skeleton with `identify_common_vars.py`.
- `process_region(region, query_vcf, reference_vcf, ...)` — L246-394. The outer two-pointer loop:
  fetch both per-contig iterators, advance the lagging pointer by `compare_variant_positions`, call
  `deal_with_same_loc_variants` on order==0, accumulate three frozensets, return `.pickable()` tuples.
  Empty-iterator edge cases L268-272.
- `modify_gt_based_on_ad_gq(record, rrecord)` — L62-112. The **matched-pair** GT-correction callback
  body. Per sample: reads ref/alt from `AD`, `GQ`, and `HPSUP` (`num_hps = len(hps[0].split(";"))`).
  Forces `GT=(1,1)` when (in order, first hit wins): `ralt/rdp ≥ 0.9`; `num_hps≥2 & ratio≥0.33`;
  `num_hps≥3 & ratio≥0.30`; `num_hps≥4 & ratio≥0.25`; `rgq<5 & ratio≥0.5 & rdp>5`. Skips if GT already (1,1).
- `merge_same_variant_rec(qrecord, rrecord, qv_tag, rv_tag, modify_gt)` — L115-135. The locus callback:
  add `qv_tag`/`rv_tag` to the REFERENCE record's FILTER, then GT-correct it. Returns the merged
  (reference-derived) record.
- `merge_vcf_headers(header1, header2)` — L433-463. Union of FILTER/INFO/FORMAT/contig/other header lines.
- The three **per-set writer finalizers** in `merge_with_priority` — L536-657. Each set gets a DIFFERENT
  GT-correction at write time: matched (L536-556, write as-is + reorder tags); query-only (L557-612,
  thresholds: `num_hps≥2 & ratio≥0.55 & dp≥5`; `num_hps≥4 & alt≥ref & dp≥5`; `ratio≥0.9 & dp≥5`; plus
  `HPSUP` list→`;`-joined string); ref-only (L613-657, threshold: `gq<5 & alt/total≥0.7`). **These are not
  one callback** — they are three post-collection per-set transforms keyed by which set the record landed in.
- `merge_with_priority(...)` — L467-676. Orchestrator: `sort_vcf` both inputs (L488/491), main-contig
  filter `main_contigs` L513 (`chr1..22,X,Y,M` + no-`chr` aliases), `multiprocessing.Pool` over contigs
  L520-529, write merged header, final `bcftools norm -d exact | bcftools sort` L671.
- `sort_vcf(vcf_file, ref_genome, ...)` — L398-429. `bcftools norm -m -both -f REF --multi-overlaps 0 -a`
  (split multiallelics + **left-align**) → `norm -d exact` (dedup) → `filter -i 'ALT[0]!="*" && COUNT(GT="alt")>0'`
  → `sort`. This is the only step that needs a reference FASTA.

**`identify_common_vars.py`** (the inhouse-common path — IDENTICAL skeleton):
- `compare_variant_positions` L124-141, `VariantRecordWrapper` L62-96, `deal_with_same_loc_variants`
  L144-214, `process_region` L222-370 — **byte-for-byte the same control flow** as the merge file
  (only variable names `reference`→`cohort` differ). This is the duplication the engine collapses.
- `determine_common_per_pysam_record(record, AC_tag, AN_tag, inhouse_common_cutoff, conf_level)` — L18-59.
  The locus callback math: read `AC`/`AN` from INFO (AC list→`AC[0]`), `na_value` guard, then
  `stat_power = binom.cdf(AC, AN, inhouse_common_cutoff)`; common iff `stat_power > conf_level`.
- `process_target_records(qrecord, crecord, added_filter="INHOUSE_COMMON", ...)` — L100-111. The locus
  callback: if `is_common AND "SDrecall" in qrecord.filter` → add `INHOUSE_COMMON` to the QUERY record's
  FILTER. Returns the query record (cohort record is only a data source).
- `annotate_inhouse_common(...)` — L373-486. Orchestrator. **Crucially writes only merged + query-only**
  (cohort-only records are *dropped*, L426-469 has no `non_overlap_crecs` writer). Contig filter is the
  regex `^chr[0-9MTXY]+$` L418 (slightly different from merge's `main_contigs` set — see Risks).

**`realign_recall/annotate_HP_tag_to_vars.py`** (optional 3rd callback — a BAM→VCF pileup op, NOT a co-iteration):
- `get_supporting_tags(bam, chrom, pos, ref, alts, tag, min_mapq, min_bq)` — L7-45. Per-variant pileup
  (`samtools pileup` truncate), classify each read as del/ins/SNV, collect the read's `HP` tag value into
  `supporting_tags[allele]`.
- `annotate_vcf(input_vcf, output_vcf, bam, tag, ...)` — L49-76. Per-record loop adding a `{tag}SUP`
  FORMAT field. **Different shape** (BAM pileup, not VCF×VCF) — recommend defer to T9 (it depends on the
  `sdrecall-io` pileup helpers, not yet specified). Listed here for completeness only.

Control-flow shape: outer per-contig parallelism (`multiprocessing.Pool` → **rayon `par_iter` over contigs**);
inside a contig, a strictly sequential single-pass two-pointer merge (NOT parallelizable — order-dependent
buffers). Each contig is independent → embarrassingly parallel across the contig axis.

### 2. Python → Rust crate mapping

| Python operation / idiom | Rust crate::api | confidence |
|---|---|---|
| `pysam.VariantFile(path)` (read) | `rust_htslib::bcf::Reader::from_path` / `IndexedReader::from_path` (via `sdrecall_io::read_vcf`) | verified-docs |
| `vcf.fetch(region)` per-contig | `bcf::IndexedReader::fetch(&mut self, rid: u32, start: u64, end: Option<u64>) -> Result<()>` then `.records()`; `rid` via `header().name2rid(contig.as_bytes())` | verified-docs |
| `pysam.VariantFile(path,'w',header)` + `.write(rec)` | `bcf::Writer` (via `sdrecall_io::write_vcf(path, &header, recs)`) | verified-docs |
| `compare_variant_positions` (start/stop order) | inline `cmp` on `(rec.rid(), rec.pos(), rec.end())`; `pos()->i64`, `end()->i64` | verified-docs |
| `VariantRecordWrapper.__eq__` (chrom,pos,ref,alts) | `LocusKey { rid:u32, pos:i64, alleles:Vec<Box<[u8]>> }` from `rec.alleles() -> Vec<&[u8]>` | verified-docs |
| `record.samples[s]['AD']` (FORMAT int) | `rec.format(b"AD").integer() -> Result<BufferBacked<Vec<&[i32]>>>` (one `&[i32]` per sample) | verified-docs |
| `record.samples[s]['GQ']` | `rec.format(b"GQ").integer()` (per-sample `&[i32]`) | verified-docs |
| `record.samples[s]['HPSUP']` (string, `;`-split) | `rec.format(b"HPSUP").string() -> Result<BufferBacked<Vec<&[u8]>>>`; split on `b';'` | verified-docs |
| `record.samples[s]['GT']` read | `rec.genotypes()?.get(s) -> Genotype` (Vec<GenotypeAllele>) | verified-docs |
| set `GT=(1,1)` | `rec.push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)])` | plausible |
| `record.filter.add(tag)` | `rec.push_filter(tag_id)` where `tag_id: &str`/id impl `FilterId`; reorder via `set_filters(&[&str])` | verified-docs |
| `tag in record.filter` / `"SDrecall" in filter` | `rec.has_filter(b"SDrecall") -> bool` | verified-docs |
| `record.info.get('AC'/'AN')` | `rec.info(b"AC").integer() -> Result<Option<BufferBacked<&[i32]>>>` | verified-docs |
| `binom.cdf(AC, AN, cutoff)` (scipy `cdf(k,n,p)`) | `statrs::distribution::Binomial::new(cutoff /*p*/, AN /*n:u64*/)?.cdf(AC /*x:u64*/)` via `DiscreteCDF` | verified-docs |
| `na_value(x)` | `sdrecall_utils::na_value(&str)`; for htslib missing use `bcf::record::*::is_missing()` / sentinel check | plausible |
| `multiprocessing.Pool(threads).imap_unordered(per_contig)` | `contigs.par_iter().map(process_contig).collect()` (`rayon` 1.8) | verified-docs |
| `merge_vcf_headers` (union FILTER/INFO/FORMAT/contig) | iterate `HeaderView` records, `Header::push_record` for missing ids on the writer header | plausible |
| `sort_vcf`: `bcftools norm -m -both -f REF` (split + **left-align**) | **decision** — leaf `bcftools norm` subprocess OR rust-htslib + `bio`/manual left-align (see Risks) | uncertain |
| `bcftools sort` / `norm -d exact` (final) | `sdrecall_io::concat_sort_vcfs(inputs, out, dedup_exact=true, threads)` (in-process) | plausible |
| `frozenset(...).pickable()` cross-process transport | **eliminated** — rayon shares memory; per-contig returns owned `Vec<bcf::Record>` directly, no (de)serialization | verified-docs |

### 3. Crate file layout (`rust_modules/vcf-ops/`)

`lib + bin`, matching the `phasing` crate style (thin `lib.rs`, `main.rs` CLI, one module per job):

```
vcf-ops/
├─ Cargo.toml                 # [lib] name="vcf_ops"; [[bin]] name="vcf-ops"; deps: rust-htslib, statrs, rayon, sdrecall-utils, sdrecall-io, clap, log
└─ src/
   ├─ lib.rs                  # pub use; re-exports coiterate + the two public entry fns
   ├─ coiterate.rs            # THE one engine (see §4). Generic over a LocusOp callback + SetSink.
   ├─ priority_merge.rs       # merge-priority LocusOp + the 3 per-set GT-correction finalizers
   ├─ inhouse_common.rs       # inhouse-common LocusOp (binomial) + SDrecall-filter guard
   ├─ gt_rules.rs             # ONE versatile GT-correction unit (see §4) — shared AD/GQ/HPSUP threshold evaluator
   └─ main.rs                 # clap: `merge` / `inhouse-common` subcommands (files in → files out)
examples/
   └─ diff_vcf_ops.rs         # differential harness: run a path, compare to a dumped Python VCF
```

**One-versatile-unit-per-job mapping (collapsing Python duplication):**
- `coiterate.rs::coiterate_sorted_vcfs` — the SINGLE engine replacing the two byte-identical
  `process_region` + `deal_with_same_loc_variants` copies. Differs only by the injected `LocusOp`.
- `gt_rules.rs::apply_gt_correction` — the SINGLE GT-rule evaluator. The Python has FOUR separate
  threshold ladders (`modify_gt_based_on_ad_gq` + three per-set finalizers) that all read the same
  `(ref_dp, alt_dp, gq, num_hps)` quadruple and decide `force_hom`. Collapse into one unit that takes
  a `&[GtRule]` rule-table (each `GtRule` = a predicate over the quadruple); the four call-sites pass
  four different rule tables. This kills the largest duplication in the Python without broadening scope.
- `priority_merge.rs` / `inhouse_common.rs` are thin: each provides a `LocusOp` impl + (merge only) the
  set-routing rule tables; they are genuine downstream orchestration, not overlapping helpers.

### 4. Core data structures + key fn signatures (explicit borrow/owner choices)

```rust
// --- the one engine ------------------------------------------------------
/// Result of one matched-or-unmatched locus, owned because records cross the
/// rayon boundary and are mutated independently per set.
pub struct CoiterSets { pub matched: Vec<bcf::Record>, pub query_only: Vec<bcf::Record>, pub ref_only: Vec<bcf::Record> }

/// The per-locus callback. `&mut q` because the merge path mutates the kept record's
/// FILTER/GT in place; returns which set the (possibly-mutated) record(s) land in.
pub trait LocusOp {
    /// Called once per (query, ref) pair that share BOTH location AND alleles.
    /// May mutate either record; returns the record to keep + its destination set.
    fn on_match(&self, q: &mut bcf::Record, r: &mut bcf::Record) -> Kept;
}
pub enum Kept { Matched(bcf::Record), DropMatch } // Matched carries the chosen record

/// THE engine. Borrows readers mutably (htslib fetch needs &mut), borrows the op,
/// owns nothing it returns by reference. One contig per call → rayon-parallel.
pub fn coiterate_sorted_vcfs(
    query: &mut bcf::IndexedReader,   // &mut: fetch advances internal state
    refr:  &mut bcf::IndexedReader,
    contig: &str,
    op: &dyn LocusOp,                 // &dyn: monomorphization not needed, one call per contig
) -> Result<CoiterSets>;

// --- the one GT-rule unit ------------------------------------------------
/// A single threshold predicate over the AD/GQ/HPSUP quadruple. `Copy`, cheap.
#[derive(Clone, Copy)]
pub struct GtRule { pub min_num_hps: u32, pub min_ratio: f64, pub min_dp: u32, pub max_gq: Option<i32>, pub alt_ge_ref: bool }
/// Evaluate a rule table against per-sample stats; `&` everywhere — read-only, no alloc.
/// Returns true on the first matching rule (Python's first-hit-wins ladder).
pub fn should_force_hom(stats: &SampleStats, rules: &[GtRule]) -> bool;
/// Cheaply extracted per sample; borrows the buffer-backed htslib slices.
pub struct SampleStats { pub ref_dp: i32, pub alt_dp: i32, pub gq: i32, pub num_hps: usize }
pub fn sample_stats(rec: &bcf::Record, sample: usize) -> SampleStats; // & — read-only view

// --- public entry points (files in → files out) -------------------------
pub fn merge_with_priority(p: MergeParams<'_>) -> Result<()>;        // borrows paths via &Path in params
pub fn annotate_inhouse_common(p: InhouseParams<'_>) -> Result<()>;
```

WHY each choice:
- `CoiterSets` returns **owned** `Vec<bcf::Record>` — records are produced inside a rayon worker and
  later concatenated/written on the main thread; owning them is the natural move (no lifetime tangle,
  no `Arc`). This replaces the Python `.pickable()` tuple round-trip entirely (zero serialization).
- `coiterate_sorted_vcfs` takes `&mut IndexedReader` because htslib `fetch` mutates reader state; `&dyn
  LocusOp` (not generic `<O: LocusOp>`) keeps `coiterate.rs` non-generic and compile-fast — there are
  exactly two ops, dynamic dispatch cost is negligible vs the I/O.
- `GtRule: Copy` + `&[GtRule]` — the rule tables are tiny `static` arrays; passing `&[GtRule]` is
  zero-copy and lets all four Python ladders share `should_force_hom`.
- `sample_stats(&rec)` borrows; `SampleStats` holds plain `i32`/`usize` copied out of the buffer-backed
  slice so the htslib `Buffer` borrow ends immediately (avoids holding `format()`'s `&'a self` across the
  GT mutation, which would conflict with `&mut rec`).
- `&Path` in `*Params<'_>` — paths are caller-owned and read-only.

### 5. Performance optimizations from ownership/borrowing

- **Kill the pickling boundary**: Python ships every record across the `multiprocessing` pickle as a
  `.pickable()` tuple (the dominant cost per the bottleneck note). Rayon shares heap memory → records
  move as owned `bcf::Record`, never serialized. This is the headline win.
- **Zero-copy FORMAT/INFO reads**: `rec.format(b"AD").integer()` returns a `BufferBacked<Vec<&[i32]>>`
  borrowing htslib's internal buffer — read `ref_dp`/`alt_dp` directly off the slice, copy only the two
  `i32` into `SampleStats`, drop the borrow before mutating. No per-record `to_dict()` clone (Python L24).
- **Reuse one `Buffer`**: use `format_shared_buffer`/`info_shared_buffer` with a per-thread reusable
  `Buffer` to avoid a fresh allocation per field per record.
- **rayon over the independent axis**: `contigs.par_iter()` — contigs are independent; inside a contig the
  merge is sequential (order-dependent buffers, must stay single-threaded). Matches the Python contig pool.
- **FxHashMap / ahash for keys**: contig→rid lookups and the co-located allele-match use `rustc_hash::FxHashMap`
  for the `u32` rid keys; the small per-locus allele match (usually 1×1) stays a linear scan (faster than a map
  for n,m≤2, per the perf guidance "HashMap for small sets → overhead").
- **`with_capacity`** on the three per-contig `Vec`s sized from the contig's record count estimate.
- **Streaming write**: collect per-contig `CoiterSets`, then stream-write in contig order to one `bcf::Writer`
  — no global materialization of all records (only one contig's worth resident at a time per worker).

### 6. Risks / open decisions

- **Indel left-normalization (the one genuine external dependency).** Python's `sort_vcf` runs
  `bcftools norm -m -both -f REF --multi-overlaps 0 -a` (split multiallelics + left-align against the ref
  FASTA) BEFORE the co-iteration, and the differential pass criterion is "record-identical *after*
  `bcftools norm`". Options: **(A)** isolated leaf `bcftools norm` subprocess (deliberate
  external-ALGORITHM choice per the plan's External-tools policy — flag, not a missing-dep fallback);
  **(B)** reimplement split+left-align on `rust-htslib` + a `bio`/`faidx` ref reader. **Recommendation:
  (A) for T5** (parity-critical, bcftools' left-align edge cases are subtle) and revisit (B) at T9 when
  the FASTA reader exists. **Needs user's call.** This is the single biggest parity hazard.
- **`set_filters` ordering parity.** Python L538-541/L559-562 explicitly *reorders* FILTER tags
  (`[f for f in filters if f != tag] + [tag]` — move the source tag to the end). htslib `set_filters`
  must replicate this exact order or the differential VCF diff (string-level) fails even when semantically
  equal. Mitigation: build the ordered `Vec<&str>` then one `set_filters`.
- **Contig-filter divergence between the two paths.** Merge uses the `main_contigs` set (L513, includes
  `chrM`/`MT`/no-`chr` aliases); inhouse uses regex `^chr[0-9MTXY]+$` (L418, `chr`-prefixed only). Keep
  them as two separate constants in the respective op modules — do NOT unify (would change behavior).
- **`HPSUP` get-default + missing-field semantics.** Python uses `.get('HPSUP', '.')` then
  `len(hps[0].split(';'))`; missing AD/GQ default to `[0,0]`/`0`. htslib returns `Err`/missing-sentinel,
  not a Python default — `sample_stats` must map missing → the same defaults exactly (off-by-one in
  `num_hps` flips a GT call). The HPSUP value is a `;`-joined string written back as one token on the
  query-only path (L596) — preserve that join.
- **`Genotypes::get` is diploid-only** (verified). All SDrecall samples are single-sample diploid, so OK;
  assert `sample_count()==1` and fail clearly otherwise (dependency-availability spirit: no multi-sample
  fallback path).
- **statrs vs scipy binomial CDF.** Both are exact (regularized incomplete beta); cross-check a handful of
  `(AC,AN,cutoff)` boundary points against scipy in the differential test. `Binomial::new` rejects p∉[0,1]
  / returns `BinomialError` — propagate as `SdError::Vcf`.
- **`AC` is a list in INFO** (Python takes `AC[0]`, L46) — read `info(b"AC").integer()` and index `[0]`.
- **Version skew**: rust-htslib 0.47 pinned in `[workspace.dependencies]`; statrs 0.18 is new to the
  workspace — add to `[workspace.dependencies]` at scaffold time (T5 coding, not this design pass).
- **HP-tag callback deferred**: `annotate_HP_tag_to_vars` is a BAM-pileup op, not VCF×VCF; it needs
  `sdrecall-io` pileup helpers that aren't in the T0 contract. Defer to T9 unless the user wants it folded
  in (would add a `bcf::Record`-from-pileup unit, out of the co-iteration engine's scope).

### 7. Test plan delta (concrete fixtures + assertions)

**Unit (tier 1), in `#[cfg(test)]` with hand-built `bcf::Record`s on a tiny in-memory header:**
- `coiterate_sorted_vcfs` cases: query-only run; ref-only run; single same-locus same-allele match;
  same-locus DIFFERENT alleles (→ both go to their *_only sets, NOT matched); multiple co-located records
  at one position (drains both, O(n·m) allele match, buffer carries the first downstream record);
  overlapping indels with differing `stop` (order by stop). Assert exact `(matched, query_only, ref_only)`
  membership by `LocusKey`.
- `should_force_hom`: table-drive each Python ladder. Priority-matched ladder boundary points:
  `ralt/rdp` = 0.89/0.90 (≥0.9 rule), `num_hps=2 & ratio` 0.32/0.33, `=3 & 0.29/0.30`, `=4 & 0.24/0.25`,
  `rgq=4/5 & ratio 0.5 & rdp 6`. Query-only ladder: `num_hps=2 & ratio 0.54/0.55 & dp 4/5`; `num_hps=4 &
  alt<ref vs alt≥ref & dp≥5`; `ratio 0.89/0.90 & dp≥5`. Ref-only ladder: `gq 4/5 & alt/total 0.69/0.70`.
  Assert `force_hom` boolean matches Python at each boundary.
- `inhouse_common` binomial: `Binomial::new(cutoff, AN).cdf(AC) > conf_level` at AC/AN boundaries
  (e.g. cutoff=0.01, conf=0.999: find the AC where cdf crosses 0.999 for AN=100, assert flip); plus the
  `is_common AND has_filter("SDrecall")` guard (common-but-no-SDrecall → no INHOUSE_COMMON).
- FILTER reorder: assert the written FILTER list equals Python's `[...others] + [tag]` order.

**Differential (tier 2), on HG006 (`examples/diff_vcf_ops.rs`):**
- Dump Python outputs: run `merge_with_priority` and `annotate_inhouse_common` on the HG006 query VCF +
  a cohort VCF (data per the doc), keep the Python output VCFs as golden.
- Run the Rust `merge` / `inhouse-common` subcommands on the SAME sorted inputs (apply the same
  pre-`norm` so the only variable is the engine, isolating the normalization decision from Risks).
- **Assertion**: after a common `bcftools norm` pass on both, `bcftools view` body is byte-identical
  (sort by `(chrom,pos,ref,alt)`; compare CHROM/POS/REF/ALT/FILTER/INFO + the GT/AD/GQ FORMAT fields).
  Report the first differing record. Pass = zero differing records on both paths.
- Log at `RUST_LOG=debug` to `/paedyl01/disk1/yangyxt/test_tmp/diff_vcf_ops_<path>.log` per the project's
  test-verification requirement.
