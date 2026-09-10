# Code Review Findings — Migrated Rust Crates

**Date:** 2026-06-11
**Scope:** the three already-migrated crates — `read_extraction`, `build_phasing_graph`, `haplotype_inspection` (production `src/` only; example/test harnesses excluded).
**Method:** four parallel review passes (code reuse, code quality/idiom, efficiency/hot-path, FFI/unsafe-safety), then the headline and subtle items were spot-verified by reading the actual code.

This file is the **single tracker** for issues found in that review. Each task doc (`tasks/T*.md`) carries a short "Review findings" section that points back here for the items it owns.

---

## How to read this

- **Severity** — practical impact on this pipeline (not generic severity). `HIGH` = worth real effort; `MED` = should fix but bounded; `LOW` = cosmetic/cleanup.
- **Owner** — the migration task where the fix naturally belongs, so we don't churn validated code out of band.
- **Verified** — ✅ = I read the code and confirmed it personally; 🔎 = reported by a review pass with `file:line`, high-confidence but not line-by-line re-confirmed here.
- **Status** — `open`, `fixed (2026-06-11)`, or `verify` (needs a check against the Python original before any change).

---

## Summary table

| ID | Finding | Sev | Category | Where | Owner | Verified | Status |
|----|---------|-----|----------|-------|-------|----------|--------|
| DUP-1 | BAM read→filter→interval implemented 3× | HIGH | Duplication | `bam_reading.rs`, `bam_lappers.rs`, `read_extraction/lib.rs` | T0 (sdrecall-io) | 🔎 | open |
| DUP-2 | Hap/error-vector + base-encoding toolkit duplicated across 2 crates | HIGH | Duplication | `haplotype_determination.rs` ↔ `pairwise_read_inspection.rs` | T0 (sdrecall-utils) | 🔎 | open |
| DUP-2a | The two `extract_hap_vector` copies have **diverged** in 2 ways (M-op AND insertion-marker placement) | HIGH | Correctness risk | `pairwise_read_inspection.rs:83`, `haplotype_determination.rs:147` | T0 + T4 | ✅ | open |
| DUP-3 | Interval-overlap logic implemented 3 ways | MED | Duplication | `structs.rs:115`, `bam_lappers.rs:663`, `graph_builder.rs` | T0 | 🔎 | open |
| DUP-4 | Shared-SNV / position→base mapping duplicated | MED | Duplication | `haplotype_determination.rs:1090` ↔ `identify_misaligned_haps.rs:389` | T0 | 🔎 | open |
| DUP-5 | Minor interval-merge / sort-dedup helpers overlap | LOW | Duplication | `identify_misaligned_haps.rs:270`, `graph_builder.rs:504` | T0 | 🔎 | open |
| PERF-1 | Cached vectors cloned on every cache hit (hot path) | HIGH | Performance | `identify_misaligned_haps.rs:607-710` | T4 | ✅ | open |
| PERF-2 | No parallelism — per-haplotype/per-region loops are single-threaded | HIGH | Performance | `identify_misaligned_haps.rs:2000,2153` | T4 | 🔎 | open |
| PERF-3 | Default (slow-hash) maps/sets on hot paths instead of fast-hash | MED | Performance | `bam_lappers.rs`, `identify_misaligned_haps.rs` | T4 | 🔎 | open |
| PERF-4 | Every BAM record cloned into the per-island store | MED | Performance | `bam_lappers.rs:519` | T4 | 🔎 | open |
| PERF-5 | "Binary search" comment over a linear scan | MED | Performance | `structs.rs:115-128` | T0/T2 | 🔎 | open |
| PERF-6 | Mate recovered by re-reading the same region (N+1); unused locking | MED | Performance | `read_extraction/lib.rs:170-200` | T7/T9 | 🔎 | open |
| ROB-1 | `panic="abort"` → any panic hard-crashes the whole process | MED | Crash-safety | `Cargo.toml:10`, `read_extraction/Cargo.toml:19` | T9 | ✅ | **fixed (2026-06-11)** |
| ROB-2 | Non-text read-name crashes read-extraction instead of erroring | MED | Crash-safety | `read_extraction/lib.rs:223-249` | T9 | ✅ | **fixed (2026-06-11)** |
| ROB-3 | `M`-CIGAR input crashes inspection (documented invariant) | MED | Crash-safety | `pairwise_read_inspection.rs:91` | T1/T9 | ✅ | **mitigated (2026-06-11)** — non-fatal via ROB-1 unwind; panic→Result deferred to T0/T4 |
| ROB-4 | One `unwrap` that cannot actually fail (cosmetic) | LOW | Crash-safety | `build_phasing_graph/python_bindings.rs:172` | T2 | ✅ | open |
| VER-1 | Comment thresholds disagreed with code | MED | Verify | `read_extraction/lib.rs:52,93`; `bam_reading.rs:100,111` | T1 | ✅ | **comments fixed (2026-06-11)**; read_extraction code-vs-Python drift → T1 |
| HYG-1 | No lint gate (`clippy`) anywhere | MED | Hygiene | workspace + 3 crates | — | ✅ | **fixed (2026-06-11)** |
| HYG-2 | Stray `println!` writing to the program's output stream | LOW | Hygiene | `build_phasing_graph/python_bindings.rs:76` | — | ✅ | **fixed (2026-06-11)** |
| HYG-3 | Large blocks of "what this does" / Rust-tutorial comments | LOW | Hygiene | `graph_builder.rs:15-159`, `bam_reading.rs` | T2 | ✅ | open |
| HYG-4 | `get_` prefix on getters (non-idiomatic) | LOW | Hygiene | `structs.rs`, `haplotype_determination.rs` | T2 | 🔎 | open |
| HYG-5 | `new()` without `Default` | LOW | Hygiene | `build_phasing_graph/structs.rs` | T2 | 🔎 | open |
| HYG-6 | One leaked file handle per island (deliberate but unnecessary) | LOW | Hygiene | `bam_lappers.rs:103` | T0/T4 | ✅ | open |
| HYG-7 | Verbose `match` where `.is_ok()` reads cleaner | LOW | Hygiene | `read_extraction/lib.rs:54,64` | — | ✅ | **fixed (2026-06-11)** |
| CORR-1 | `is_sequencing_error` allele-depth indexed by ASCII (always ad=0) + missing-pos returns `true` | HIGH | Correctness | `phasing/haplotype_determination.rs:587-619` (+ prod `build_phasing_graph`) | 2nd pass | ✅ | **fixed (2026-06-14)** — re-validate (caveat in §F) |
| CORR-2 | `find_uncovered_regions` assumed sorted inspected intervals (Bug 4) | LOW | Correctness | `phasing/graph_builder.rs:564` (+ prod) | 2nd pass | ✅ | **fixed (2026-06-14)** |
| API-1 | `phase_bam` ignored its `reference` arg, used empty `PhaserParams.reference_genome` (Bug 1) | MED | Correctness | `phasing/lib.rs:200` | 2nd pass | ✅ | **fixed (2026-06-14)** |
| ROB-5 | HP `push_aux` aborts on BAMs already carrying an HP tag (Bug 3) | MED | Crash-safety | `phasing/hp_writer.rs:58`, `sdrecall/bam_filter.rs:90` | 2nd pass | ✅ | **fixed (2026-06-14)** |
| HYG-8 | READMEs install/list the deleted `build_phasing_graph` crate (Bug 5) | LOW | Hygiene/docs | `README.md`, `rust_modules/README.md` | 2nd pass | ✅ | **fixed (2026-06-14)** |

---

## A. Duplication across crates → handled by T0 foundation crates

The dominant redundancy is **between** crates, not within them. `build_phasing_graph` was migrated first; `haplotype_inspection` re-implemented the same low-level jobs rather than sharing. This is largely **expected debt** that the two planned foundation crates (`sdrecall-io`, `sdrecall-utils`) are designed to absorb — so most of this is "do it as part of T0", not "fix now".

### DUP-1 — Three copies of "read the BAM, drop noisy reads, build an interval index" (HIGH)
- `build_phasing_graph/src/bam_reading.rs:16-365`, `haplotype_inspection/src/bam_lappers.rs:39-618`, and a lighter third in `read_extraction/src/lib.rs:42-103`.
- They share, almost line-for-line: the `samtools collate` invocation, the qname-group streaming, and the "is this read noisy?" filter (MAPQ / alignment length / mate / median base-quality / soft-clip).
- **Plain language:** the same "load and clean a BAM file" code exists three times. If we ever change a filter rule, we'd have to change it in three places and hope they stay in sync.
- **Fix:** one BAM reader + one noisy-read filter in `sdrecall-io`. Prime candidate for that crate.

### DUP-2 — The whole per-read "vector" toolkit exists twice (HIGH)
- `build_phasing_graph/src/haplotype_determination.rs` vs `haplotype_inspection/src/pairwise_read_inspection.rs`: `extract_hap_vector`, `extract_error_vector`, the variant-counting helpers, the A/T/C/G/N encoder (**three** copies of the same map), the read-ID helper, and the caching wrappers.
- **Plain language:** the math that turns a read into the numeric vectors the algorithm works on is written twice, once per crate.
- **Fix:** one copy in `sdrecall-utils` (pure logic, no file I/O).

### DUP-2a — Those two copies have quietly DIVERGED (HIGH, correctness risk)
The two functions are `haplotype_inspection/src/pairwise_read_inspection.rs:83` (call it **STRICT**) and `build_phasing_graph/src/haplotype_determination.rs:147` (call it **LENIENT**). They were supposed to be the same port of Python's `get_hapvector_from_cigar` (`fp_control/pairwise_read_inspection.py:165-264`). They differ in **two** ways, verified against the Python on 2026-06-11:

1. **`M` (CIGAR op 0) handling** — STRICT `panic!`s on `M`; LENIENT treats `M` as a match (`push 1`). The Python asserts `operation != 0` (line 201) → **STRICT matches Python**; LENIENT does not. This divergence is *loud* (a panic), so it can't pass silently.
2. **Insertion-marker placement** — the dangerous, *silent* one. On an `I` op the marker value is `len*4` in both, but:
   - STRICT/Python **defer** the marker to the next reference-consuming op and write it at *that* op's first position (Python sets an `insertion` flag, then `hapvector[index] = insertion` on the next op). The marker lands on the reference base **after** the insertion.
   - LENIENT writes the marker **immediately, in place**, *overwriting the previous element*: `hap_vector[last_idx] = len*4`. The marker lands on the reference base **before** the insertion, and **clobbers whatever was there** (e.g. an adjacent `-4` SNV is destroyed).
   - Worked example, CIGAR `10= 1X 2I 5=`: STRICT/Python → `[1×10, -4, 8, 1,1,1,1]` (marker at index 11, the SNV at index 10 preserved); LENIENT → `[1×10, 8, 1,1,1,1,1]` (marker at index 10, **SNV lost**). Different vectors, off-by-one marker, and an erased variant.

- **Why LENIENT's validation still passed:** with `--eqx` alignment `M` never appears (so #1 never fires), and an insertion immediately adjacent to another variant is rare — so the divergent cases were not exercised. The divergence is real at the code level regardless.
- **Fix:** when T0/T4 merges these into one `sdrecall-utils` function, adopt the **STRICT** behavior (it is the Python-faithful one) — don't just pick a copy. LENIENT is effectively buggy relative to Python; verify on any read with an insertion adjacent to an SNV/indel during the T1 differential run.
- **Decisions (2026-06-11) — full record in [`module_vector_encoding.md`](module_vector_encoding.md) § Golden encoding:**
  - **DivA:** reject `M` (strict, Python-faithful) — minimap2 always emits `--eqx` CIGAR, so `M` is malformed input. Replace the `panic!` with a typed `Result` error that early-stops the job (fail-fast, clear message — ROB-3's deeper fix).
  - **DivB:** the *overwrite* of an adjacent base by the insertion marker loses information in **both** copies and in Python. Target encoding flips to **summation** with indel unit 10 (deletion `-10`/pos, insertion `+10*L` added to the prepending base; decode via nearest-10 delta). Entails: hap **padding sentinel moves `-10`→`-20`** (`HAP_PAD`, filter `>= -8`→`!= HAP_PAD`) to dodge the new deletion value; `count_snv` gains a compound pass (`v > 1 && v % 10 == 6`); `count_continuous_indel_blocks` predicate → `v == -10 || v > 1`. **Caveat:** diverges from Python numerically → keep current encoding through T1; switch in T4 validated by the **CIGAR/pileup oracle** (Python retired for these fns) + output-level set parity, all consumers updated in lockstep.

### DUP-3 / DUP-4 / DUP-5 (MED/LOW)
- DUP-3: interval-overlap is done three ways — a hand-rolled search in `build_phasing_graph/src/structs.rs:115`, `rust-lapper` in `bam_lappers.rs:663`, and endpoint sweeps in `graph_builder.rs`. Pick **one** interval backend in the foundation crate.
- DUP-4: "find shared variant positions, look up the alternate base, compare" is implemented in both `haplotype_determination.rs:1090` and `identify_misaligned_haps.rs:389`.
- DUP-5: small interval-merge and sort-then-dedup helpers overlap; fold into the interval module.

---

## B. Performance → the agenda for the 4–6× target (T4 + a benchmarked pass)

These are the levers for the headline speedup. They live in validated code, so they belong in a **measured** pass (before/after per-island timing), not piecemeal edits.

### PERF-1 — Big arrays copied on every cache hit, on the hottest path (HIGH) ✅ verified
- `identify_misaligned_haps.rs:607-710` (`stat_refseq_similarity`): the cached genome vector is `.clone()`d every time it's looked up (line 608); the per-read sequence data (four arrays) is cloned for every member read and the homolog (lines 686, 705); and slices are copied into owned arrays (`.to_owned()`, lines 643-644) only to be read.
- **Plain language:** we keep making full copies of data we only read. This runs for every read, in every region, for every candidate haplotype — i.e. the busiest loop in the pipeline.
- **Caveat:** these copies look like they were added to satisfy Rust's borrow checker (the same lookup tables are modified later in the same scope). Removing them needs a small restructure (compute-then-store, or an index/handle), **not** a blind find-and-replace.

### PERF-2 — The busiest loops run on a single core (HIGH)
- The per-haplotype loop (`:2000`) and per-region loop (`:2153`) are sequential; there's no parallelism anywhere despite the units being independent.
- **Plain language:** this is the biggest single speedup available, but the shared caches (PERF-1) currently force one-at-a-time execution. Fixing PERF-1's copying and giving each worker its own cache unlocks running haplotypes/regions in parallel.

### PERF-3 / PERF-4 / PERF-5 / PERF-6 (MED)
- PERF-3: hot maps/sets use the default (cryptographic, slower) hash; project convention is the fast-hash variants for these keys.
- PERF-4: `bam_lappers.rs:519` copies every read into the per-island store; check whether it can move them instead.
- PERF-5: `structs.rs:115` says "binary search" in the comment but does a linear scan — use a real binary search.
- PERF-6: `read_extraction` recovers a read's mate by re-scanning the same region (extra work per read), and wraps the writer in a lock it never needs (single-threaded).

---

## C. Crash-safety → matters much more once everything is one process (T9 + a robustness pass)

### ROB-1 — A crash is a *hard* crash (MED) ✅ verified — **fixed 2026-06-11**
- `panic="abort"` was set (`rust_modules/Cargo.toml:10`, `read_extraction/Cargo.toml:19`). This meant: if the Rust code hit an unexpected condition, the **entire Python process died instantly** with no error message and no traceback — it could not be caught and turned into a normal error.
- **Plain language:** today each genomic "island" runs in its own subprocess, so one bad island only kills that subprocess. But the north-star plan fuses everything into **one** Rust process — at that point, one bad island would take down the *whole run*.
- **Fix:** removed `panic="abort"` from the workspace root profile (panic strategy now defaults to `unwind`). For a PyO3 `cdylib` this is also the correct default — an unwinding panic is caught at the `#[pyfunction]` boundary and surfaced as a Python exception instead of aborting. Also deleted the redundant `[profile.release]` block from `read_extraction/Cargo.toml` (Cargo ignores profiles in non-root members; it was dead config that still implied `abort` and emitted a build warning). This single change makes **every** panic site — including ROB-2 and ROB-3 — catchable rather than fatal.

### ROB-2 / ROB-3 — Two specific "just crash" spots (MED) ✅ verified
- ROB-2 (`read_extraction/lib.rs`, FASTQ-write block) — **fixed 2026-06-11**: the read-name `std::str::from_utf8(qname).unwrap()` (R1 and R2) now returns a clean `PyValueError` (`"R{1,2} read name is not valid UTF-8"`) instead of panicking; the sequence `from_utf8(...).unwrap()` was changed to `String::from_utf8_lossy(...)` (sequence bytes are always ASCII, so output is unchanged but the panic path is gone).
- ROB-3 (`pairwise_read_inspection.rs:91`) — **mitigated 2026-06-11, deeper fix deferred**: a BAM aligned without the `=`/`X` CIGAR style still hits the `panic!` guard, but with ROB-1's `unwind` it is now caught at the PyO3 boundary and reported as a Python exception (mirrors the Python original, which `assert`s on the same condition). A full panic→`Result` conversion is **intentionally deferred to T0/T4**: this is the exact `extract_hap_vector` function being consolidated under DUP-2a, and changing its signature now would churn several validated hot-path call sites right before they are merged. Our pipeline always aligns with `--eqx`, so the guard does not fire on normal inputs.

### ROB-4 — A non-issue, listed for completeness (LOW) ✅ verified
- `build_phasing_graph/python_bindings.rs:172`'s `unwrap` can never actually fail (the value always exists by construction). Purely cosmetic — use `expect("…")` with a reason.

---

## D. Verify against Python (don't just "fix")

### VER-1 — Comments and code disagreed on thresholds (MED) ✅ verified against Python 2026-06-11
Checked each against the Python original before touching it; the two crates landed in **opposite** places:

- **`bam_reading.rs` (vs `fp_control/bam_ncls.py:178-190`) — code is correct, only the log *messages* were stale.** Code uses `low_qual_count >= 75` and `soft_clip >= 75`, both of which **match the Python** (`num_low >= 75`, `softclip >= 75`). The Rust debug strings said "≥ 50" and "≥ 20". Notably the soft-clip "≥ 20" was *already wrong in the Python message* (Python checks `>= 75` but prints "≥ 20") and was copied verbatim. **Fixed:** both Rust messages now print `>= 75` (lines 100, 111). Pure cosmetic; no behavior change.
- **`read_extraction/lib.rs` (vs `realign_recall/read_extraction.py:75-78`) — comments updated to match code, but the *code itself* genuinely diverges from Python. This stays a T1 item.** Per the user's instruction the comments were corrected to the code (`|AS − XS| ≤ 10` at line 52; fallback note `MAPQ < 60` at line 93). **But** the Rust filter logic does *not* match the Python production filter:
  - Rust: `mapq < 60` (always) **AND**, when `multi_aligned`, `![SA] && [XA] && |AS−XS| ≤ 10`.
  - Python `filter_expr`: `![SA] && ([XA] || mapq < 50)` — an **OR** between `[XA]` and `mapq < 50`, no `AS−XS` term, threshold 50 not 60.
  - (The Python source's own inline note at `read_extraction.py:77` describing the Rust as `![SA] && [XA] && abs(AS-XS) <= 5` is *also* stale — it predates the `≤ 10` value.)
  - These are different read-selection semantics and can select different read sets. **Action:** T1 differential validation must reconcile which is authoritative and align code (not just comments). Editing the comment here only removes the *internal* Rust contradiction; it does **not** resolve the Rust↔Python disagreement.

---

## E. What's already good (so we don't "fix" it)

- **POS-1:** zero `unsafe` blocks across all three crates; no raw pointers, no manual C-memory handling. The C libraries (htslib, the HiGHS solver) stay behind safe wrappers. There is **no memory-safety risk to clear before the fuse.**
- **POS-2:** `haplotype_inspection`'s Python entry point routes every fallible step into a clean Python error — zero panics on that boundary.
- **POS-3:** `bilc_solver.rs` and `pairwise_read_inspection.rs` are the cleanest modules (single-pass scans, pre-sized buffers, lookup tables, no hot copies).
- **POS-4:** no outdated idioms (`lazy_static`/`once_cell`/`try!`); pre-sizing, `format!`, and modern one-time-init are used correctly.

---

## Fixed in this review pass (2026-06-11, behavior-preserving)

These were safe to change without touching any logic:

1. **HYG-2** — removed a stray `println!` in `build_phasing_graph/src/python_bindings.rs` that wrote to the program's output stream and ignored the logging level.
2. **HYG-7** — simplified two verbose `match … { Ok => true, Err => false }` blocks to `.is_ok()` in `read_extraction/src/lib.rs`.
3. **HYG-1** — added a workspace lint gate: `[workspace.lints.clippy] all = "warn"` plus `[lints] workspace = true` in each crate. Deliberately conservative — only `clippy::all` (not `pedantic`, which would flood validated code with warnings), and **not** `unsafe_code = "deny"` (which would fight the PyO3 macros). Clippy lints don't affect normal builds, so this can't break compilation.

## Fixed in the follow-up pass (2026-06-11, at user request — `cargo check --workspace` clean, log `test_tmp/cargo_check_rob_ver_fixes.log`)

> **Committed + pushed** to `origin/rust-migration`: ROB/VER fixes in `b18c907` (panic/read_extraction) and `55dbfbd` (bam_reading messages). The contemporaneous haplotype_inspection WIP was also committed — `bd9f23b` (harness reorg `src/bin`→`examples/` + clippy opt-in), `3970373` (internals: `fast_median`→f32, `qname_to_node`→i32, dead-code removal, density cross-check test), `ab6fa3d` (README/copilot refresh). The migration docs + encoding/perf memory remain uncommitted by design (`docs/` is gitignored).

4. **ROB-1** — removed `panic="abort"` from the workspace root `[profile.release]`; deleted the dead (overridden) `[profile.release]` block from `read_extraction/Cargo.toml`. Panics now unwind and are catchable at the PyO3 boundary instead of aborting the host process. Also silences the prior "profile ignored for non-root member" build warning.
5. **ROB-2** — `read_extraction/src/lib.rs`: read-name `from_utf8(...).unwrap()` → clean `PyValueError`; sequence `from_utf8(...).unwrap()` → `String::from_utf8_lossy(...)`. No more panic on a non-UTF-8 qname.
6. **VER-1 (comments only)** — corrected the stale thresholds *in comments/log strings* to match the code: `bam_reading.rs` lines 100/111 (`>= 50`/`>= 20` → `>= 75`, confirmed matching Python); `read_extraction/lib.rs` lines 52/93 (`≤ 5` → `≤ 10`, `MAPQ < 50` → `MAPQ < 60`). **The read_extraction code↔Python logic divergence is NOT resolved by this — it remains a T1 item (see section D).**

Everything else above is **left in place** and tracked here for its owning task.

---

## Deliberately NOT changed (and why)

- **ROB-3 panic→`Result` conversion:** the M-CIGAR guard is now non-fatal (ROB-1 unwind), so the safety concern is addressed. The full conversion to a typed error is deferred to **T0/T4**, where this exact `extract_hap_vector` is consolidated under DUP-2a — converting its signature now would churn validated hot-path callers right before the merge.
- **read_extraction code↔Python filter divergence (VER-1):** the Rust filter (`mapq<60`, `|AS−XS|≤10`) genuinely differs from the Python production filter (`![SA] && ([XA] || mapq<50)`). Only the *Rust internal comments* were corrected; reconciling the *logic* is a **T1** differential-validation decision, not an out-of-band edit.
- **Tutorial/ASCII comment blocks (HYG-3):** large but harmless; bulk-deleting comments from validated code is opinionated and best done with the owning crate's parity re-confirm (T2).
- **`get_` renames (HYG-4), the duplications (DUP-\*), the performance items (PERF-\*):** all change real structure or behavior in validated code, so they belong to their owning tasks (T0/T2/T4) where they can be re-validated, not done out of band. See the PERF-1/PERF-2 fix plan handed to the user (owned by T4) for the agreed approach.

---

## F. Second review pass (2026-06-14) — `phasing` crate (post `build_phasing_graph`→`phasing` fuse)

A second external review of the migrated `phasing` crate (which absorbed `build_phasing_graph`) + the `sdrecall` orchestrator surfaced five issues. All **fixed + pushed** on `rust-migration`: `6bf5f5a` (API-1), `583e4b1` (CORR-1), `c036bad` (ROB-5), `3186d88` (HYG-8), `0d12963` (CORR-2); a sixth item (phasing-wiring duplication) was de-duplicated in `ee31d18`. `cargo test` clean — phasing 23 + fp-control 2, 0 warnings (`/paedyl01/disk1/yangyxt/test_tmp/cargo_test_bug4_unify_20260614.log`).

### CORR-1 — `is_sequencing_error` allele-depth lookup was a no-op + missing-position polarity flip (HIGH) ✅
The consequential one. The per-base allele-depth lookup used `get_allele_depth(pos_data, target_base as usize)`, but `target_base` is the **raw ASCII** base from `record.seq().as_bytes()` (65..=84) and `get_allele_depth` returns 0 for any index ≥ 5 — so `ad` was **always 0**, `af` always 0.0, the `af ≤ 0.02` artifact test always true, and tolerance collapsed to "is the base Q<13", ignoring the real allele frequency. Separately, a position absent from the pileup returned `true` (tolerate) where Python `seq_err_det_stacked_bases` returns `False`. Both made the Rust over-permissive → it merged reads carrying a **real low-quality variant** onto one haplotype.
- **Fix:** encode the base via the canonical `base_to_index` (the same map `build_allele_depth_map` fills the array with) before indexing; flip the missing-position return to `false`. Now matches the original Python algorithm and the pileup ground truth.
- **⚠ Validation caveat (read before trusting the old parity numbers):** T4 **246/246** and T2/T3 used a **bug-derived reference** — `diff_vs_hybrid` feeds a Python partition dumped from a run that used the *buggy* production `build_phasing_graph`, and the fused Rust graph was a byte-identical copy with the same bug, so they agreed by construction. With CORR-1 fixed the fused graph now **diverges** from that frozen dump on islands with low-quality SNV mismatches — that divergence **is the correction**, not a regression. To re-validate apples-to-apples: mirror CORR-1 into the production crate (below), regenerate the dump, re-run `diff_vs_hybrid`; and/or check directly against a samtools-pileup oracle on the mismatch positions. Live SDrecall precision/recall reflect the buggy behavior → a GIAB re-benchmark is warranted (expected: precision ↑).
- T3's **270/270 is unaffected** — it feeds a *dumped weight matrix* into `phase()`, never running the graph builder, so `is_sequencing_error` is off that path.

### CORR-2 — `find_uncovered_regions` assumed sorted inspected intervals (Bug 4, LOW) ✅
The per-edge weight dedup subtracted already-inspected sub-regions with a forward cursor sweep — correct only if the inspected intervals are start-sorted. But `get_overlap_intervals` emits the ≤4 read×read overlaps in read-combination order and `FastIntervals::add` stores them verbatim, so with overlapping mates they can arrive unsorted → a region re-counted or skipped → one edge's weight perturbed. **Fix:** filter to the query-overlapping intervals and sort by start before the sweep (no merge needed; the cursor uses `max()`). Narrow trigger + quantitative impact (a weight nudge, not a Same/Different flip), so unlike CORR-1 it won't meaningfully move the partition. Unit tests added.

### API-1 — `phase_bam` ignored its `reference` argument (Bug 1, MED) ✅
`phase_bam(bam, _reference, …)` used `params.reference_genome` (empty under `PhaserParams::default()`) for the bcftools-mpileup map — the documented library call would feed bcftools an empty `--fasta-ref` and fail; the CLI worked only because it set the field. **Fix:** use the `reference` argument; delete the redundant `PhaserParams.reference_genome` field (single source of truth).

### ROB-5 — HP `push_aux` aborts on already-tagged BAMs (Bug 3, MED) ✅
rust-htslib 0.47.1 `push_aux` returns `BamAuxTagAlreadyPresent` on a duplicate tag, so writing `HP:Z:…` to a BAM already carrying an HP tag (pre-phased / re-annotated / idempotent rerun) aborted the whole BAM. **Fix:** `remove_aux(b"HP")` (ignore not-present) before `push_aux`, at both sites (`hp_writer.rs`, `sdrecall/bam_filter.rs`).

### HYG-8 — READMEs referenced the deleted `build_phasing_graph` crate (Bug 5, LOW) ✅
`README.md` (`pip install build-phasing-graph`) and `rust_modules/README.md` (crate table) still presented the absorbed crate as separately installable. **Fix:** removed; role folded into the `phasing` crate's row/description.

### Also done — phasing-wiring de-duplication (commit `ee31d18`)
`phase_bam` and `fp_control::run_fp_control` both wired BAM→allele-depth→graph→phase (with an inlined `build_vertex_qname`). Extracted `phasing::build_and_phase(bam, reference, &PhaserParams) → Option<PhasedReads>` as the single phasing core; both callers now apply only their own downstream gates (the phaser writes the HP-tagged/unphased BAM; fp-control applies the ≤2-vertex/≤2-haplotype shortcuts + inspects). Behavior-preserving (unit tests pass); full fp-control equivalence rides on the same `diff_vs_hybrid` re-run as CORR-1.

### Production mirror (branch `rust-impl`) — staged for the user, not committed out of band
CORR-1 and CORR-2 exist **verbatim** in the live PyO3 crate `rust_modules/build_phasing_graph/src/{haplotype_determination,graph_builder}.rs` (API-1/ROB-5 are N/A there — the old crate has no `phase_bam`/`hp_writer`). Mirroring them is needed for the apples-to-apples differential and improves live precision, but it changes **deployed** behavior (wheel rebuild + GIAB benchmark) and the `rust-impl` worktree currently holds unrelated uncommitted work — so the patch is handed to the user to apply, not committed here.
