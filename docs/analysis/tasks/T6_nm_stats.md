# T6 — `nm-stats` crate (`cal_edge_NM_values`)

**Crate:** `nm-stats` (new — lib + bin)
**Status (2026-07-17):** Algorithm validated, but production integration is incomplete. The orchestrator currently omits the NM cutoff scan rather than applying `nm_stats::nm_distribution_poisson`.
**Depends on:** T0
**Replaces (Python):** `realign_recall/cal_edge_NM_values.py` (`calculate_NM_distribution_poisson`).

## Goal & scope boundary

Compute the NM (edit-distance) distribution cutoff used to flag noisy edges. Stream proper-pair, MAPQ-60, non-`XA` reads; for each compute `NM_tag − max_indel_gap` ("scattered edit distance"); take the mean; then raise an integer `cutoff` until `poisson.cdf(cutoff, mean) ≥ 1 − conf_level`. Self-contained and easy — a clean first standalone-crate exercise.

Out of scope: anything beyond the (cutoff, mean) scalars.

## Data flow

```
BAM ─→ iterate proper-pair MAPQ60 non-XA reads (up to ~3M)
     ─→ per read: NM_tag − max_indel_gap  → f64 array
     ─→ mean
     ─→ poisson cdf loop → integer cutoff
     ─→ (cutoff, mean)
```

## Dependencies

- Crates: `rust-htslib` 0.47.0 (fetch + tag/CIGAR parse), `statrs` 0.18.0 (Poisson), `sdrecall-utils`.
- External tools: none.

## Performance bottleneck / rationale

The cost is the multi-million-record pysam fetch + tag parse; the scipy Poisson part is negligible. `rust-htslib` accelerates the per-read scan well, and the task is fully self-contained (no cross-crate data), making it a good early win and a template for the lib+bin pattern.

## Tests

### Unit (tier 1)
- Mean + cutoff on a synthetic record set with known NM/indel values.
- Poisson cutoff loop against hand-computed cdf thresholds.

### Differential vs Python (tier 2)
- Run on a real BAM; compare (cutoff, mean) to the Python function.

**Pass criterion:** identical integer `cutoff`; `mean` within `1e-9`.

**Data:** a production realigned BAM (HG002/HG006).

## Progress
- [x] Scaffold `nm-stats` crate (lib + `--bam --conf-level --sample-size --threads` bin) — `rust_modules/nm-stats/` (+ minimal `sdrecall-utils` foundation crate: `SdError`/`Result`, `configure_parallelism`, `init_console_logger`)
- [x] Per-read NM−indel scan (rust-htslib) — `passing_scatter_dist(&Record)` (one fused filter+map unit; zero per-record alloc, reused `Record` buffer)
- [x] Poisson cutoff loop — `poisson_cutoff(mean, conf_level)` via `statrs::Poisson` + `DiscreteCDF`; `max(_,4)` floor; `mean<=0` short-circuit (statrs rejects λ≤0 where scipy returns cdf=1)
- [x] Unit tests — **10 pass** (9 nm-stats + 1 sdrecall-utils), clippy clean. Logs: `/paedyl01/disk1/yangyxt/test_tmp/nm_stats_test_release_20260611.log`, `…/nm_stats_clippy_20260611.log`
- [x] Differential vs Python on real BAM — **PASS**
- [ ] Wire the returned Poisson cutoff into the orchestrator FP-control/noisy-edge path and add an end-to-end regression proving the filter is applied.

### Differential result (2026-06-11) — PASS (bit-identical)

Ran the Python `calculate_NM_distribution_poisson` and the Rust `nm-stats` CLI on the same real BAM `example/HG002.CMRG.hg38.test.bam` (coord-sorted, 229,943 reads → 54,904 usable; both sides use the **full** passing set since n < 3M sample_size, so the fetch()-vs-read() ordering risk R1 is moot):

- **cutoff:** py `4` == rust `4`
- **mean:** py `0.6937928019816407` == rust `0.6937928019816407` — **bit-identical f64** (`0x1.6338cf656c363p-1`, |diff| = 0), exceeding the "within 1e-9" criterion.
- Driver: `/paedyl01/disk1/yangyxt/test_tmp/diff_nm.py` (calls the real Python fn); compare log: `/paedyl01/disk1/yangyxt/test_tmp/diff_nm_HG002_20260611.log`.
- Note: the SDrecall Python logger writes to **stdout**, so a differential harness must take the last stdout line (the `cutoff\tmean`); the Rust CLI keeps logs on stderr and stdout clean.
- Remaining for full sign-off: HG006 (needs Zenodo BAM) — the design's second differential sample.

**Template established.** This is the reusable lib+bin pattern (streaming `read(&mut rec)` scan, one fused filter+map unit, thin numeric tail, typed `SdError`, clap bin printing the Python-dump format) that T5/T7–T9 copy.

---

## Migration design — frontier propagation (2026-06-11)

This crate is the **template** every later lib+bin stage (T5, T7–T9) copies: a streaming
record scan with **zero per-record allocation**, a single versatile filter+map unit,
a thin CLI, and a typed error via `sdrecall-utils::SdError`. Keep it spartan.

### 1. Python logic inventory

Single function, `realign_recall/cal_edge_NM_values.py`:

| fn | lines | what it does | control-flow shape |
|----|-------|--------------|--------------------|
| `calculate_NM_distribution_poisson(bam, conf_level=0.01, sample_size=3_000_000, logger)` | 7–45 | streams BAM records, filters to proper-pair / primary / non-dup / non-qcfail / no-`XA` / MAPQ==60, computes per-read `NM − max_indel_gap`, takes the mean, then linear-searches the smallest integer `cutoff` with `poisson.cdf(cutoff, mean) ≥ 1−conf_level`, floors it at 4 | one **sequential** `for read in bam_handle.fetch()` pass (no region, no parallelism — `fetch()` with no args = whole-file linear scan), early `break` at `n >= sample_size`, then a tiny `while` cdf loop |

Line-level semantics to port faithfully:
- **L11** `nm_array = np.empty(sample_size)` — pre-allocated `f64` buffer of length `sample_size`; only the first `n` slots are written, then sliced (L36) before `np.mean`.
- **L14** `bam_handle.fetch()` — *no region argument*. For an indexed BAM this iterates **mapped reads** in coordinate order; functionally equivalent to a plain linear "read every record" pass. (See risk R1 on unmapped tail — irrelevant here because every filter read is mapped.)
- **L15–21** filter conjunction, in this exact order: `is_proper_pair` AND NOT `is_secondary` AND NOT `is_supplementary` AND NOT `is_duplicate` AND NOT `is_qcfail` AND `len([t for t in get_tags() if t[0]=="XA"])==0` (i.e. **no `XA` tag present**) AND `mapping_quality == 60`.
- **L22–23** `gap_sizes = [t[1] for t in cigartuples if t[0] in (1,2) and t[1] > 1]; max_gap = max(gap_sizes) if gap_sizes else 0` — over CIGAR ops, op code `1`=Ins, `2`=Del; collect lengths **strictly > 1**; take the max (`0` if none). This is the largest single indel gap, i.e. "the read's edit distance minus its one biggest contiguous indel".
- **L24–26** `edit_dist = read.get_tag("NM"); scatter_edit_dist = edit_dist − max_gap` → written to `nm_array[n]`. **`NM` is assumed present** on every passing read (it is, for BWA/minimap output); a `KeyError` in Python = hard failure. We mirror that: missing `NM` on a passing read is a typed error, not a silent skip.
- **L28–29** `if n >= sample_size: break` — cap the sample at `sample_size` reads (default 3M).
- **L34–36** under-sample warning + slice.
- **L38** `nm_mean = np.mean(nm_array)` over the `n` written values.
- **L39–41** `cutoff = 0; while poisson.cdf(cutoff, nm_mean) < 1 − conf_level: cutoff += 1` — smallest integer cutoff whose Poisson CDF reaches the `1−conf_level` quantile.
- **L44** `cutoff = 4 if cutoff < 4 else cutoff` — **floor the cutoff at 4** (port this; it is a real business rule, not a Python artefact).
- **L45** returns `(cutoff, nm_mean)`.

Note the current production caller (`misalignment_elimination.py:132`) has this **commented out** — so this is a faithful logic port with no live consumer to break; T9 wiring will decide where it re-enters.

### 2. Python → Rust crate mapping

| Python operation / idiom | Rust crate::api | confidence |
|---|---|---|
| `pysam.AlignmentFile(bam,"rb")` | `rust_htslib::bam::Reader::from_path(&Path)` + `Read::set_threads` | verified-docs |
| `for read in bam_handle.fetch()` (whole-file linear) | `while let Some(r) = Read::read(&mut reader, &mut rec) { r?; … }` — reuses one `bam::Record`, returns `Option<Result<()>>`, `None`=EOF | verified-docs |
| `read.is_proper_pair` / `.is_secondary` / `.is_supplementary` / `.is_duplicate` / `.is_qcfail` | `Record::is_proper_pair()` / `is_secondary()` / `is_supplementary()` / `is_duplicate()` / `is_quality_check_failed()` — all `(&self)->bool` | verified-docs |
| `read.mapping_quality == 60` | `Record::mapq() -> u8` (compare `== 60`) | verified-docs |
| `len([t for t in get_tags() if t[0]=="XA"])==0` (no XA) | `Record::aux(b"XA").is_err()` (htslib returns `Err(Errors::BamAuxTagNotFound)` when absent) | verified-docs |
| `read.get_tag("NM")` | `Record::aux(b"NM")` → match `Aux::U8/U16/U32/I8/I16/I32(v) => v as i64`; any other/absent ⇒ `SdError::Htslib` | verified-docs |
| `read.cigartuples` (op,len) iteration | `Record::cigar() -> CigarStringView`, iterate `Cigar::Ins(len)|Cigar::Del(len)` (zero-copy alt: decode `raw_cigar() -> &[u32]` inline, op=`c & 0xf`, len=`c >> 4`) | verified-docs |
| `t[0] in (1,2) and t[1] > 1` → `max(...)` else 0 | fold over cigar ops: `max_indel_gap = ins/del lens where len>1, else 0` (`u32`) | verified-docs |
| `np.empty(sample_size)` + write-then-slice | `Vec::<f64>::with_capacity(sample_size)`; `push` only on pass; len = `n` (no slice needed) | verified-docs |
| `np.mean(nm_array[:n])` | `sum / n as f64` (Welford unnecessary; values are small ints) | verified-docs |
| `scipy.stats.poisson.cdf(cutoff, mean)` | `statrs::distribution::{Poisson, DiscreteCDF}`: `Poisson::new(mean)?.cdf(cutoff_u64)` | verified-docs |
| `1 - conf_level` quantile loop | `while pois.cdf(cutoff) < 1.0 - conf_level { cutoff += 1 }` (`cutoff: u64`) | verified-docs |
| `logger.warning/info` | `log::{warn,info}` (front-end set by `sdrecall_utils::init_console_logger`) | plausible |

External-dependency replacement summary: **pysam → rust-htslib** (in-process, no subprocess), **scipy/numpy → statrs + plain f64 arithmetic** (no ndarray needed — a single mean over a `Vec<f64>` doesn't justify it). No pybedtools / networkx / graph-tool in this task.

### 3. Crate file layout

```
nm-stats/
├─ Cargo.toml          # lib + bin; deps: rust-htslib 0.47, statrs 0.18, sdrecall-utils,
│                      #   clap 4.6 (bin only, behind no feature — it's a binary dep), log
├─ src/
│  ├─ lib.rs           # pub fn nm_distribution_poisson(...) + pub struct NmCutoff; re-exports
│  └─ bin/
│     └─ nm_stats.rs   # clap CLI: --bam --conf-level [--sample-size] → prints "cutoff\tmean"
└─ examples/
   └─ diff_nm.rs       # differential harness: runs lib on a real BAM, prints (cutoff,mean)
                       #   in the same format the Python dump emits, for set/scalar compare
```

**One-versatile-unit-per-job.** The whole task is a *filter → scalar-map → reduce* over a
record stream. There is exactly **one** load-bearing unit:

- `fn passing_scatter_dist(rec: &bam::Record) -> Result<Option<f64>>` — encodes the L15–26
  filter+map for a single record: returns `Ok(None)` if the read is filtered out, `Ok(Some(d))`
  with `d = NM − max_indel_gap` if it passes, `Err` if a passing read is missing `NM`. This is
  the single place the predicate lives; the streaming driver and every unit test reuse it.

The Python has **no duplicated path** to collapse here (it's already one loop), so the
discipline is *not to introduce* a second helper: do **not** split "is this read eligible"
and "compute its distance" into two functions that each re-walk the record — fuse them into
`passing_scatter_dist` so the CIGAR/tag access happens once per record. The driver
`nm_distribution_poisson` is genuine upstream orchestration (open reader, loop, mean, cdf),
and the `bin` is the file-in/scalar-out wrapper — both legitimately distinct from the unit.

### 4. Core data structures + key fn signatures

```rust
use rust_htslib::bam::{self, Read, record::{Aux, Cigar}};
use std::path::Path;
use sdrecall_utils::{Result, SdError};   // typed error + Result alias

/// (cutoff, mean). Newtype-light: a plain struct keeps the bin/Python-dump format obvious.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NmCutoff { pub cutoff: u64, pub mean: f64 }

/// THE unit. Borrows the record (read-only, no clone); returns an owned scalar or None.
/// WHY &Record: we never mutate and never keep the record past the call — a shared borrow
/// is zero-cost and lets the driver reuse one Record buffer across the whole scan.
/// WHY Option inside Result: filtered-out (normal) is `Ok(None)`; missing-NM (a real data
/// defect on a read we accepted) is `Err`, matching Python's KeyError-is-fatal semantics.
fn passing_scatter_dist(rec: &bam::Record) -> Result<Option<f64>>;

/// Driver / orchestration. Takes the path by &Path (caller owns the buffer); scalars by value.
/// WHY &Path not PathBuf: we only read it to open the reader — no ownership needed.
/// WHY threads param: rust-htslib BGZF decompression parallelizes the I/O-bound scan; pass
/// the budget from sdrecall_utils::configure_parallelism upstream.
pub fn nm_distribution_poisson(
    bam_path: &Path,
    conf_level: f64,        // default 0.01 at the call site, not here (Rust has no kwargs)
    sample_size: usize,     // default 3_000_000 at the call site
    threads: u8,
) -> Result<NmCutoff>;

/// Pure numeric tail, split out ONLY because it is independently unit-testable against
/// hand-computed cdf thresholds and has zero I/O. Reused by the driver.
/// WHY owned f64 in / NmCutoff out: trivially-copyable scalars; no borrow story.
fn poisson_cutoff(mean: f64, conf_level: f64) -> Result<u64>;  // returns max(cutoff, 4)-able value; floor applied in driver or here — see R4
```

`max_indel_gap` helper is **inlined** into `passing_scatter_dist` (a 6-line fold over
`rec.cigar()`); it is not a public function — extracting it would be a near-duplicate of the
one place CIGAR is walked and violates the no-overlap rule.

Driver skeleton (the reusable streaming template):

```rust
let mut reader = bam::Reader::from_path(bam_path)
    .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?;
reader.set_threads(threads as usize).ok();          // best-effort decompression threads
let mut rec = bam::Record::new();                    // ONE buffer, reused every iteration
let mut acc = Vec::<f64>::with_capacity(sample_size.min(1 << 20)); // cap pre-alloc
let mut n = 0usize;
while let Some(r) = reader.read(&mut rec) {           // zero-alloc streaming
    r.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
    if let Some(d) = passing_scatter_dist(&rec)? {
        acc.push(d);
        n += 1;
        if n >= sample_size { break; }
    }
}
if n == 0 { return Err(SdError::InsufficientPairs(0)); } // mean of empty = NaN; fail loudly
if n < sample_size { log::warn!("BAM {} only had {n} usable reads", bam_path.display()); }
let mean = acc.iter().sum::<f64>() / n as f64;
let cutoff = poisson_cutoff(mean, conf_level)?;
Ok(NmCutoff { cutoff, mean })
```

### 5. Performance optimizations from ownership / borrowing

- **Zero per-record allocation.** `read(&mut rec)` decodes into one reused `bam::Record`
  (vs pysam, which mints a Python object per read). Over 3M records this removes 3M
  allocations — the dominant win and the reason this is the template.
- **Borrow, never clone, in the unit.** `passing_scatter_dist(&rec)` only reads flags/tags/
  cigar; nothing escapes, so no `.clone()` ever appears in the hot path.
- **Cheap-predicate-first short-circuit.** Order the filter so the *cheapest* bitflag checks
  (`is_proper_pair`, `is_secondary`, …, `mapq()==60`) run before the comparatively pricier
  `aux(b"XA")` / `aux(b"NM")` lookups and the CIGAR walk — matching Python's order but now
  the `&&` short-circuit actually skips work on the ~majority of rejected reads.
- **`raw_cigar()` over `cigar()` for the gap fold (optional micro-opt).** `cigar()` allocates
  a fresh `CigarStringView` (docs: "creates a fresh copy"); decoding `raw_cigar() -> &[u32]`
  inline (`op = c & 0xf`, `len = c >> 4`, Ins=1/Del=2) is fully zero-copy. Start with the
  readable `cigar()` form; switch to `raw_cigar()` only if profiling flags the alloc.
- **Capped pre-allocation.** `Vec::with_capacity(sample_size.min(1<<20))` avoids reserving a
  24 MB f64 buffer up front when most BAMs hold far fewer usable reads; it grows amortized-O(1)
  to the real `n`. (Python's `np.empty(3M)` always reserves the full buffer.)
- **Single sequential pass, no rayon.** The scan is I/O/decompression-bound and inherently
  ordered (early-break at `sample_size`); the parallelism that matters is htslib's internal
  BGZF threads via `set_threads`, not record-level rayon. Mean + cdf loop are O(n)/O(cutoff)
  and negligible — do **not** over-engineer with a parallel reduce.
- **`f64` accumulation, no ndarray.** A single mean doesn't justify pulling `ndarray`/`numpy`;
  `acc.iter().sum()` is one tight loop. (FxHashMap/ahash not needed — there are no maps here.)

### 6. Risks / open decisions

- **R1 — `fetch()` vs linear `read()` parity.** Python's `fetch()` with no region iterates
  **mapped** reads in coord order and **skips the unmapped tail**; our `read()` loop reads
  *every* record including the unmapped tail. This is parity-safe because every read that
  passes the L15–21 filter is necessarily mapped (proper-pair + MAPQ 60), so the unmapped tail
  contributes nothing — but the *order* of sampling differs from `fetch()` only if a BAM is
  not coordinate-sorted. **Decision needed:** confirm production BAMs here are always
  coord-sorted (they are post-`samtools sort`); if so, `read()` and `fetch()` sample the same
  first-`N` set and `(cutoff, mean)` match exactly. If not, the sampled subset (hence mean)
  could differ. Flag for the differential run.
- **R2 — NM aux integer type.** SAM spec stores `NM` as an unsigned int, but BWA/minimap may
  emit `i` (signed). The match must accept **all** of `Aux::U8/U16/U32/I8/I16/I32`; a `Float`/
  absent `NM` on a passing read ⇒ `SdError::Htslib` (mirrors Python's fatal `KeyError`). Do
  **not** add a "skip read if NM missing" fallback (dependency-availability rule: fail clearly).
- **R3 — `cutoff` type for `cdf`.** `statrs::Poisson` implements `DiscreteCDF<u64, f64>`, so the
  cutoff loop variable is `u64`, and the `DiscreteCDF` trait must be `use`d for `.cdf()` to be
  in scope. Return as `u64` (or cast to the int type T9's caller expects). Verified on docs.rs.
- **R4 — where the `cutoff = max(cutoff, 4)` floor lives.** Python applies it *after* the cdf
  loop (L44). Keep it in **one** place — apply inside `poisson_cutoff` so the unit returns the
  business-final value, OR in the driver right before constructing `NmCutoff`. Decide once;
  don't apply it twice. (Recommendation: inside `poisson_cutoff`, so the scalar tail is
  self-contained and testable as one unit.)
- **R5 — empty-sample behaviour.** Python computes `np.mean([])` = `nan` and then loops the cdf
  forever / errors. We **diverge intentionally** for safety: `n == 0 ⇒ Err(SdError::InsufficientPairs(0))`.
  Confirm no production BAM legitimately yields zero proper-pair MAPQ-60 reads (it shouldn't).
- **R6 — Poisson cdf numeric parity vs scipy.** statrs and scipy both compute the regularized
  incomplete gamma; agreement to ~1e-12 is expected, and since `cutoff` is the *integer* where
  the cdf first crosses `1−conf_level`, sub-1e-9 differences cannot flip the integer except on
  a measure-zero boundary. Differential test asserts integer `cutoff` equality and `mean`
  within 1e-9 (already the doc's pass criterion).
- **R7 — `set_threads` return.** `Read::set_threads` returns `Result`; treat failure as
  non-fatal (`.ok()` / log) — thread-pool unavailability must not abort the scan, but it is
  *not* a missing-dependency fallback (the scan still runs single-threaded by htslib default).

### 7. Test plan delta

**Unit (tier 1) — hand-built `bam::Record`s, no file I/O:**
- `passing_scatter_dist`:
  - *passes*: proper-pair, primary, non-dup, non-qcfail, no `XA`, MAPQ 60, `NM=5`,
    CIGAR `10=2I3=4D5=` → max indel gap = `max(2,4)=4` → returns `Ok(Some(1.0))` (5−4).
  - *indel len == 1 excluded*: CIGAR with a `1I`/`1D` and `NM=3`, no gap>1 → gap=0 →
    `Ok(Some(3.0))` (confirms the strict `>1` from L22).
  - *each filter rejection*: flip exactly one of {not proper-pair, secondary, supplementary,
    dup, qcfail, MAPQ 59, has `XA`} → `Ok(None)` (7 cases, one assert each).
  - *missing NM on a passing read* → `Err(SdError::Htslib(_))` (matches `assert!(matches!(...))`).
- `poisson_cutoff`:
  - `mean=1.0, conf_level=0.01` → hand-compute: cdf(0)=0.3679, cdf(1)=0.7358, cdf(2)=0.9197,
    cdf(3)=0.9810, cdf(4)=0.9963 ≥ 0.99 → raw cutoff = 4 → final 4.
  - `mean=0.05` → raw cutoff resolves to 1, **floored to 4** → assert `==4` (locks the L44 floor).
  - `mean=8.0, conf_level=0.001` → larger cutoff (> 4) → assert the floor does **not** lower it.
- `nm_distribution_poisson` end-to-end on a tiny **synthetic indexed BAM** written in the test
  (3–4 reads with known NM/CIGAR): assert `NmCutoff { cutoff, mean }` equals the
  hand-computed pair; assert under-sample `warn!` path when `n < sample_size`.

**Differential (tier 2) — `examples/diff_nm.rs`:**
- Run `nm_distribution_poisson` on a production realigned BAM (HG002 + HG006), print
  `cutoff\tmean` to stdout; run the Python `calculate_NM_distribution_poisson` on the same BAM
  via the existing dump harness; assert **identical integer `cutoff`** and **`|mean_rs − mean_py| < 1e-9`**
  (the doc's stated pass criterion). Capture `RUST_LOG=debug` output to
  `/paedyl01/disk1/yangyxt/test_tmp/diff_nm_<sample>.log` per the project's test-verification rule.
- Sanity probe for R1: if integer cutoffs ever disagree, dump the first-`N` sampled NM values
  from both sides to confirm it is a sampling-order (sort) issue, not a filter/encoding bug.
