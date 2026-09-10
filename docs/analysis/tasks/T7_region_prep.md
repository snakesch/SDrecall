# T7 — `region-prep` crate (`prepare_masked_align_region` projection)

**Crate:** `region-prep` (new — lib + bin)
**Status (2026-07-17):** Production-integrated. The focused HG002 differential produced byte-identical FC and 576/576 NFC BEDs; the stage is exercised by the validated three-assembly orchestrator runs.
**Depends on:** T0
**Replaces (Python):** `realign_recall/prepare_masked_align_region.py` (`extract_and_pad_segments` and the per-RG fc/nfc region projection).

## Goal & scope boundary

Per realignment group (RG), compute the **fc** (functional/target) and **nfc** (non-functional/counterpart) realignment BEDs via strand-aware **relative-coordinate segment projection + padding + merge**. This is genuine in-Python compute (not a bedtools wrapper), which is why it's worth a dedicated crate rather than folding into the orchestrator.

Out of scope: the surrounding pybedtools intersect/merge bookkeeping that has direct interval-library equivalents (handled via `sdrecall-io` interval ops).

## Data flow

```
query BED + counterpart BED + ref coordinates (per RG)
   → project segments into relative coordinates (strand-aware)
   → pad + merge
   → fc_bed, nfc_bed (per RG subgroup)
```

These BEDs feed Phase 2a realignment (read extraction region + masked-align region).

## Dependencies

- Crates: `bedrs` 0.2.26 (interval ops, via `sdrecall-io`); interval math on `sdrecall-utils::GenomicInterval`.
- External tools: none (replaces the custom Python projection; bedtools-equivalent ops use the interval library).

## Performance bottleneck / rationale

Moderate per-RG CPU in Python coordinate projection. Not the top hotspot, but it's real compute (not a thin wrapper), so it belongs in Rust for the full-pure-Rust goal and removes a pybedtools dependency.

## Tests

### Unit (tier 1)
- Segment projection on +/− strand cases; padding at contig boundaries; merge of adjacent/overlapping projected segments.

### Differential vs Python (tier 2)
- For each RG, compare Rust fc/nfc BEDs to the Python output.

**Pass criterion:** identical fc/nfc BED intervals (sorted) per RG.

**Data:** per-RG BEDs from a prepared HG002/HG006 run.

## Progress
- [x] Scaffold `region-prep` crate (lib + bin) — `rust_modules/region-prep/` (`all_region_bed.rs`, `project.rs`, `per_rg.rs`, bin, `examples/diff_region_prep.rs`)
- [x] Strand-aware relative-coordinate projection — `project::extract_and_pad_segments` (pure, borrow-in/owned-out; the two Python strand branches collapsed into one `match same_strand`)
- [x] Pad + merge — `NFC_PAD=600` / `FC_SLOP=500` separate consts; book-ended-fusing merge
- [x] Unit tests (strand/boundary/merge) — **26 pass**, clippy clean (`/paedyl01/disk1/yangyxt/test_tmp/region_prep_test4_20260612.log`)
- [x] Differential fc/nfc BED vs Python per RG — **PASS, byte-identical**

### Result (2026-06-12) — 26 unit tests + real-data differential PASS (byte-identical)

- **Differential on real HG002 RG0 (all 576 subgroups):** FC target identical (77 ivs, cov 1042172) + **576/576 NFC BEDs byte-identical** vs the actual Python `prepare_masked_align_region` module. Log `/paedyl01/disk1/yangyxt/test_tmp/region_prep_diff_full2_20260612.log`; reference generator `/paedyl01/disk1/yangyxt/test_tmp/gen_region_prep_ref.py`.
- **Verbatim-port hazards handled:** R1 OR-predicate (verbatim), R2 unstranded final merge, R3 mtime+coverage reuse, R5 basename-derived `rg_label`, R6 distinct pad/slop consts, R7 required `.fai`, R8 disjoint-FC → empty Vec (documented divergence from Python's `IndexError`).
- **Two findings the differential forced (both now resolved):**
  - **R2-bis → fixed in `sdrecall-io`:** the differential proved pybedtools `.merge()` (= `bedtools merge -d 0`) **fuses book-ended** intervals. `sdrecall-io::sort_merge_bed` had been wrongly made *strict-overlap* (a propagated mis-spec); **corrected 2026-06-12** to fuse book-ended (`next.start <= prev.end`), matching pybedtools. All 52 sdrecall-io tests + the corrected bookended pin pass. (T7's local `merge_bookended` workaround is now redundant but harmless.)
  - **R9 (new):** Python writes the per-subgroup FC bed as **3 columns** → `main_interval.strand` is `'.'`, so the opposite-strand (reverse-complement) projection branch runs for nearly every NFC row. Matched by feeding the FC interval as `Strand::Unknown`. Without R2-bis+R9, coverage matched but 21/576 subgroups' interval positions diverged.

## Migration design — frontier propagation (2026-06-11)

> Coding-grade design for the `region-prep` crate. Builds on the T0 `sdrecall-utils`/`sdrecall-io`
> interface (see CLAUDE.md / the prompt's interface contract). No code is written here; this makes the
> port mechanical. The **only genuine compute** is the strand-aware relative-coordinate projection inside
> `extract_and_pad_segments`; everything else (intersect / slop / sort+merge / coverage / freshness reuse)
> maps 1:1 onto already-built `sdrecall-io` + `sdrecall-utils` units and is **not** re-implemented here.

### 0. Exact input contract (must read before porting)

The producer of the per-RG "whole region" BED is **`preparation/build_beds_and_masked_genomes.py:174-189`**
(`establish_beds_per_RG_cluster`). It is the only thing that defines the column meaning T7 consumes.
The `all_homo_regions_bed` (a.k.a. `whole_region_bed`) is a **7-column, tab-separated, headerless** BED:

| col | 0 `chrom` | 1 `start` | 2 `end` | 3 `col4` | 4 `col5` | 5 `strand` | 6 `tag` |
|-----|-----------|-----------|---------|----------|----------|------------|---------|
| **FC** rows (build_beds:178/180) | chrom | abs start | abs end | `"."` | `"."` | `+`/`-` | `FC:{label}_{idx}` |
| **NFC** rows (build_beds:189) | chrom | abs start | abs end | `fc_node_rela_start` (int) | `fc_node_rela_end` (int) | `+`/`-` | `NFC:{label}_{idx}` |

For NFC rows, **col4/col5 are the NFC interval's coordinates *projected into the FC node's own coordinate
frame*** (computed by `HOMOSEQ_REGION.qnode_relative_region`, `preparation/homoseq_region.py:56-104` — that
projection is **T8's** job, not T7's; T7 only *reads* the already-written integers). This is exactly what
`extract_and_pad_segments` reads as `rel_start_interval = int(interval[3])`, `rel_end_interval = int(interval[4])`
(`prepare_masked_align_region.py:84-85`).

`pandas.read_table` is applied with explicit names
`["chrom","start","end","col4","col5","strand","tag"]` (`prepare_masked_align_region.py:162`), and tag-based
selection uses regex `^FC.*` (line 163), exact `FC:{rg}_{sub}` / `NFC:{rg}_{sub}` (lines 180-181). **Parity
note:** `rg_label` is *re-derived* from the filename — `os.path.basename(whole_region_bed).split("_")[0]`
(line 160) — overriding the passed-in `rg_label`. The Rust port must mirror this (label = first `_`-split
token of the BED basename), or feed the same value, to keep tag matching identical.

### 1. Python logic inventory

Source: `realign_recall/prepare_masked_align_region.py` (238 lines) + the format producer above.

| # | Python fn | lines | 1-line description | control-flow shape |
|---|-----------|-------|--------------------|--------------------|
| P1 | `extract_and_pad_segments` | 10-141 | **THE compute.** Intersect single FC interval with target → relative-coord pad+merge of overlapping segments → strand-aware slice of each NFC interval → sort+merge → coverage; with mtime-based output reuse. | 3 sequential loops: (a) per overlapping segment → relative coords (53-60); (b) sorted-merge of relative segments (65-76); (c) per NFC interval × per merged relative segment → strand-aware absolute projection (81-121). |
| P2 | `prepare_masked_align_region_per_RG_subgroup` | 200-236 | Per (RG, subgroup): split the 7-col frame into the single FC row + ≥1 NFC rows for that subgroup; write temp BEDs; call P1 for the **nfc** output; emit a CSV record string. | straight-line; 2 asserts (FC unique == 1, NFC ≥ 1). |
| P3 | `prepare_masked_align_region_per_RG` | 151-196 | Per RG: read whole BED once; build the **fc** target BED via `fc.intersect(target).slop(500).sort().merge()` + coverage; loop subgroups calling P2; string-substitute the shared fc path/size into each subgroup record. | 1 loop over `rg_subids`; FC target computed **once per RG** and shared. |
| P4 | `imap_prepare_masked_align_region_per_RG` | 146-147 | tuple-arg shim for `imap_unordered`. | trivial. |

Caller / parallel axis: `realign_and_recall.py:58-65` — `ctx.Pool(num_jobs).imap_unordered` over `uniq_rg_labels`
(one process per RG). In Rust this becomes **`rayon::par_iter` over RGs**, in-process (no pool/pickle).

#### P1 control flow, ported verbatim (the load-bearing arithmetic — cite, do not "improve")

1. **Overlap with target** (line 45): `overlapping_segments = single_interval.intersect(target_regions)` →
   `sdrecall-io::intersect(&[fc_iv], &target)`. `main_interval` = the *one* FC interval (line 48); capture
   `main_start, main_end, main_strand` (line 49).
2. **Relative padded segments** (53-60), per overlapping segment `seg`:
   `start = max(seg.start - padding, 0)`; `end = seg.end + padding` (**no upper clamp here** — clamping to the
   contig happens later via slop on the FC side and via interval-boundary clamps on the NFC side);
   `relative_start = start - main_start`; `relative_end = end - main_start`; push 4-tuple
   `(relative_start, relative_end, start, end)`. `padding` default **600** (line 14).
3. **Sort + merge of relative segments** (65-76): sort by **original absolute start `x[2]`**; fold:
   if `current[2] <= last[3]` (abs start ≤ last abs end) merge → `merged_end = max(last[3], current[3])`,
   `new = (last[0], merged_end - main_start, last[2], merged_end)`. **Parity quirk to preserve:** the merged
   relative-start is `last[0]` (carried from the *first* segment of the run, i.e. relative to `main_start`),
   and the merged relative-end is recomputed as `merged_end - main_start` (consistent with rel = abs − main_start).
4. **Strand-aware projection into each NFC interval** (81-121), per NFC `interval` with
   `rel_start_interval=int(col4)`, `rel_end_interval=int(col5)`:
   - overlap-filter merged relative segments with `[rel_start_interval, rel_end_interval)` using the predicate
     **`rel_start < rel_end_interval or rel_end > rel_start_interval`** (line 91 — note: an **OR**, see Risk R1)
     and clamp each to `(max(rs, rsi), min(re, rei))` (89-90).
   - per clamped `(rel_start, rel_end)` with `rel_start < rel_end` (98):
     - **same strand** (`interval_strand == main_strand`): `rel_small_start = rel_start - rel_start_interval`;
       `rel_small_end = rel_end - rel_start_interval` (103-104).
     - **opposite strand**: `rel_small_start = rel_end_interval - rel_end`;
       `rel_small_end = rel_end_interval - rel_start` (107-108) — the reverse-complement flip.
     - `abs_start = interval.start + rel_small_start`; `abs_end = interval.start + rel_small_end` (111-112);
       clamp to interval bounds `abs_start = max(interval.start, abs_start)`,
       `abs_end = min(interval.end, abs_end)` (115-116); emit `(chrom, abs_start, abs_end, strand=interval_strand, name=interval.name)`
       iff `abs_start < abs_end` (119-120).
5. **Finalize** (123-141): `BedTool(segments).sort().merge()` → `total_coverage()`; if `output_bed` given and is
   **fresher than all three inputs** (mtime, 126-129) and the existing coverage equals the freshly-computed
   coverage, **reuse** the existing file; else save and return `(path, coverage)`.

### 2. Python → Rust crate mapping

External deps replaced: **pybedtools/bedtools** (intersect/slop/sort/merge/total_coverage), **pandas**
(table read + tag filtering), **`os.path.getmtime`** (freshness). All in-process.

| Python operation / idiom | Rust crate::api | confidence |
|---|---|---|
| `pb.BedTool(path)` read 3-col / `pd.read_table(...names=7cols)` | `sdrecall-io::read_bed` (3-col target) + a thin **7-col tagged reader** in this crate (see §3 `AllRegionRow`) | verified-docs (sdrecall-io contract) |
| `single.intersect(target)` | `sdrecall-io::intersect(a, b)` → `Vec<GenomicInterval>` (wraps `bedrs::IntervalContainer` / `Intersect` trait) | verified-docs (bedrs re-exports `traits::Intersect`, `IntersectIter`; `Bed3::intersect(&b)->Option<Bed3>`) |
| `.slop(b=500, g=fai)` | `sdrecall-io::slop(ivs, 500, &chrom_sizes)` (clamps to `[0, contig_len]`) | verified-docs (sdrecall-io contract); bedrs `extend_left`/`extend_right` exist on `Coordinates` |
| `.sort().merge()` + `.total_coverage()` | `sdrecall-io::sort_merge_bed(ivs, stranded=false)` then `ivs.iter().map(GenomicInterval::len).sum::<i64>()` (replaces `total_coverage`) | verified-docs (bedrs re-exports `MergeIter`; `Coordinates::len()` confirmed) |
| `BedTool(segments).sort().merge()` for **stranded** nfc output | `sort_merge_bed(ivs, stranded=true)` — Python's final `.merge()` is **unstranded** (line 124/139); see Risk R2 | plausible |
| `df.loc[tag.str.contains("^FC.*")]` / `== "FC:{rg}_{sub}"` | filter `Vec<AllRegionRow>` by parsed `RgTag` enum (`Fc{label,sub}` / `Nfc{label,sub}`) — one parse, reused | verified-docs (std) |
| `int(interval[3])`, `int(interval[4])` | already-parsed `i64` fields on `AllRegionRow` (col4/col5) | verified-docs |
| `interval.strand` (`+`/`-`) vs `main_interval.strand` | `sdrecall-utils::Strand` (`Forward`/`Reverse`/`Unknown`); equality compare | verified-docs (interface contract) |
| `prepare_tmp_file(tmp_dir, ".bed")` | `sdrecall-io::prepare_tmp_file(&tmp_dir, ".bed")` | verified-docs (contract) |
| `os.path.getmtime` freshness + coverage-equality reuse | `sdrecall-utils::is_file_up_to_date(out, &[deps])` **+** explicit coverage compare (read existing → sum len) | verified-docs (contract); reuse needs coverage check too (see Risk R3) |
| `ref_genome.replace(".fasta",".fasta.fai")` + read `g` for slop | read `.fai` → `AHashMap<String,i64>` chrom_sizes (tiny helper in `sdrecall-io`, or local fai parse) | plausible |
| `ctx.Pool(...).imap_unordered` over RGs | `rayon` `par_iter()` over RG labels (in-process) | verified-docs (rayon) |
| CSV record string `f"{rg},{sub},...,{nfc},...,{size}"` | typed `RgSubgroupRecord` struct → `serde`/`write_tsv` at orchestrator, or `Display` for byte-compat (see Risk R4) | plausible |

### 3. Crate file layout (lib + bin/examples)

```
region-prep/
├─ Cargo.toml                    # deps: sdrecall-utils, sdrecall-io, rayon, thiserror, log, clap(bin); ahash, rustc-hash
├─ src/
│  ├─ lib.rs                     # pub re-exports; #![warn(clippy::all)]; pub use of the 4 units below
│  ├─ all_region_bed.rs          # AllRegionRow + RgTag parse; read_all_region_bed(); split_subgroup()  [I/O bookkeeping]
│  ├─ project.rs                 # THE unit: extract_and_pad_segments() — pure, no file I/O                [COMPUTE]
│  ├─ per_rg.rs                  # prepare_masked_align_region_per_rg(): FC target once + subgroup loop    [orchestration]
│  └─ bin/region_prep.rs         # thin CLI (clap): files in → fc/nfc BEDs out per RG
└─ examples/
   └─ diff_region_prep.rs        # differential harness: run on a dumped HG002/HG006 RG → assert vs Python BEDs
```

**One-versatile-unit-per-job mapping (no overlapping helpers — user rule #1):**

- `project::extract_and_pad_segments` — **the single compute unit.** Pure function over borrowed slices, *no
  file I/O, no temp files, no mtime logic*. Both the FC-side "intersect+slop+merge" and the NFC-side projection
  are distinct jobs; only the **NFC strand-aware projection** is genuine compute, so only *it* lives here. The FC
  target (`intersect → slop → sort_merge → coverage`) is **pure `sdrecall-io` calls** and stays in `per_rg.rs`
  (it is interval bookkeeping, explicitly out-of-scope per the task's scope boundary).
- **Collapse of Python's two strand branches:** lines 101-108 are two near-duplicate code paths differing only in
  how `(rel_small_start, rel_small_end)` is computed. In Rust they collapse into **one** match-on-`bool`
  (`same_strand`) returning a `(i64,i64)` tuple — *no* second helper fn. This is the only place the Python had
  duplicated logic; it becomes a single 4-line `match`.
- `all_region_bed::read_all_region_bed` + `split_subgroup` — the *only* tag-parsing/row-splitting unit; reused by
  both the FC-target path and every NFC subgroup (replaces the three separate `df.loc[...]` slices at
  prepare_masked_align_region.py:163/180/181 with one parsed-row filter).
- `per_rg::prepare_masked_align_region_per_rg` — pure **upstream orchestration**: reads the BED once, computes the
  shared FC target once, `rayon`-maps subgroups onto `extract_and_pad_segments`, writes outputs + records. It adds
  no compute of its own beyond sequencing existing units.

### 4. Core data structures + key fn signatures (explicit borrow/owner choices)

```rust
// --- all_region_bed.rs ---------------------------------------------------------
/// One parsed row of the 7-col all_homo_regions BED. Owns its chrom String once
/// (read once per RG); col4/col5 are pre-parsed i64 (sentinel for "." FC rows).
#[derive(Clone, Debug)]
pub struct AllRegionRow {
    pub chrom: String,        // owned: each row keeps its chrom; reused across the per-subgroup loop
    pub start: i64,
    pub end: i64,
    pub col4: i64,            // FC rows: i64::MIN sentinel ("."); NFC rows: rel_start_interval
    pub col5: i64,            // FC rows: i64::MIN sentinel ("."); NFC rows: rel_end_interval
    pub strand: Strand,       // sdrecall-utils::Strand
    pub tag: RgTag,           // parsed once (no repeated regex)
}
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum RgTag { Fc { label: String, sub: String }, Nfc { label: String, sub: String } }

/// Read the whole-region BED ONCE per RG. Owned Vec because every subgroup filters it.
pub fn read_all_region_bed(path: &Path) -> Result<Vec<AllRegionRow>>;
//  WHY owned Vec: read once, then borrowed (&[AllRegionRow]) by every subgroup — single alloc, many borrows.

/// Borrow the rows; return borrowed views (the one FC row + the NFC rows) for a subgroup.
pub fn split_subgroup<'a>(rows: &'a [AllRegionRow], label: &str, sub: &str)
    -> Result<(&'a AllRegionRow, Vec<&'a AllRegionRow>)>;
//  WHY &'a refs: no row data copied; asserts FC count==1 (Python:213) & NFC>=1 (Python:214) → SdError on violation.

// --- project.rs (THE compute unit) --------------------------------------------
/// Strand-aware relative-coordinate projection + pad + merge. PURE: no I/O, no mtime.
/// Ports extract_and_pad_segments (prepare_masked_align_region.py:10-121); the file-write +
/// freshness-reuse tail (123-141) stays in the orchestrator, not here.
pub fn extract_and_pad_segments(
    fc_interval:  &GenomicInterval,      // the single FC interval (borrow — used read-only)
    nfc_intervals:&[NfcInterval],        // borrowed slice — iterated once, not retained
    target:       &[GenomicInterval],    // borrowed — intersected against fc_interval
    padding:      i64,                   // default 600 (Python:14)
) -> Vec<GenomicInterval>;               // OWNED result: freshly-built nfc segments the caller writes/merges
//  WHY &GenomicInterval / &[..]: inputs are read-only; WHY owned Vec out: these are new intervals with no
//  borrow tie to inputs — caller owns + sorts/merges/writes them.

/// Minimal carrier for an NFC interval's projection inputs (avoids passing whole AllRegionRow into pure code).
pub struct NfcInterval<'a> {
    pub chrom: &'a str, pub start: i64, pub end: i64,
    pub strand: Strand, pub name: &'a str,        // Python keeps interval.name (line 120)
    pub rel_start_interval: i64, pub rel_end_interval: i64,  // col4/col5
}
//  WHY borrowed &str fields: chrom/name are reused verbatim from AllRegionRow; zero copy until the output
//  GenomicInterval is constructed.

// internal (project.rs) — the relative pad+merge, factored as ONE step (not a public helper):
//   fn relative_padded_merged(fc:&GenomicInterval, target:&[GenomicInterval], padding:i64)
//       -> Vec<(i64,i64)>;   // returns (relative_start, relative_end) merged segments only
//   WHY return only the 2 relative fields: the abs start/end (Python tuple slots 2,3) are needed solely as the
//   merge sort key & overlap test inside this fn — they don't escape, so don't widen the public type.

// --- per_rg.rs (orchestration) -------------------------------------------------
pub struct RgSubgroupRecord {
    pub rg_label: String, pub subgroup_id: String,
    pub fc_bed: PathBuf, pub nfc_bed: PathBuf,
    pub fc_bed_size: i64, pub nfc_bed_size: i64,
}
/// One RG: FC target computed ONCE, then rayon over subgroups. Mirrors prepare_masked_align_region_per_RG.
pub fn prepare_masked_align_region_per_rg(
    rg_label: &str,
    rg_subids: &[String],
    tmp_dir: &Path, target_region_bed: &Path, whole_region_bed: &Path,
    ref_genome: &Path,
) -> Result<Vec<RgSubgroupRecord>>;
//  WHY &str/&[String]/&Path borrows: caller (orchestrator) owns these; this fn only reads them.
//  WHY owned Vec<Record> out: records are produced fresh and handed up to be tabulated.
```

`fc_bed_size` / `nfc_bed_size` replace `total_coverage()` via `ivs.iter().map(GenomicInterval::len).sum()`.

### 5. Performance optimizations from ownership/borrowing

- **Read the 7-col BED once per RG, borrow everywhere.** `read_all_region_bed → Vec<AllRegionRow>` (one alloc);
  `split_subgroup` returns `&AllRegionRow` views per subgroup — Python re-`.loc[]`-filtered the DataFrame 3× per
  subgroup and re-wrote temp BEDs; we avoid all of that. No per-subgroup chrom/name `String` clones.
- **FC target computed once per RG, shared by all subgroups** (matches Python:171-175 hoist) — `intersect→slop→
  sort_merge→sum` runs a single time; its `PathBuf`/size are `Clone`-copied into each record (cheap), not
  recomputed.
- **`rayon::par_iter` over subgroups within an RG, and over RGs at the orchestrator** — the independent axis is
  (RG, subgroup); `extract_and_pad_segments` is pure and `Send`, so it parallelizes with zero shared mutable
  state. Replaces the Python `mp.Pool` (process fork + pickle + temp-file IPC) with in-process work-stealing.
- **`NfcInterval<'a>` borrows `chrom`/`name` `&str`** straight from `AllRegionRow`; the only allocation in the
  hot loop is the final `GenomicInterval` push (unavoidable — it's the output). Pre-`Vec::with_capacity` the
  output to `nfc_intervals.len()` (lower-bound) to cut reallocs.
- **`FxHashMap<String,i64>` for chrom_sizes** (string keys → per the project rule, `ahash`/`AHashMap`; int keys
  → `FxHashMap`). Here chrom→size is string-keyed → `AHashMap<String,i64>` (matches `slop`'s `chrom_sizes` param
  type in the T0 contract).
- **No temp BED round-trip for the FC/NFC intermediate.** Python wrote `fc_bed_path`/`nfc_bed_path` temp files
  purely to hand to bedtools; in-process we pass `&[GenomicInterval]` slices directly. Temp files are written
  only for the final `nfc` output (the durable artifact downstream stages read) and the shared `fc` target.
- **Streaming sort key without widening types:** the merge step keeps `(rel,rel)` public but sorts on the
  abs-start computed inline, so the abs coords never enter a returned struct (smaller cache footprint).

### 6. Risks / open decisions

- **R1 — the OR overlap predicate (line 91) is suspicious but must be ported verbatim.**
  `if rel_start < rel_end_interval or rel_end > rel_start_interval` is an **OR**, which is true for almost every
  segment (a real interval-overlap test would be the AND `rel_start < rei AND rel_end > rsi`). The subsequent
  `max/min` clamp (89-90) + the `rel_start >= rel_end → continue` guard (98) make most non-overlapping cases
  collapse to empty anyway, so the OR is *probably* harmless in practice — but **it is load-bearing for parity**.
  Port the OR exactly; flag for the user whether to "fix" it (would diverge from Python → must be a deliberate
  golden-encoding-style decision, validated against pileup/IGV, not silently).
- **R2 — final merge strandedness.** Python's closing `BedTool(corresponding_segments).sort().merge()`
  (line 124/139) is **unstranded** (default `merge`), even though each emitted interval carries a strand
  (line 120). So opposite-strand NFC segments that become adjacent **will be merged across strands**, dropping
  strand on the output. Port as `sort_merge_bed(ivs, stranded=false)` to match; note the strand column in the
  nfc BED output is therefore effectively collapsed/merged-away by Python — confirm downstream
  (Phase-2a read-extraction) does not depend on per-segment strand in the nfc BED.
- **R3 — freshness reuse needs both mtime AND coverage equality.** Python reuses the existing output only if it's
  newer than all 3 inputs **and** `ori_bed_size == return_bed_size` (lines 126-134). `is_file_up_to_date` covers
  only the mtime half; the coverage half must be added explicitly (read existing nfc BED → sum len → compare),
  otherwise we'd reuse a stale-but-newer file. Decision: keep both checks for byte-parity, or (cleaner) drop the
  coverage check and rely on content hashing via `update_file_on_md5` — needs user sign-off since it changes the
  reuse condition.
- **R4 — output record format.** Python returns a comma-joined string
  `"{rg},{sub},target_fc_region_bed,{nfc},target_fc_region_size,{nfc_size}"` with literal placeholder tokens that
  `prepare_masked_align_region_per_RG` then string-substitutes (lines 191-192). The caller parses these with
  `pd.read_csv(names=[...])` (`realign_and_recall.py:84`). In the pure-Rust world this becomes a typed
  `RgSubgroupRecord` consumed in-process. The production orchestrator now uses the typed record directly; exact
  6-field CSV formatting is retained only where a historical differential/CLI fixture requires it. The old
  PyO3 placeholder-substitution boundary is gone.
- **R5 — `rg_label` re-derivation from filename** (line 160) overrides the argument; if the Rust orchestrator
  passes a label that disagrees with `basename.split("_")[0]`, tag matching silently yields 0 FC rows → the
  `len(fc)==1` assert fires. Mirror the filename derivation, or assert the two agree.
- **R6 — `padding=600` vs `slop b=500` are two different constants** (P1 default 600 on the NFC side; P3 uses
  500 on the FC target side). Do not unify them; keep as two named consts (`NFC_PAD=600`, `FC_SLOP=500`).
- **R7 — `.fai` dependency for slop.** `slop` needs contig sizes from `<ref>.fasta.fai`
  (`prepare_masked_align_region.py:172`). Per the dependency-availability rule: **require** the `.fai`, fail
  clearly (`SdError::Io`) if absent — no fallback that skips boundary clamping.
- **R8 — empty-input edge.** Python indexes `relative_segments[0]` (line 66) with no guard: if the FC interval
  doesn't intersect the target, `overlapping_segments` is empty → `IndexError`. The Rust port should return an
  empty `Vec` (the correct, non-crashing behavior) — a **deliberate, documented** divergence (more correct), not
  a silent fallback; flag for user confirmation.

### 7. Test plan delta (refines the doc's tier-1 / tier-2 above)

**Tier 1 — unit (in `project.rs`, hand-built fixtures, exact assertions):**
- `proj_same_strand`: FC `chr1:1000-2000 +`; target `chr1:1200-1400`; padding 100 → relative segment
  `[1100-100, 1500-1000]` → for an NFC `chr1:5000-6000 +` with col4/col5 `0..1000`, assert emitted
  `chr1:5100-5500 +` (same-strand: `abs = nfc.start + (rel - rel_start_interval)`).
- `proj_opposite_strand`: same FC/target; NFC `chr1:5000-6000 -`, col4/col5 `0..1000` → assert the
  reverse-complement flip (`rel_small_start = rei - rel_end`, `rel_small_end = rei - rel_start`) yields the
  mirrored coordinates; pin exact numbers from the Python arithmetic (lines 107-108, 111-116).
- `pad_clamp_at_zero`: overlapping segment starting < `padding` from contig 0 → `start = max(seg.start-pad, 0)`
  (line 55) clamps to 0; assert relative_start becomes `-main_start` not negative-of-pad.
- `merge_adjacent_relative`: two target overlaps whose padded abs ranges touch (`current[2] <= last[3]`) merge
  into one with `merged_end = max(...)` and relative-start carried from the first (the line-72 quirk).
- `overlap_predicate_OR`: a merged relative segment fully outside `[rsi,rei)` but satisfying the OR (line 91) →
  assert the clamp+`rel_start>=rel_end` guard drops it to empty (locks in R1 behavior).
- `interval_bound_clamp`: projection that would exceed `interval.end` → assert `abs_end = min(interval.end,..)`
  (line 116) and the `abs_start<abs_end` emit guard (line 119).
- `empty_no_intersect`: FC disjoint from target → assert `extract_and_pad_segments` returns `vec![]` (R8).

**Tier 2 — differential vs Python (`examples/diff_region_prep.rs`):**
- Inputs: dump a real `all_homo_regions_bed` (7-col) + `target_bed` + `ref.fai` from a prepared **HG002 + HG006**
  run (the prompt/task says per-RG BEDs from a prepared run). For each `(rg_label, subgroup_id)`:
  1. run Python `prepare_masked_align_region_per_RG_subgroup` → reference nfc BED + size + the shared fc BED + size;
  2. run Rust `per_rg::prepare_masked_align_region_per_rg` → same artifacts;
  3. assert **sorted-merged interval set-equality** on fc and nfc BEDs (the task's pass criterion) AND
     `total_coverage`/`sum(len)` equality.
- Run with `RUST_LOG=region_prep=debug`, redirect to `/paedyl01/disk1/yangyxt/test_tmp/diff_region_prep_<date>.log`
  per the project's test-verification requirement; report the log path.
- Parity-hazard assertions to surface early: (a) any RG where the OR predicate (R1) makes Rust≠Python; (b) any
  nfc BED where stranded-vs-unstranded final merge (R2) diverges; (c) freshness-reuse path (R3) — run twice, assert
  second run reuses iff coverage+mtime both hold.
