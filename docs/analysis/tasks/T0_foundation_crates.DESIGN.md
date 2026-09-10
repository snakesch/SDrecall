> Historical companion design appendix for `T0_foundation_crates.md`, recorded before implementation. T0 is production-integrated as of 2026-07-17. Preserve this file for interface rationale; the live source and the current task status document implementation deviations such as orchestrator-local `Paths` and allowed external-tool wrappers.

# T0 DESIGN — `sdrecall-utils` + `sdrecall-io` (coding-grade)

Scope reminder: T0 defines the **inter-crate interface contract** that T5–T9 consume. This appendix records the exact Rust type proposals, function signatures with explicit borrow/owner choices, file layout, duplicate-collapse map (DUP-1/DUP-2), verified crate API mapping, and original test plan. The implementation has since landed; proposed signatures here are historical unless they match the live API.

Verified crate APIs (docs.rs, 2026-06-11, via SOCKS proxy):
- rust-htslib 0.47.0 — `IndexedReader::from_path<P: AsRef<Path>>(P)->Result<Self>`; `Read` trait: `read(&mut Record)->Option<Result<()>>`, `records()->Records`, `rc_records()->RcRecords`, `set_threads(usize)->Result<()>`; `IndexedReader::fetch<T: Into<FetchDefinition>>(T)`; `Writer::from_path(P, header:&Header, format:Format)->Result<Self>`, `Writer::set_threads`, `write(&Record)`.
- bedrs 0.2.26 — types `Bed3`, `StrandedBed3`, `IntervalContainer`, `Strand`; traits `Coordinates`, `Intersect`, `Overlap`, `Subtract`, `Merge`; `StrandedBed3::try_from(record)`.
- petgraph-graphml 5.0.0 — `GraphMl::new(G)`, `.pretty_print(bool)`, `.export_node_weights(Box<dyn Fn(&NW)->Vec<(Cow<str>,Cow<str>)>>)`, `.export_edge_weights(...)`, `.to_writer<W: Write>(&self, W)->Result<()>`.

---

## 1. Python logic inventory

### `sdrecall-utils` origins

| Python fn | File:line | 1-line description | control-flow shape |
|-----------|-----------|--------------------|--------------------|
| `SDrecallPaths.__init__` | `src/const.py:69-127` | derive `assembly`, `sample_id`, `target_tag`, `work_dir`, frag-size stats; mkdir tree | one-shot derivation; calls `get_insert_size_distribution` (I/O — moves to io crate) |
| `_extract_assembly_version` | `src/const.py:130-164` | substring-match hg19/hg38/chm13 from sd-map then ref name | cascade of `if substr in name`; raises on miss |
| `_extract_sample_id` | `src/const.py:166-170` | basename split on `.`, take `[0]` | pure string |
| `_extract_target_tag` | `src/const.py:172-181` | "exome" default / `default_target` → "exome" / basename `[0]` | pure string |
| `_normalize_rg_label` | `src/const.py:301-330` | int/numeric-str/`RG\d+` → `RG{n}`; else raise | match + regex; **pure** (drop the side-effecting mkdir from path getters) |
| all `*_path()` getters | `src/const.py:230-474` | ~35 derived path strings from `work_dir`/`recall_results_dir`/`rg_dir`/`repo_dir` + `sample_id`/`basename`/`assembly`/`target_tag` | pure `os.path.join`; **two side-effecting ones** (`ref_genome_fai_path` shells `samtools faidx`, `rg_dir` mkdirs) split out |
| `register_realign_group` / `_discover_and_register_realign_groups` | `src/const.py:266-298` | track `RG{n}` set; discover from FS by `^RG\d+$` dir scan | set insert + `read_dir` loop |
| `check_*_validity` | `src/const.py:476-665` | mtime freshness + size + BED-field/VCF-header sanity | freshness loops; shells `bash shell_utils.sh`/`bcftools` (→ io crate) |
| `configure_parallelism` | `src/utils.py:158-167` | `ceil(total/per_job)` → `(num_jobs, per_job)` | pure arithmetic |
| `is_file_up_to_date` | `src/utils.py:39-41` | `mtime(target) > all mtime(deps)` | pure-ish (stat) |
| `na_value` | `src/utils.py:44-66` | 4 regexes for NA-like strings + None/NaN | match + regex |
| `prepare_tmp_file` | `src/utils.py:34-36` | mkdir + `NamedTemporaryFile(delete=False)` | one-shot (FS) |
| `update_plain_file_on_md5` | `src/utils.py:69-91` | md5 old vs new; replace-or-delete | hash compare (FS) |
| logging | `src/log.py:11-134` | `ColoredFormatter`, `init_logger`, `configure_logger`, `log_command`/`error_handling_decorator`/`log_decorator` → `(ok, result, log)` tuples for `imap_*` workers | decorator wrap; per-worker `StringIO` capture |
| warning silencing | `src/suppress_warning.py:1-9` | mute Numba/NumPy/Future warnings | no Rust analog needed (Numba/NumPy gone); record as N/A |

### `sdrecall-io` origins

| Python fn | File:line | description | control-flow |
|-----------|-----------|-------------|--------------|
| `migrate_bam_to_ncls` | `fp_control/bam_ncls.py:288-465` | **canonical** BAM→interval-index: collate-by-qname, drop noisy reads, build per-chrom interval tree + read_dict/qname maps | single streaming pass over collated records, group-by-qname flush |
| `is_read_noisy` | `fp_control/bam_ncls.py:115-193` | **canonical** noisy-read predicate (primary-only; unmapped/qcfail/MAPQ/ref_end/seq; paired: span<75, mate-tid, proper-pair warn; SE: dup; if `filter_noisy`: median-baseQ, #Q<cut≥75, softclip≥75) | linear checks |
| `_collate_bam_file` | `fp_control/bam_ncls.py:196-216` | `samtools collate -f` to temp BAM (paired only) | subprocess |
| `_process_qname_group` | `fp_control/bam_ncls.py:219-285` | per-qname: noisy→skip; paired needs R1+R2; assign qname_idx; merge interval per chrom | grouping |
| `overlapping_reads_iterator` / `overlap_qname_idx_iterator` | `fp_control/bam_ncls.py:13-89` | overlap query → reads / qname_idx (dedup, true-overlap re-check) | generator |
| `get_insert_size_distribution` | `src/insert_size.py:10-49` | `num_runs`×sample reads (flag `0x2 && !0x90C`, same ref, `|TLEN|≤8000`, prob 0.001), 99th-pct trim, mean/median/std, average across runs | nested loop, early-break at `num_samples` |
| `merge_bams` | `src/utils.py:170-205` | header-fix + validity + `samtools merge | sort` + index | subprocess |
| `combine_vcfs` | `src/utils.py:137-155` | `bcftools concat -a -d exact | sort -Oz` + tabix | subprocess |
| `sortBed_and_merge` | `src/utils.py:95-115` | `bed.sort().merge(s=True, c=…, o=first/distinct)` | pybedtools |
| `merge_bed_files` | `src/utils.py:117-133` | dedup paths, `cat(postmerge=False)` then `sort` | pybedtools |

DUP-1 (three Rust copies to collapse): `build_phasing_graph/src/bam_reading.rs:16-365`, `haplotype_inspection/src/bam_lappers.rs:39-618`, `read_extraction/src/lib.rs:42-103`.

---

## 2. Python → Rust crate mapping

| Python operation / idiom | Rust crate::api | confidence |
|--------------------------|-----------------|------------|
| `pysam.AlignmentFile(p,"rb")` indexed `.fetch(chrom,s,e)` | `rust_htslib::bam::IndexedReader::from_path(p)?` + `.fetch((chrom, s, e))` (`T: Into<FetchDefinition>`) | verified-docs |
| `pysam` stream all records (collated, no index) | `rust_htslib::bam::Reader::from_path(p)?` + `Read::records()` / `Read::rc_records()` | verified-docs |
| `read.mapping_quality / is_unmapped / is_secondary / reference_start / reference_end / cigartuples / query_qualities / next_reference_name` | `bam::Record::{mapq,is_unmapped,is_secondary,pos,cigar().end_pos(),qual,cigar(),mtid}` + `HeaderView::tid2name` (already used in `bam_lappers.rs`) | verified-docs |
| write BAM | `bam::Writer::from_path(p, &header, Format::Bam)?` + `set_threads` + `write(&rec)` | verified-docs |
| `samtools merge | sort | index` (`merge_bams`) | rust-htslib: copy records from each input `Reader` into one `Writer`; then `bam::index::build`. **Header SQ-line reconcile is non-trivial — see Risks.** | plausible |
| `bcftools concat -a -d exact | sort` (`combine_vcfs`) | `bam::Reader`-style `bcf::Reader`/`bcf::Writer` (`rust_htslib::bcf`) for read+write; **sort + exact-dedup implemented in Rust** (k-way merge of coordinate-sorted inputs) | plausible |
| `pybedtools sort().merge(s=True,…)` (`sortBed_and_merge`) | `bedrs::IntervalContainer` `.sort()` + `Merge` trait `.merge()` (stranded via `StrandedBed3`) | verified-docs |
| `pybedtools cat(postmerge=False).sort()` (`merge_bed_files`) | concat `Vec<Bed3>` + `IntervalContainer::new(...).sort()` | verified-docs |
| `bedtools intersect / slop / complement` (used in region-prep/sd-prep) | bedrs `Intersect`/`Subtract`/`complement` + manual grow (`start-=n; end+=n` clamped to chrom size) | verified-docs (intersect/subtract); plausible (slop) |
| `NCLS(starts,ends,idx)` + `all_overlaps_both` | `rust_lapper::Lapper<u32,u32>::new(intervals)` + `.find(s,e)` (already in `bam_lappers.rs:663`) | verified-docs |
| graph-tool `.save(graphml)` (`prepare_recall_regions.py`/`graph_query.py`) | `petgraph_graphml::GraphMl::new(&g).export_node_weights(...).export_edge_weights(...).to_writer(w)?` | verified-docs |
| graph-tool `.load(graphml)` (only if a stage reads back) | `quick_xml::Reader` event parse → `petgraph::Graph` (no graphml *reader* in petgraph-graphml) | plausible |
| `np.mean/median/std` on insert sizes; `np.percentile(...,99)` | `statrs` / hand-rolled: mean = sum/n, std = population (matches `np.std` ddof=0), median = mid-of-sorted, 99th pct = `np.percentile` *linear interpolation* (must replicate, see Risks) | plausible |
| `tempfile.NamedTemporaryFile(delete=False)` | `tempfile::Builder::new().suffix(..).tempfile_in(dir)?` (or `.keep()` for delete=False semantics) | verified-docs |
| `hashlib.md5(file)` compare | `md-5` crate (`Md5::digest`) or `seahash` for non-crypto; compare hex | plausible |
| `os.path.getmtime` freshness | `std::fs::metadata(p)?.modified()?` (`SystemTime`) | verified-docs |
| `logging` colored + per-worker capture tuples | `log` 0.4 facade + a custom `Logger` impl (colored to stderr) + a buffer-capturing sink for the `(ok,result,log)` worker pattern | plausible |

---

## 3. Crate file layout

```
sdrecall-utils/                    (lib only; deps: serde, thiserror, log, rustc-hash, ahash)
├─ src/lib.rs                       # re-export the public vocabulary; #![warn(clippy::all)]
├─ src/interval.rs                  # GenomicInterval, Strand, RegionKey + overlap/contain/merge (the ONE interval value-type; bedrs is the engine in io)
├─ src/ids.rs                       # HapId(i32), QnameIdx(u32) newtypes (+ From/Display)
├─ src/error.rs                     # SdError (thiserror) — crate-wide error enum
├─ src/paths.rs                     # Paths (port of SDrecallPaths): pure derivation + getters (NO I/O)
├─ src/parallel.rs                  # configure_parallelism (the one budget unit)
├─ src/freshness.rs                 # is_file_up_to_date, na_value (the one NA classifier)
└─ src/logging.rs                   # init_logger/configure_logger + WorkerLog capture (the one logging unit)

sdrecall-io/                       (lib + optional bin; deps: sdrecall-utils, rust-htslib, bedrs, petgraph(+graphml), quick-xml)
├─ src/lib.rs                       # re-export readers/writers; #![warn(clippy::all)]
├─ src/bam_read.rs                  # ★ DUP-1 home: open→collate→drop-noisy→interval-index (the ONE reader) + is_read_noisy (the ONE predicate)
├─ src/bam_write.rs                 # write_bam / sort / index / merge (the ONE bam writer unit)
├─ src/bed.rs                       # read_bed/write_bed + merge/intersect/slop/complement (bedrs wrappers, the ONE bed unit)
├─ src/vcf.rs                       # read_vcf (sorted cursor) / write_vcf / concat / sort (the ONE vcf unit)
├─ src/graphml.rs                   # write_graphml (petgraph-graphml) + read_graphml (quick-xml) 
├─ src/tsv.rs                       # read_tsv / write_tsv (the ONE tsv unit; serde rows)
├─ src/insert_size.rs               # get_insert_size_distribution (the ONE frag-stat unit)
├─ src/tmp.rs                       # prepare_tmp_file + update_on_md5 (the ONE atomic-replace unit)
└─ examples/reemit_bed_vcf.rs       # differential harness: re-emit HG002 fixtures, diff vs Python (links the lib)
```

**Duplicate-collapse decisions (enforcing the #1 coding rule):**
- **DUP-1 → `bam_read.rs::IndexReadBam` (one reader).** The three current copies differ only in (a) output container (FASTQ pairs vs `Lapper` vs raw records) and (b) whether they build an index. Collapse by making the reader *yield filtered qname-groups* through one callback/iterator; the **output shape is the caller's job**, not a second reader. Signature in §4 (`for_each_qname_group`). `read_extraction`'s lighter mate-fetch path becomes a thin caller of the same group iterator. The `is_read_noisy` predicate is the canonical `bam_ncls.py:115-193` one (single source).
- **DUP-3 → one interval backend.** `sdrecall-utils::GenomicInterval` is the *value type* (cheap, serde, hashable); `bedrs::IntervalContainer` is the *set-ops engine* inside io; `rust_lapper::Lapper` is the *point-query index* inside `bam_read.rs`. Three roles, no overlap — the hand-rolled `structs.rs:115` sweep and `graph_builder.rs` endpoint sweeps are retired in favor of these.
- **HYG-6:** `bam_read.rs` reaps the collate child with a plain `drop`/`child.wait()`, **never** `mem::forget` (the `bam_lappers.rs:103` fd-leak); htslib already dup'd the fd. Add the pipe-read + fd-count regression test.
- **DUP-2 NOTE (cross-task):** the per-read **vector toolkit** (extract_hap/error_vector, base encoder, var counting) is *compute*, not I/O — it does **not** live here. Its single home is `haplotype_inspection` today and `fp-control` (T4) after fusion. T0 only defines the **seam**: `bam_read.rs` hands out `&Record` (or owned `Record`) groups; the toolkit borrows from them. The M-CIGAR decision is recorded as a flag (§6), applied at T4, not T0.

---

## 4. Core data structures + key fn signatures

### `sdrecall-utils` types (THE interface contract)

```rust
// ── ids.rs ─────────────────────────────────────────────────────────────
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, serde::Serialize, serde::Deserialize)]
pub struct HapId(pub i32);     // WHY i32: matches Python hap labels incl. -1 sentinel; Copy → pass by value
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, serde::Serialize, serde::Deserialize)]
pub struct QnameIdx(pub u32);  // WHY u32: dense 0..N interval-index from bam_read; Copy, fits Lapper<_,u32> val

// ── interval.rs ────────────────────────────────────────────────────────
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, serde::Serialize, serde::Deserialize)]
pub enum Strand { #[default] Unknown, Forward, Reverse }   // WHY Copy enum: 1 byte, no alloc

#[derive(Clone, PartialEq, Eq, Hash, Debug, serde::Serialize, serde::Deserialize)]
pub struct GenomicInterval {
    pub chrom: String,   // WHY owned String: intervals outlive any single BAM header; cross-crate move-friendly
    pub start: i64,      // WHY i64: htslib pos is i64; half-open BED [start,end)
    pub end: i64,
    pub strand: Strand,
}
impl GenomicInterval {
    pub fn overlaps(&self, other: &Self) -> bool;          // &self,&other: read-only, no alloc
    pub fn contains_point(&self, chrom: &str, pos: i64) -> bool;
    pub fn len(&self) -> i64;                              // end-start
    pub fn region_key(&self) -> RegionKey;                 // borrow→owned key, see below
}

/// Hashable map key for per-region dicts (the `(chrom,start,end)` Python tuple key).
#[derive(Clone, PartialEq, Eq, Hash, Debug, serde::Serialize, serde::Deserialize)]
pub struct RegionKey { pub chrom: String, pub start: i64, pub end: i64 }
// WHY a distinct type (not reuse GenomicInterval): strand-insensitive key; smaller; intent-revealing.

// ── error.rs ───────────────────────────────────────────────────────────
#[derive(thiserror::Error, Debug)]
pub enum SdError {
    #[error("I/O error on {path}: {source}")]
    Io { path: String, #[source] source: std::io::Error },
    #[error("htslib error: {0}")]
    Htslib(String),                                        // WHY String: rust_htslib::errors::Error not Clone/Send-clean across boundary
    #[error("BAM CIGAR contains 'M' op but =/X (--eqx) mode is required for {qname}")]
    CigarMOp { qname: String },                            // the DUP-2a seam, surfaced as a typed error not a panic
    #[error("read name is not valid UTF-8")]
    NonUtf8ReadName,
    #[error("could not determine assembly from sd_map={sd_map} or ref={reference}")]
    UnknownAssembly { sd_map: String, reference: String },
    #[error("invalid RG label or index: {0}")]
    InvalidRgLabel(String),
    #[error("BED parse error at line {line}: {msg}")]
    BedParse { line: usize, msg: String },
    #[error("VCF/BCF error: {0}")]
    Vcf(String),
    #[error("GraphML error: {0}")]
    GraphMl(String),
    #[error("insufficient paired-end data ({0} qnames)")]
    InsufficientPairs(usize),
}
pub type Result<T> = std::result::Result<T, SdError>;       // crate-wide alias

// ── paths.rs ───────────────────────────────────────────────────────────
#[derive(Clone, Debug)]
pub struct Paths {
    pub ref_genome: PathBuf, pub input_bam: PathBuf, pub reference_sd_map: PathBuf,
    pub target_bed: Option<PathBuf>, pub output_dir: PathBuf, pub repo_dir: PathBuf,
    pub assembly: String, pub sample_id: String, pub target_tag: String,
    pub basename: String, pub work_dir: PathBuf,
    pub recall_results_dir: PathBuf, pub realign_groups_dir: PathBuf, pub tmp_dir: PathBuf,
    // frag-size stats are computed in sdrecall-io and SET on Paths (keeps utils I/O-free):
    pub avg_frag_size: Option<f64>, pub median_frag_size: Option<f64>, pub frag_size_std: Option<f64>,
}
impl Paths {
    /// Pure derivation — NO mkdir, NO samtools, NO insert-size (caller wires those from io).
    pub fn derive(
        ref_genome: &Path, input_bam: &Path, reference_sd_map: &Path, output_dir: &Path,
        target_bed: Option<&Path>, sample_id: Option<&str>, target_tag: Option<&str>,
        ref_genome_tag: Option<&str>, repo_dir: &Path,
    ) -> Result<Self>;                                      // borrows &Path, returns owned Paths
    // every *_path() getter returns an OWNED PathBuf (cheap join; callers store/compare):
    pub fn pooled_raw_bam_path(&self) -> PathBuf;
    pub fn recall_raw_vcf_path(&self) -> PathBuf;
    pub fn final_recall_vcf_path(&self) -> PathBuf;
    pub fn realign_meta_table_path(&self) -> PathBuf;
    pub fn qnode_grouping_graph(&self) -> PathBuf;
    pub fn multi_align_bed_path(&self) -> PathBuf;
    pub fn multiplex_graph_path(&self) -> PathBuf;
    pub fn annotated_graph_path(&self) -> PathBuf;          // .replace(".graphml",".trim.annoPC.graphml")
    pub fn directed_graph_path(&self, chrom: &str) -> PathBuf;
    // RG getters take the same flexible label and normalize internally:
    pub fn rg_query_bed_path(&self, rg: RgRef<'_>) -> PathBuf;
    pub fn masked_genome_path(&self, rg: RgRef<'_>) -> PathBuf;
    pub fn minimap_index_path(&self, rg: RgRef<'_>) -> PathBuf;
    pub fn rg_dir_path(&self, rg: RgRef<'_>) -> PathBuf;    // PURE join (mkdir is an io-crate concern)
    // ... full 1:1 set with src/const.py:230-474 ...
    pub fn normalize_rg_label(rg: RgRef<'_>) -> Result<String>;  // RG{n}, the one normalizer
}
/// One input accepted everywhere an RG is named — collapses the int/str/RG-str Python overloads.
pub enum RgRef<'a> { Index(u32), Label(&'a str) }          // WHY borrow &str: no alloc at call sites

// ── parallel.rs / freshness.rs ─────────────────────────────────────────
pub fn configure_parallelism(total_threads: usize, threads_per_job: f64) -> (usize, usize);
// WHY (usize,usize) owned tuple: tiny Copy values, matches Python return.
pub fn is_file_up_to_date(target: &Path, deps: &[&Path]) -> std::io::Result<bool>; // &[&Path]: borrow list, no Vec move
pub fn na_value(s: &str) -> bool;                          // &str: read-only classify, the one NA unit

// ── logging.rs ─────────────────────────────────────────────────────────
pub fn init_console_logger(level: log::LevelFilter);      // colored stderr (ColoredFormatter port)
/// Per-worker capture replacing log_command/error_handling_decorator's (ok,result,log) tuple.
pub fn run_captured<T, E: std::fmt::Display>(
    f: impl FnOnce() -> std::result::Result<T, E>,
) -> (bool, std::result::Result<T, String>, String);
// WHY FnOnce + owned return: one-shot worker; the captured log String is owned (sent to imap_* collector).
```

### `sdrecall-io` API surface

```rust
use sdrecall_utils::{GenomicInterval, QnameIdx, Result};
use rust_htslib::bam::Record;

// ── bam_read.rs — the ONE reader (DUP-1) ───────────────────────────────
pub struct NoisyFilter { pub mapq: u8, pub basequal_median: u8, pub paired: bool, pub filter_noisy: bool }
impl Default for NoisyFilter { /* mapq=10, basequal_median=15, paired=true, filter_noisy=true */ }

/// THE consolidated noisy predicate (canonical = bam_ncls.py:115-193). Borrows the record.
pub fn is_read_noisy(rec: &Record, hv: &rust_htslib::bam::HeaderView, f: &NoisyFilter) -> bool;

/// Streaming primitive: open (collate if paired) → drop secondary/supp/dup → group by qname
/// → drop noisy groups → hand each clean group to `sink`. ONE function; callers decide the shape
/// (FASTQ pairs, Lapper index, raw record vec). Replaces all 3 DUP-1 copies.
pub fn for_each_qname_group<F>(bam_path: &Path, f: &NoisyFilter, threads: u8, mut sink: F) -> Result<()>
where F: FnMut(&str /*qname*/, &[Record] /*clean group, borrowed*/) -> Result<()>;
// WHY &[Record] borrow: sink reads/clones selectively; reader keeps ownership → no per-read clone in the common path.

/// Convenience built on `for_each_qname_group`: builds the per-chrom point-query index
/// (the haplotype_inspection use-case). Returns owned maps (callers store them for the island's lifetime).
pub struct BamIndex {
    pub lapper: rustc_hash::FxHashMap<String, rust_lapper::Lapper<u32, QnameIdx>>,
    pub reads:  rustc_hash::FxHashMap<QnameIdx, Vec<Record>>,   // owned: index outlives the reader
    pub qname:  rustc_hash::FxHashMap<QnameIdx, String>,
    pub qname_idx: ahash::AHashMap<String, QnameIdx>,
    pub noisy:  ahash::AHashSet<String>,
}
pub fn build_bam_index(bam_path: &Path, f: &NoisyFilter, threads: u8) -> Result<BamIndex>;

/// Indexed region fetch for callers that don't need the full in-memory index.
pub struct RegionReader { /* wraps bam::IndexedReader */ }
impl RegionReader {
    pub fn open(bam_path: &Path, threads: u8) -> Result<Self>;
    pub fn fetch(&mut self, iv: &GenomicInterval) -> Result<()>;   // &GenomicInterval: borrow, no move
    pub fn records(&mut self) -> impl Iterator<Item = Result<Record>> + '_;  // owned Record per yield (htslib semantics)
}

// ── bam_write.rs ───────────────────────────────────────────────────────
pub fn write_bam(out: &Path, header: &rust_htslib::bam::Header, recs: impl IntoIterator<Item = Record>, threads: u8) -> Result<()>;
pub fn sort_bam(input: &Path, out: &Path, threads: u8) -> Result<()>;     // in-Rust k-way/merge sort, no samtools
pub fn index_bam(bam: &Path) -> Result<()>;                              // rust_htslib::bam::index::build
pub fn merge_bams(inputs: &[&Path], out: &Path, ref_fasta: &Path, threads: u8) -> Result<()>;
// WHY &[&Path] inputs: borrow the list; merge copies records into one Writer (replaces merge_bams shell).

// ── bed.rs ─────────────────────────────────────────────────────────────
pub fn read_bed(path: &Path) -> Result<Vec<GenomicInterval>>;            // owned Vec: caller mutates/moves
pub fn write_bed(path: &Path, ivs: &[GenomicInterval]) -> Result<()>;    // borrow slice, read-only emit
pub fn sort_merge_bed(ivs: &[GenomicInterval], stranded: bool) -> Vec<GenomicInterval>; // sortBed_and_merge
pub fn merge_bed_files(paths: &[&Path]) -> Result<Vec<GenomicInterval>>; // dedup paths + cat + sort
pub fn intersect(a: &[GenomicInterval], b: &[GenomicInterval]) -> Vec<GenomicInterval>;
pub fn complement(ivs: &[GenomicInterval], chrom_sizes: &ahash::AHashMap<String,i64>) -> Vec<GenomicInterval>;
pub fn slop(ivs: &[GenomicInterval], by: i64, chrom_sizes: &ahash::AHashMap<String,i64>) -> Vec<GenomicInterval>;
// WHY borrow-in / owned-out everywhere: bedrs consumes an IntervalContainer internally; callers chain freely.

// ── vcf.rs ─────────────────────────────────────────────────────────────
pub fn read_vcf(path: &Path) -> Result<rust_htslib::bcf::Reader>;        // sorted cursor; caller streams
pub fn write_vcf(path: &Path, header: &rust_htslib::bcf::header::HeaderView, recs: impl IntoIterator<Item = rust_htslib::bcf::Record>) -> Result<()>;
pub fn concat_sort_vcfs(inputs: &[&Path], out: &Path, dedup_exact: bool, threads: u8) -> Result<()>; // combine_vcfs

// ── graphml.rs ─────────────────────────────────────────────────────────
pub fn write_graphml<N, E>(path: &Path, g: &petgraph::Graph<N, E>,
    node_attrs: impl Fn(&N) -> Vec<(String, String)>,
    edge_attrs: impl Fn(&E) -> Vec<(String, String)>) -> Result<()>;     // wraps GraphMl::new(&g)...to_writer
pub fn read_graphml(path: &Path) -> Result<petgraph::Graph<ahash::AHashMap<String,String>, ahash::AHashMap<String,String>>>; // quick-xml

// ── tsv.rs ─────────────────────────────────────────────────────────────
pub fn read_tsv<T: serde::de::DeserializeOwned>(path: &Path) -> Result<Vec<T>>;
pub fn write_tsv<T: serde::Serialize>(path: &Path, rows: &[T]) -> Result<()>;

// ── insert_size.rs ─────────────────────────────────────────────────────
pub struct FragStats { pub mean: f64, pub median: f64, pub std: f64 }
pub fn get_insert_size_distribution(bam: &Path) -> Result<Option<FragStats>>; // None == Python (None,None,None)

// ── tmp.rs ─────────────────────────────────────────────────────────────
pub fn prepare_tmp_file(tmp_dir: &Path, suffix: &str) -> Result<tempfile::NamedTempFile>;
pub fn update_file_on_md5(old: &Path, new: &Path) -> Result<bool>;       // replace-if-different; bool = replaced
```

---

## 5. Performance optimizations from ownership/borrowing

1. **Zero per-record clone in the hot read path.** `for_each_qname_group` hands the sink `&[Record]`; the FASTQ/region callers read fields by reference and only `clone()` the rare records they keep. The current `read_extraction` path clones every mate (`lib.rs:201-204`) and `bam_lappers.rs:519` clones the whole group — the consolidated reader clones only when the caller (`build_bam_index`) genuinely needs owned records for the island lifetime.
2. **`rc_records()` for the index build.** htslib's `RcRecords` reuses one allocation per record; `build_bam_index` uses it and clones into the owned `reads` map exactly once per retained read.
3. **`FxHashMap` for the `QnameIdx`/`HapId` integer-keyed maps** (`reads`, `lapper`) and **`ahash` for the `String`-keyed maps** (`qname_idx`, `noisy`) — the project convention, and measured-faster than SipHash for these key shapes.
4. **Borrow-in / owned-out BED & interval ops.** `&[GenomicInterval]` inputs avoid moving the caller's Vec; bedrs builds its `IntervalContainer` once, runs the set op, returns a fresh Vec. No intermediate clones.
5. **`Copy` newtypes (`HapId`, `QnameIdx`, `Strand`, `RgRef::Index`)** pass by value — no refs, no lifetimes threaded through call graphs.
6. **`Cow`-free, owned `PathBuf` getters.** Path getters are pure `join` (one alloc each) and are called O(stages) times, not in any loop — owning is simpler than threading a `'a` and the alloc is negligible. (Documented trade-off, not an oversight.)
7. **Drop, not leak, the collate fd (HYG-6).** Plain `drop(child_stdout)` after `from_path` (htslib dup'd it) eliminates the per-island fd leak; lets `fp-control` (T4) open hundreds of islands in one long-running process.
8. **Streaming VCF concat/sort.** `concat_sort_vcfs` does a k-way merge over already-coordinate-sorted inputs (each island VCF is sorted) → O(N) merge instead of load-all-then-sort.
9. **rayon over the independent axis only at the orchestrator.** T0 keeps libs single-threaded-by-default (`threads: u8` forwarded to htslib's own pool); the per-island `rayon` parallelism lives in T4/T9 where islands are the independent axis — avoids nested oversubscription.

---

## 6. Risks / open decisions

- **DUP-2a M-CIGAR seam (deferred to T4, flagged now).** `build_phasing_graph` tolerates `M`; `haplotype_inspection` panics on `M` (`pairwise_read_inspection.rs:83-95`). Decision recorded in CLAUDE.md/golden-encoding: **reject `M` via a typed error** (`SdError::CigarMOp`) because minimap2 `--eqx` always emits `=`/`X`. T0 only defines the error variant + the `&Record` seam; the per-read vector toolkit stays in compute crates. **Do not apply the golden encoding switch (DivB) here** — it diverges numerically from Python and lands at T4.
- **`merge_bams` header reconcile.** Python runs `modify_bam_sq_lines` (`shell_utils.sh`) to fix SQ lines against the ref before merge. The in-Rust merge must reproduce this SQ-line reconciliation or downstream tools mis-map tids. **Needs the user's call:** port the SQ-rewrite in Rust vs. keep a single leaf `samtools merge` subprocess. Recommend porting (header is small) but flag as parity risk.
- **`np.percentile(...,99)` interpolation.** NumPy uses linear interpolation between order statistics; a naive "element at 0.99·n" differs at small N. `get_insert_size_distribution` must replicate NumPy's `linear` method to match frag-size stats byte-for-byte. Also: it depends on Python `random.random()` sampling — **the differential test must tolerate Monte-Carlo variance** (assert mean within tolerance, not equality). Flag.
- **`bcftools concat -a -d exact` semantics.** `-d exact` dedups records identical on CHROM/POS/REF/ALT (and `-a` requires indexed inputs / allows overlaps). The Rust k-way merge must replicate the exact-dedup key and ordering (bcftools sort is lexicographic-by-header-contig-order, not string order). Parity risk for VCF re-emit test.
- **GraphML round-trip vs graph-tool.** graph-tool writes typed `<key>` attrs; petgraph-graphml writes string attrs. If any stage *reads back* a graph-tool-written GraphML (T8 sd-prep), the `quick-xml` reader must map graph-tool's `<data>` keys. If graphs are only ever written by Rust and read by Rust, this is moot — **confirm no cross-tool GraphML read survives the migration.**
- **`rust_htslib::errors::Error` is not cleanly `Send`/`Clone`.** `SdError::Htslib(String)` stringifies at the boundary (loses structured kind). Acceptable for a typed library error; note for callers that want to match on htslib error kinds.
- **`suppress_warning.py` has no Rust analog** (Numba/NumPy retired). Record as intentionally dropped, not migrated.
- **`Paths` side-effects split.** Python getters mkdir (`rg_dir`) and shell `samtools faidx` (`ref_genome_fai_path`). T0 keeps `Paths` pure; the mkdir/faidx move to an io-crate `ensure_layout(&Paths)` + `ensure_faidx(&Path)`. Parity: ensure call sites that relied on the getter's mkdir side-effect are updated (audit at each consuming task).
- **Dependency-availability:** no fallbacks designed. If rust-htslib/bedrs fail to build (LIBCLANG/OpenSSL env), the crate fails to compile — per the user rule, fix the env, don't add a degraded mode.

---

## 7. Test plan delta

### Unit (tier 1) — concrete fixtures
- **`Paths` string parity (golden table).** Hard-code one known run config (the HG002 example: `ref=…hg38.fasta`, `input_bam=HG002….bam`, `sd_map=…hg38…`, `output_dir=/tmp/x`, `target_bed=…CMRG…`). Assert each of the ~35 getters equals the exact string `SDrecallPaths` produces (capture the Python strings once into a `paths_golden.json` fixture). Cover `_extract_assembly_version` cascade (hg19/hg38/chm13/t2t + raise) and `normalize_rg_label` (`0`→`RG0`, `"3"`→`RG3`, `"RG7"`→`RG7`, `"x"`→Err).
- **`na_value` table** — port the 4 regex branches: `"nan"`, `"NA"`, `";nan;"`, `"na;nan"`, `".-*_ "`, `""`, `"NaN_NA"` → true; `"chr1"`, `"0"`, `"NAN1"` → false. Exact parity with `src/utils.py:44-66`.
- **`GenomicInterval`** — overlaps/contains/merge edge cases: touching `[10,20)`/`[20,30)` (no overlap, half-open), nested, identical, empty.
- **`configure_parallelism`** — `(16, 4.0)→(4,4)`, `(15,4.0)→(4,4)` (ceil), `(1,4.0)→(1,4)`.
- **DUP-1 `is_read_noisy` parity** — build synthetic `Record`s hitting each branch (unmapped, MAPQ<cut, span<75, mate-other-chrom, median-baseQ≤cut, #lowQ≥75, softclip≥75) and assert the Rust predicate == the `bam_ncls.py:115-193` truth table.
- **HYG-6 fd test** — after `build_bam_index` on a small paired BAM through the collate pipe, assert the process fd count returns to baseline (no leaked pipe fd) and `child.wait()` succeeded.
- **Round-trip** read→write→read for BED, TSV, VCF, GraphML on tiny fixtures; assert structural equality (BED/TSV byte-identical; VCF record-identical; GraphML node/edge set + attrs equal).

### Differential vs Python (tier 2) — on HG002 fixtures
- **BED re-emit:** run `sort_merge_bed` on `example/…CMRG….bed`, diff byte-for-byte against `sortBed_and_merge` output (stranded path).
- **VCF re-emit:** `concat_sort_vcfs` on the example sdrecall VCFs vs `combine_vcfs` output, compared after `bcftools norm`-free record-key normalization (CHROM,POS,REF,ALT,GT). Account for `-d exact` dedup + contig-order sort.
- **`build_bam_index` parity:** run on HG002 chr1:1633000-1635000 (the existing validated island) and assert the retained-qname set + per-chrom interval set == Python `migrate_bam_to_ncls` (already 82/82 in `bam_lappers` validation — re-assert through the consolidated reader).
- **`get_insert_size_distribution`:** assert mean/median/std within a tolerance band of the Python output (Monte-Carlo: same BAM, fixed seed if feasible; else ±2% on mean, exact on median after percentile-trim).

**Pass criterion (unchanged + sharpened):** byte-identical BED/TSV re-emit; record-identical VCF (normalized) re-emit; `Paths` getters string-identical to `SDrecallPaths`; `is_read_noisy` and `build_bam_index` set-equal to Python on the HG002 island; frag-stats within tolerance.
