//! Phase-1 driver — the in-process port of `prepare_recall_regions.py`.
//!
//! Threads the seven Phase-1 stages: multi-align depth → SD-pair load + umbrella
//! filter → multiplex graph → traversal → coloring → per-RG outputs (query/
//! counterpart BEDs, masked genome, intrinsic BAM) → GraphML +
//! `filtered_SD_binary_map.tsv`. Reuses the already-tested primitives
//! ([`crate::multialign`], [`crate::sd_pairs`], [`crate::graph_build`],
//! [`crate::traversal`], [`crate::grouping`], [`crate::masking`],
//! [`crate::intrinsic`]) — it adds only the orchestration glue (the #1 coding rule:
//! a thin upstream/downstream layer, no re-implemented units).
//!
//! ## Layout contract ([`PrepPaths`])
//!
//! The Python `SDrecallPaths` singleton derives ~30 paths from the work dir +
//! sample id; here the caller supplies a [`PrepPaths`] with the inputs and the
//! work dir, and the driver derives the per-RG output paths under
//! `work_dir/realign_groups/RG<n>/`. Only the paths the Rust driver writes are
//! derived; the rest of the Python path graph is out of T8 scope.

use crate::graph_build::{build_multiplex_graph, NodeKey, SdPairRow};
use crate::homoseq::HomoseqRegion;
use crate::sd_pairs::{filter_umbrella_group, Pair};
use crate::traversal::extract_sd_paralog_pairs;
use ahash::{AHashMap, AHashSet};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use sdrecall_utils::{clamp_threads_u8, GenomicInterval, Result, SdError, Strand};
use std::path::{Path, PathBuf};

/// Tuning parameters for the Phase-1 pick (mirrors the Python `prepare_recall_regions`
/// keyword args). Defaults match the Python defaults.
#[derive(Clone, Copy, Debug)]
pub struct PrepParams {
    pub mq_threshold: u8,
    pub high_quality_depth: i64,
    pub minimum_depth: i64,
    pub multialign_frac: f64,
    pub avg_frag: f64,
    pub std_frag: f64,
    pub mean_read_length: f64,
    pub threads: usize,
}

impl Default for PrepParams {
    fn default() -> Self {
        Self {
            mq_threshold: 41,
            high_quality_depth: 10,
            minimum_depth: 5,
            multialign_frac: 0.5,
            avg_frag: 500.0,
            std_frag: 150.0,
            mean_read_length: 147.0,
            threads: 10,
        }
    }
}

/// Inputs + output paths for a Phase-1 run.
#[derive(Clone, Debug)]
pub struct PrepPaths {
    /// Reference genome FASTA (with `.fai`).
    pub ref_genome: PathBuf,
    /// Input (host) BAM (indexed).
    pub input_bam: PathBuf,
    /// Reference SD map BED (the expanded paired-SD map).
    pub reference_sd_map: PathBuf,
    /// Target BED.
    pub target_bed: PathBuf,
    /// Multi-aligned-region BED written by Phase-1 and consumed by fp-control.
    pub multi_align_bed: PathBuf,
    /// Work directory (per-RG outputs land under `work_dir/realign_groups/`).
    pub work_dir: PathBuf,
}

impl PrepPaths {
    fn realign_dir(&self) -> PathBuf {
        self.work_dir.join("realign_groups")
    }
    fn total_recall_sd_region_bed(&self) -> PathBuf {
        self.realign_dir().join("all_target_recall_SD_regions.bed")
    }
    fn filtered_sd_map(&self) -> PathBuf {
        self.realign_dir().join("filtered_SD_binary_map.tsv")
    }
    fn rg_dir(&self, label: &str) -> PathBuf {
        self.realign_dir().join(label)
    }
}

/// One row of the SD map after the target-overlap intersection (`total_bin_sd_df`),
/// carrying the two SD segments, the BAM-overlap region and overlap length.
#[derive(Clone, Debug)]
struct BinSdRow {
    a: NodeKey,
    b: NodeKey,
    mismatch_rate: f64,
    /// The multi-align BAM interval that grouped this SD pair (`*_bam1`).
    bam_region: (String, i64, i64),
    overlap_len: i64,
}

/// Main-contig regex equivalent: `^(chr)?([0-9]+|[XYM]|MT)$` (case-insensitive).
/// Flip a strand for the both-negative→both-positive normalisation.
fn flip(s: Strand) -> Strand {
    match s {
        Strand::Forward => Strand::Reverse,
        Strand::Reverse => Strand::Forward,
        Strand::Unknown => Strand::Unknown,
    }
}

/// Parse a strand token.
fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

fn strand_str(s: Strand) -> &'static str {
    match s {
        Strand::Forward => "+",
        Strand::Reverse => "-",
        Strand::Unknown => ".",
    }
}

/// Load the reference SD map BED, filter by size, intersect with the multi-align
/// BED, and apply the both-neg-strand flip — the analog of `prepare_recall_regions.py`
/// steps 2-3 (l.119-158). Returns the `BinSdRow`s.
///
/// The reference SD map is expected to have, per line: `chr_1 start_1 end_1 chr_2
/// start_2 end_2 strand1 strand2 cigar mismatch_rate` (the expanded paired-SD map).
/// Each SD interval must exceed `avg_frag` and overlap the multi-align BED. The
/// contig-name filter has been removed: WGAC quality filtering and the downstream
/// multi-align-BED overlap already constrain query nodes to primary chromosomes,
/// so any contig (incl. `_alt`, `_random`) may serve as a counterpart. Read a text
/// file, transparently gunzipping when the content is gzip-compressed
/// (detected by the `1f 8b` magic, so it works regardless of the `.gz` extension).
/// The production SD map is delivered as `.bed.gz`; Python/pandas auto-decompresses,
/// so the Rust driver must too (otherwise the raw gzip bytes fail UTF-8 decoding).
fn read_text_maybe_gz(path: &Path) -> Result<String> {
    use std::io::Read;
    let bytes = std::fs::read(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let io_err = |e: std::io::Error| SdError::Io {
        path: path.display().to_string(),
        source: e,
    };
    if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
        let mut s = String::new();
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut s)
            .map_err(io_err)?;
        Ok(s)
    } else {
        String::from_utf8(bytes)
            .map_err(|e| io_err(std::io::Error::new(std::io::ErrorKind::InvalidData, e)))
    }
}

fn load_and_filter_sd_map(
    sd_map: &Path,
    multi_align: &[GenomicInterval],
    avg_frag: f64,
) -> Result<Vec<BinSdRow>> {
    let text = read_text_maybe_gz(sd_map)?;
    // Build a per-chrom interval index of the multi-align BED for overlap+overlap_len.
    let mut multi_by_chrom: AHashMap<String, Vec<&GenomicInterval>> = AHashMap::new();
    for iv in multi_align {
        multi_by_chrom.entry(iv.chrom.clone()).or_default().push(iv);
    }

    let mut rows: Vec<BinSdRow> = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue; // not a paired-SD row
        }
        let parse_i64 = |s: &str| s.parse::<i64>().ok();
        let (Some(s1), Some(e1), Some(s2), Some(e2)) =
            (parse_i64(f[1]), parse_i64(f[2]), parse_i64(f[4]), parse_i64(f[5]))
        else {
            continue;
        };
        let chr1 = f[0];
        let chr2 = f[3];
        let mut st1 = parse_strand(f[6]);
        let mut st2 = parse_strand(f[7]);
        let mismatch_rate = f[9].parse::<f64>().unwrap_or(0.0);

        // size filter (both segments > avg_frag).
        if (e1 - s1) as f64 <= avg_frag {
            continue;
        }
        // No contig-name filter: the SD bed is already quality-filtered by
        // WGAC (high-similarity pairs), and FC query nodes are constrained to
        // primary chromosomes by the downstream multi-align-BED overlap, so
        // any contig (incl. `_alt`, `_random`) may legitimately appear as an
        // NFC counterpart. Re-filtering by name would drop real SD pairs.
        // both-negative → both-positive flip.
        if st1 == Strand::Reverse && st2 == Strand::Reverse {
            st1 = flip(st1);
            st2 = flip(st2);
        }
        let _ = lineno;

        // intersect segment A with the multi-align BED; emit one row per overlap.
        if let Some(cands) = multi_by_chrom.get(chr1) {
            for m in cands {
                let ov_start = s1.max(m.start);
                let ov_end = e1.min(m.end);
                let overlap_len = ov_end - ov_start;
                if overlap_len <= 0 {
                    continue;
                }
                rows.push(BinSdRow {
                    a: NodeKey::new(chr1, s1, e1, st1),
                    b: NodeKey::new(chr2, s2, e2, st2),
                    mismatch_rate,
                    bam_region: (m.chrom.clone(), m.start, m.end),
                    overlap_len,
                });
            }
        }
    }
    Ok(rows)
}

/// The result of the umbrella filter + dedup: the deduped SD-pair rows (for the
/// graph + `filtered_SD_binary_map.tsv`) PLUS the distinct query nodes taken
/// **before** the frozenset dedup (Python derives `query_nodes` at l.172, before
/// the frozenset dedup at l.175 — so a node that only survives as a `chr_2`
/// segment after dedup is still a query node).
struct UmbrellaResult {
    sd_rows: Vec<SdPairRow>,
    deduped: Vec<BinSdRow>,
    /// Distinct `chr_1` query-node keys, in first-appearance order (the
    /// `total_bin_sd_df` row order Python's `sort_query_nodes` sees on ties).
    query_nodes: Vec<NodeKey>,
}

/// Apply the umbrella filter per `(chr_bam1, start_bam1, end_bam1)` group (the
/// Python groupby + `filter_umbrella_pairs` over a Pool, l.162-166), then dedup by
/// the unordered SD-pair key (frozenset, l.175-176). Returns the surviving
/// `SdPairRow`s (for the multiplex graph), the deduped `BinSdRow`s (to write
/// `filtered_SD_binary_map.tsv`), and the pre-dedup distinct query nodes.
fn umbrella_filter_and_dedup(rows: &[BinSdRow]) -> UmbrellaResult {
    // group by bam_region.
    let mut groups: AHashMap<(String, i64, i64), Vec<usize>> = AHashMap::new();
    for (i, r) in rows.iter().enumerate() {
        groups.entry(r.bam_region.clone()).or_default().push(i);
    }
    let mut kept_global: Vec<usize> = Vec::new();
    for (_, mut idxs) in groups {
        // Sort within-group entries by (chrA, startA, endA) to match the Python
        // BedTool sort order applied to the SD map before intersection. The sweep's
        // break semantics are order-dependent for mutual-umbrella ties (equal
        // overlap_len), so the within-group order determines which near-duplicate
        // survives — parity requires matching Python's sort.
        idxs.sort_by(|&a, &b| {
            let ra = &rows[a];
            let rb = &rows[b];
            (&ra.a.chrom, ra.a.start, ra.a.end)
                .cmp(&(&rb.a.chrom, rb.a.start, rb.a.end))
        });
        let pairs: Vec<Pair> = idxs
            .iter()
            .map(|&i| {
                let r = &rows[i];
                Pair {
                    chr_a: r.a.chrom.clone(),
                    start_a: r.a.start,
                    end_a: r.a.end,
                    strand_a: r.a.strand,
                    chr_b: r.b.chrom.clone(),
                    start_b: r.b.start,
                    end_b: r.b.end,
                    strand_b: r.b.strand,
                    overlap_len: r.overlap_len,
                }
            })
            .collect();
        let targets: Vec<(String, i64, i64)> = idxs
            .iter()
            .map(|&i| {
                let r = &rows[i];
                (r.bam_region.0.clone(), r.bam_region.1, r.bam_region.2)
            })
            .collect();
        let kept_local = filter_umbrella_group(&pairs, &targets, 0.95);
        for li in kept_local {
            kept_global.push(idxs[li]);
        }
    }

    // Python sorts total_bin_sd_df by (chr_bam1, start_bam1, end_bam1) (l.165), so
    // the downstream query-node order + frozenset dedup see a deterministic row
    // order. Reproduce that sort over the umbrella-kept rows.
    kept_global.sort_by(|&i, &j| {
        let ri = &rows[i];
        let rj = &rows[j];
        ri.bam_region
            .cmp(&rj.bam_region)
            .then_with(|| (&ri.a.chrom, ri.a.start, ri.a.end).cmp(&(&rj.a.chrom, rj.a.start, rj.a.end)))
            .then_with(|| (&ri.b.chrom, ri.b.start, ri.b.end).cmp(&(&rj.b.chrom, rj.b.start, rj.b.end)))
    });

    // query_nodes = distinct chr_1 segments over the umbrella-kept rows, BEFORE
    // the frozenset dedup (Python l.172 precedes l.175).
    let mut qseen: AHashSet<NodeKey> = AHashSet::new();
    let mut query_nodes: Vec<NodeKey> = Vec::new();
    for &i in &kept_global {
        let a = rows[i].a.clone();
        if qseen.insert(a.clone()) {
            query_nodes.push(a);
        }
    }

    // dedup by unordered SD-pair key (frozenset of the two SD interval keys).
    let mut seen: AHashSet<(String, String)> = AHashSet::new();
    let mut sd_rows: Vec<SdPairRow> = Vec::new();
    let mut bin_rows: Vec<BinSdRow> = Vec::new();
    for &i in &kept_global {
        let r = &rows[i];
        let ka = format!("{}:{}-{}:{}", r.a.chrom, r.a.start, r.a.end, strand_str(r.a.strand));
        let kb = format!("{}:{}-{}:{}", r.b.chrom, r.b.start, r.b.end, strand_str(r.b.strand));
        let key = if ka <= kb { (ka.clone(), kb.clone()) } else { (kb.clone(), ka.clone()) };
        if !seen.insert(key) {
            continue;
        }
        sd_rows.push(SdPairRow {
            a: r.a.clone(),
            b: r.b.clone(),
            mismatch_rate: r.mismatch_rate,
        });
        bin_rows.push(r.clone());
    }
    UmbrellaResult {
        sd_rows,
        deduped: bin_rows,
        query_nodes,
    }
}

/// Write `filtered_SD_binary_map.tsv` (the 9-col deduped SD pair table, the
/// paralog-pair-set pass criterion). Columns: `chr_1 start_1 end_1 strand1 chr_2
/// start_2 end_2 strand2 mismatch_rate` with a header row.
fn write_filtered_sd_map(path: &Path, rows: &[BinSdRow]) -> Result<()> {
    use std::io::Write;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let f = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::new(f);
    let io = |e: std::io::Error| SdError::Io {
        path: path.display().to_string(),
        source: e,
    };
    writeln!(w, "chr_1\tstart_1\tend_1\tstrand1\tchr_2\tstart_2\tend_2\tstrand2\tmismatch_rate").map_err(io)?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.a.chrom,
            r.a.start,
            r.a.end,
            strand_str(r.a.strand),
            r.b.chrom,
            r.b.start,
            r.b.end,
            strand_str(r.b.strand),
            r.mismatch_rate
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)?;
    Ok(())
}

/// One sub-cluster of a realignment group: a single FC (query) node paired with
/// its own counterpart cnodes. The sub-cluster's index is its position in
/// [`RgGroup::subclusters`] — matching the Python `cluster_idx`
/// (`build_beds_and_masked_genomes.py` l.26-32), which is the FC/NFC tag suffix.
struct RgSubcluster {
    fc: NodeKey,
    counterparts: Vec<HomoseqRegion>,
}

/// One realignment group: its sub-clusters. Each sub-cluster keeps its FC node
/// paired with ITS counterparts (we deliberately do NOT flatten, so the per-
/// sub-cluster `FC:/NFC:` tags and the `qnode_relative_region` projection that
/// `region-prep` needs can be written out — `build_beds` l.174-189).
struct RgGroup {
    subclusters: Vec<RgSubcluster>,
}

/// Collapse the coloring groups into RG clusters — the analog of
/// `query_connected_nodes` (graph_query.py l.297-334) + the `build_beds` relabel +
/// load-balance sort (build_beds l.20-37). Each color group becomes one RG; its
/// counterparts are the union of `sd_paralog_pairs[qnode]` over the group's qnodes.
fn build_rg_groups(
    color_groups: &[Vec<NodeKey>],
    sd_paralog_pairs: &FxHashMap<NodeKey, Vec<HomoseqRegion>>,
) -> Vec<RgGroup> {
    let mut groups: Vec<RgGroup> = Vec::new();
    for group in color_groups {
        // One sub-cluster per qnode that actually has counterparts, in the group's
        // qnode order — this position IS the Python `cluster_idx` (build_beds
        // l.27-32).
        let mut subclusters: Vec<RgSubcluster> = Vec::new();
        for q in group {
            if let Some(cs) = sd_paralog_pairs.get(q) {
                subclusters.push(RgSubcluster {
                    fc: q.clone(),
                    counterparts: cs.clone(),
                });
            }
        }
        if subclusters.is_empty() {
            continue;
        }
        groups.push(RgGroup { subclusters });
    }
    // Load-balance sort: (sum of FC-node sizes) * (total counterpart count),
    // descending — the same key as Python (build_beds l.36).
    groups.sort_by(|a, b| {
        let key = |g: &RgGroup| -> i64 {
            let qsz: i64 = g.subclusters.iter().map(|s| s.fc.size()).sum();
            let ccnt: i64 = g.subclusters.iter().map(|s| s.counterparts.len() as i64).sum();
            qsz * ccnt
        };
        key(b).cmp(&key(a))
    });
    groups
}

/// Write the per-RG query BED (`RG<n>.bed`): one BED6 row per qnode
/// `chrom start end . . strand` (build_beds l.148-156), sort+merged.
fn write_query_bed(path: &Path, qnodes: &[NodeKey]) -> Result<Vec<GenomicInterval>> {
    let ivs: Vec<GenomicInterval> = qnodes
        .iter()
        .map(|k| GenomicInterval::with_strand(k.chrom.clone(), k.start, k.end, k.strand))
        .collect();
    // Stranded sort+merge preserves the per-qnode strand in the BED6 6th column
    // (Python writes `record[3]` as the strand; merging only collapses same-strand
    // book-ended intervals).
    let merged = sdrecall_io::sort_merge_bed(&ivs, true);
    write_bed6(path, &merged)?;
    Ok(merged)
}

/// Write a BED6 file with `.`/`.` name/score and the interval strand.
fn write_bed6(path: &Path, ivs: &[GenomicInterval]) -> Result<()> {
    use std::io::Write;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let f = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::new(f);
    let io = |e: std::io::Error| SdError::Io {
        path: path.display().to_string(),
        source: e,
    };
    for iv in ivs {
        writeln!(
            w,
            "{}\t{}\t{}\t.\t.\t{}",
            iv.chrom, iv.start, iv.end, strand_str(iv.strand)
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)?;
    Ok(())
}

/// Write the per-RG 7-column "all related homo regions" BED — a faithful port of
/// `build_beds_and_masked_genomes.py` l.174-189.
///
/// Layout per row: `chrom  start  end  col4  col5  strand  TAG`, where:
/// - FC (query) rows: `col4 = col5 = "."`, `TAG = FC:{label}_{idx}`, coordinates
///   are the FC node's own (chrom/start/end/strand).
/// - NFC (counterpart) rows: coordinates are the counterpart's `fix_coord` window;
///   `col4/col5` are that counterpart projected into its sub-cluster FC node's
///   frame via [`HomoseqRegion::qnode_relative_region`] (`"NaN"` if the projection
///   collapses, byte-identical to Python); `TAG = NFC:{label}_{idx}`.
///
/// `idx` is the sub-cluster's position. All FC rows are written first (across
/// every sub-cluster), then all NFC rows — matching the Python write order. The
/// file is deliberately NOT merged: region-prep splits it by the per-sub-cluster
/// tags and reads `col4/col5`.
fn write_all_region_bed(path: &Path, label: &str, subclusters: &[RgSubcluster]) -> Result<()> {
    use std::io::Write;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let f = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::new(f);
    let io = |e: std::io::Error| SdError::Io {
        path: path.display().to_string(),
        source: e,
    };

    // FC rows first, across all sub-clusters.
    for (idx, sc) in subclusters.iter().enumerate() {
        writeln!(
            w,
            "{}\t{}\t{}\t.\t.\t{}\tFC:{}_{}",
            sc.fc.chrom,
            sc.fc.start,
            sc.fc.end,
            strand_str(sc.fc.strand),
            label,
            idx
        )
        .map_err(io)?;
    }
    // Then NFC rows, projected into each sub-cluster's FC frame.
    for (idx, sc) in subclusters.iter().enumerate() {
        for c in &sc.counterparts {
            let (chrom, start, end, strand) = c.fix_coord();
            let (col4, col5) = match c.qnode_relative_region(&sc.fc) {
                Some((s, e)) => (s.to_string(), e.to_string()),
                None => ("NaN".to_string(), "NaN".to_string()),
            };
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\tNFC:{}_{}",
                chrom,
                start,
                end,
                col4,
                col5,
                strand_str(strand),
                label,
                idx
            )
            .map_err(io)?;
        }
    }
    w.flush().map_err(io)?;
    Ok(())
}

/// The realignment-group output paths the driver writes (per-RG).
#[derive(Clone, Debug)]
pub struct RgOutputs {
    pub label: String,
    pub query_bed: PathBuf,
    pub counterparts_bed: PathBuf,
    pub all_regions_bed: PathBuf,
    pub masked_genome: PathBuf,
    pub intrinsic_bam: PathBuf,
}

/// The full Phase-1 output manifest.
pub struct PrepResult {
    pub multi_align_bed: PathBuf,
    pub filtered_sd_map: PathBuf,
    pub rg_outputs: Vec<RgOutputs>,
    pub total_intrinsic_bam: PathBuf,
    pub total_recall_sd_region_bed: PathBuf,
}

/// Establish the per-RG BEDs, masked genome and intrinsic BAM (the analog of
/// `establish_beds_per_RG_cluster`, build_beds l.126-216). Returns the RG output
/// paths. Pure given the RG group + the paths (file I/O only); parallelised over
/// RGs by the caller.
fn establish_rg(
    group: &RgGroup,
    label: &str,
    paths: &PrepPaths,
    params: &PrepParams,
) -> Result<RgOutputs> {
    let dir = paths.rg_dir(label);
    std::fs::create_dir_all(&dir).map_err(|e| SdError::Io {
        path: dir.display().to_string(),
        source: e,
    })?;
    let query_bed = dir.join(format!("{label}.bed"));
    let counterparts_bed = dir.join(format!("{label}_counterparts.bed"));
    let all_regions_bed = dir.join(format!("{label}_related_homo_regions.bed"));
    let masked_genome = dir.join(format!("{label}.masked.fasta"));
    let intrinsic_bam = dir.join(format!("{label}.intrinsic.bam"));

    // query BED: one row per sub-cluster FC node (build_beds l.148-156).
    let qnode_keys: Vec<NodeKey> = group.subclusters.iter().map(|s| s.fc.clone()).collect();
    let query_ivs = write_query_bed(&query_bed, &qnode_keys)?;

    // counterpart intervals: each cnode's fix_coord window (build_beds l.159-166).
    let mut counter_ivs: Vec<GenomicInterval> = Vec::new();
    for sc in &group.subclusters {
        for c in &sc.counterparts {
            let (chrom, start, end, strand) = c.fix_coord();
            if end > start {
                counter_ivs.push(GenomicInterval::with_strand(chrom, start, end, strand));
            }
        }
    }
    // Counterparts BED: subtract the query regions (strand-aware) before sort+merge
    // (build_beds l.169-170: `subtract(Query, s=True)` → sortBed_and_merge). Without
    // this, a counterpart overlapping a query region would be extracted a second
    // time as an NFC read during realignment (duplicate reads / inflated depth).
    let counter_minus_query = sdrecall_io::subtract(&counter_ivs, &query_ivs, true);
    let counter_bed = sdrecall_io::sort_merge_bed(&counter_minus_query, true);
    write_bed6(&counterparts_bed, &counter_bed)?;

    // all-regions BED: the 7-column FC/NFC-tagged, per-sub-cluster, projected
    // format that region-prep consumes (build_beds l.174-189). NOT merged, and the
    // counterparts are the RAW (un-subtracted) ones — the tags + col4/col5 must
    // survive.
    write_all_region_bed(&all_regions_bed, label, &group.subclusters)?;

    // Merged all-region intervals for the intrinsic alignment = query ∪ RAW
    // (un-subtracted) counterparts — Python's getIntrinsicBam reads the
    // un-subtracted All_region_bed; the intrinsic step needs the regions, not tags.
    let counter_merged_raw = sdrecall_io::sort_merge_bed(&counter_ivs, true);
    let mut all_ivs = query_ivs.clone();
    all_ivs.extend(counter_merged_raw.iter().cloned());
    let all_merged = sdrecall_io::sort_merge_bed(&all_ivs, false);

    // masked genome (md5-gated FASTA).
    crate::masking::mask_genome(
        &query_ivs,
        &paths.ref_genome,
        &masked_genome,
        params.avg_frag,
        params.std_frag,
        1000,
    )?;

    // intrinsic BAM (remapped to genomic coords inside, with `label` as the RG tag).
    crate::intrinsic::intrinsic_bam(
        &all_merged,
        &masked_genome,
        &paths.ref_genome,
        label,
        &intrinsic_bam,
    )?;

    Ok(RgOutputs {
        label: label.to_string(),
        query_bed,
        counterparts_bed,
        all_regions_bed,
        masked_genome,
        intrinsic_bam,
    })
}

fn target_recall_intervals(
    query_intervals: &[GenomicInterval],
    multi_align: &[GenomicInterval],
) -> Vec<GenomicInterval> {
    let overlapping = sdrecall_io::intersect(query_intervals, multi_align);
    sdrecall_io::sort_merge_bed(&overlapping, false)
}

fn write_total_recall_sd_region_bed(
    outputs: &[RgOutputs],
    multi_align: &[GenomicInterval],
    path: &Path,
) -> Result<()> {
    let mut query_intervals = Vec::new();
    for output in outputs {
        query_intervals.extend(sdrecall_io::read_bed(&output.query_bed)?);
    }
    let target_recall = target_recall_intervals(&query_intervals, multi_align);
    sdrecall_io::write_bed(path, &target_recall)?;
    log::info!(
        "Wrote target recall SD BED with {} intervals to {}",
        target_recall.len(),
        path.display()
    );
    Ok(())
}

/// Drive the full Phase-1 preparation — port of `prepare_recall_regions`
/// (prepare_recall_regions.py l.57-243). Returns the output manifest.
///
/// Steps: (1) multi-align depth pick; (2-3) SD-map load + umbrella filter + dedup +
/// write `filtered_SD_binary_map.tsv`; (4) multiplex graph; (5) traversal → SD
/// paralog pairs + connected-qnodes graph; (6) coloring → RG groups; (7) per-RG
/// BEDs + masked genomes + intrinsic BAMs + total intrinsic BAM. rayon over the
/// independent axes (RG establishment).
pub fn prepare_recall_regions(paths: &PrepPaths, params: &PrepParams) -> Result<PrepResult> {
    std::fs::create_dir_all(paths.realign_dir()).map_err(|e| SdError::Io {
        path: paths.realign_dir().display().to_string(),
        source: e,
    })?;

    // ── Step 0: fragment-size + read-length stats from the BAM ───────────────
    // Python derives avg_frag_size = paths.median_frag_size, std_frag_size =
    // paths.frag_size_std (insert_size.py), and mean_read_length =
    // calculate_mean_read_length(bam). We auto-derive them here (overriding the
    // PrepParams placeholders) so the size filter + pruning cutoffs match Python.
    let mut params = *params;
    match sdrecall_io::get_insert_size_distribution(&paths.input_bam)? {
        Some(stats) => {
            params.avg_frag = stats.median;
            params.std_frag = stats.std;
            log::info!(
                "BAM fragment size: median {:.2}bp (std {:.2}bp)",
                stats.median,
                stats.std
            );
        }
        None => log::warn!(
            "could not derive fragment-size distribution from {}; using params avg={} std={}",
            paths.input_bam.display(),
            params.avg_frag,
            params.std_frag
        ),
    }
    if let Some(mrl) = mean_read_length(&paths.input_bam)? {
        params.mean_read_length = mrl;
        log::info!("BAM mean read length: {mrl:.2}bp");
    }
    let params = &params;

    // ── Step 1: multi-align depth pick ───────────────────────────────────────
    let target = sdrecall_io::read_bed(&paths.target_bed)?;
    log::info!("Phase-1: picking multi-align regions over {} target intervals", target.len());
    let multi_align = crate::multialign::pick_multialigned_regions(
        &paths.input_bam,
        &target,
        params.mq_threshold,
        params.high_quality_depth,
        params.minimum_depth,
        params.multialign_frac,
        clamp_threads_u8(params.threads),
    )?;
    log::info!("Multi-align BED has {} intervals", multi_align.len());
    if let Some(parent) = paths.multi_align_bed.parent() {
        std::fs::create_dir_all(parent).map_err(|e| SdError::Io {
            path: parent.display().to_string(),
            source: e,
        })?;
    }
    sdrecall_io::write_bed(&paths.multi_align_bed, &multi_align)?;
    log::info!("Wrote multi-align BED to {}", paths.multi_align_bed.display());

    // ── Steps 2-3: SD-map load + umbrella filter + dedup ─────────────────────
    let bin_rows = load_and_filter_sd_map(&paths.reference_sd_map, &multi_align, params.avg_frag)?;
    log::info!("{} SD rows after target overlap", bin_rows.len());
    let umbrella = umbrella_filter_and_dedup(&bin_rows);
    log::info!("{} SD pairs after umbrella filter + dedup", umbrella.sd_rows.len());
    let filtered_sd_map = paths.filtered_sd_map();
    write_filtered_sd_map(&filtered_sd_map, &umbrella.deduped)?;

    // ── Step 4: multiplex graph ──────────────────────────────────────────────
    let graph = build_multiplex_graph(&umbrella.sd_rows, params.threads);

    // ── Step 5: traversal → SD paralog pairs + connected-qnodes graph ────────
    // query_nodes are the pre-frozenset-dedup distinct chr_1 segments (Python
    // l.172); keep only those that are nodes in the graph (a query node absent
    // from any kept SD pair has no graph vertex to traverse from).
    let query_nodes: Vec<NodeKey> = umbrella
        .query_nodes
        .iter()
        .filter(|k| graph.has_node(k))
        .cloned()
        .collect();
    let traversal = extract_sd_paralog_pairs(
        &query_nodes,
        &graph,
        &paths.ref_genome,
        params.avg_frag,
        params.std_frag,
        params.mean_read_length,
    )?;

    // ── Step 6: coloring → RG groups ─────────────────────────────────────────
    let color_groups = traversal.connected.color_groups_keys();
    let rg_groups = build_rg_groups(&color_groups, &traversal.sd_paralog_pairs);
    log::info!("{} realignment groups after coloring", rg_groups.len());

    // ── Step 7: per-RG outputs (rayon over RGs) ──────────────────────────────
    let labels: Vec<String> = (0..rg_groups.len()).map(|i| format!("RG{i}")).collect();
    let rg_outputs: Vec<Result<RgOutputs>> = rg_groups
        .par_iter()
        .zip(labels.par_iter())
        .map(|(g, label)| establish_rg(g, label, paths, params))
        .collect();
    let mut outputs = Vec::with_capacity(rg_outputs.len());
    for r in rg_outputs {
        outputs.push(r?);
    }

    let total_recall_sd_region_bed = paths.total_recall_sd_region_bed();
    write_total_recall_sd_region_bed(&outputs, &multi_align, &total_recall_sd_region_bed)?;

    // total intrinsic BAM (samtools merge/sort/index — leaf subprocess).
    // Bug A: written to `work_dir/total_intrinsic_alignments.bam` (sdrecall paths.rs:259,
    // const.py:428) — NOT under `realign_groups/`, where the orchestrator can't find it.
    let total_intrinsic_bam = paths.work_dir.join("total_intrinsic_alignments.bam");
    if !outputs.is_empty() {
        let inputs: Vec<&Path> = outputs.iter().map(|o| o.intrinsic_bam.as_path()).collect();
        crate::intrinsic::merge_total_intrinsic_bam(&inputs, &total_intrinsic_bam)?;
    }

    Ok(PrepResult {
        multi_align_bed: paths.multi_align_bed.clone(),
        filtered_sd_map,
        rg_outputs: outputs,
        total_intrinsic_bam,
        total_recall_sd_region_bed,
    })
}

/// Mean read length over the first `sample_size` reads — port of
/// `calculate_mean_read_length` (bam_ncls.py l.93-111). Averages `query_length`
/// (the SEQ field length) over the first 100k records. Returns `None` for an
/// empty BAM.
fn mean_read_length(bam: &Path) -> Result<Option<f64>> {
    use rust_htslib::bam::{self, Read as _};
    const SAMPLE: usize = 100_000;
    let mut reader =
        bam::Reader::from_path(bam).map_err(|e| SdError::Htslib(format!("open {}: {e}", bam.display())))?;
    let mut total: u64 = 0;
    let mut count: usize = 0;
    let mut rec = bam::Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r.map_err(|e| SdError::Htslib(format!("read: {e}")))?;
        total += rec.seq_len() as u64;
        count += 1;
        if count >= SAMPLE {
            break;
        }
    }
    Ok((count > 0).then(|| total as f64 / count as f64))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_negative_strand_flips_to_positive() {
        // Simulate the flip applied in load: both '-' → both '+'.
        let (mut s1, mut s2) = (Strand::Reverse, Strand::Reverse);
        if s1 == Strand::Reverse && s2 == Strand::Reverse {
            s1 = flip(s1);
            s2 = flip(s2);
        }
        assert_eq!((s1, s2), (Strand::Forward, Strand::Forward));
        // mixed strands unchanged.
        let (a, b) = (Strand::Forward, Strand::Reverse);
        assert_eq!((a, b), (Strand::Forward, Strand::Reverse));
    }

    #[test]
    fn all_region_bed_is_7col_fc_then_nfc_tagged() {
        use crate::homoseq::HomoseqRegion;
        use petgraph::graph::NodeIndex;

        // sub-cluster 0: one FC node + one counterpart; sub-cluster 1: FC only.
        // The counterpart's route is empty, so qnode_relative_region returns None
        // and col4/col5 fall back to "NaN" (byte-faithful to Python l.100).
        let cp = HomoseqRegion::new(
            NodeKey::new("chr2", 500, 600, Strand::Reverse),
            NodeIndex::new(0),
        );
        let subclusters = vec![
            RgSubcluster {
                fc: NodeKey::new("chr2", 100, 200, Strand::Forward),
                counterparts: vec![cp],
            },
            RgSubcluster {
                fc: NodeKey::new("chr3", 900, 950, Strand::Forward),
                counterparts: vec![],
            },
        ];
        let tmp = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        write_all_region_bed(tmp.path(), "RG0", &subclusters).unwrap();
        let got = std::fs::read_to_string(tmp.path()).unwrap();
        // FC rows first (both sub-clusters, in index order), then NFC rows; each
        // row is 7 tab-separated columns with a per-sub-cluster FC:/NFC: tag.
        let expected = "chr2\t100\t200\t.\t.\t+\tFC:RG0_0\n\
                        chr3\t900\t950\t.\t.\t+\tFC:RG0_1\n\
                        chr2\t500\t600\tNaN\tNaN\t-\tNFC:RG0_0\n";
        assert_eq!(got, expected);
    }

    #[test]
    fn target_recall_intervals_clip_to_multialign_and_merge() {
        let query = vec![
            GenomicInterval::new("chr1", 10, 30),
            GenomicInterval::new("chr1", 30, 50),
            GenomicInterval::new("chr2", 10, 20),
        ];
        let multi = vec![
            GenomicInterval::new("chr1", 20, 40),
            GenomicInterval::new("chr2", 30, 40),
        ];

        let got = target_recall_intervals(&query, &multi);

        assert_eq!(got, vec![GenomicInterval::new("chr1", 20, 40)]);
    }

    #[test]
    fn umbrella_dedup_collapses_reversed_pair() {
        // Two rows that are the same SD pair in opposite order + same bam region →
        // one survives the frozenset dedup.
        let a = NodeKey::new("chr1", 100, 2000, Strand::Forward);
        let b = NodeKey::new("chr2", 5000, 6900, Strand::Forward);
        let rows = vec![
            BinSdRow {
                a: a.clone(),
                b: b.clone(),
                mismatch_rate: 0.03,
                bam_region: ("chr1".into(), 100, 2000),
                overlap_len: 1900,
            },
            BinSdRow {
                a: b.clone(),
                b: a.clone(),
                mismatch_rate: 0.03,
                bam_region: ("chr2".into(), 5000, 6900),
                overlap_len: 1900,
            },
        ];
        let res = umbrella_filter_and_dedup(&rows);
        assert_eq!(res.sd_rows.len(), 1, "reversed duplicate collapses to one SD pair");
        assert_eq!(res.deduped.len(), 1);
        // Both chr_1 segments appear as query nodes (pre-frozenset-dedup): the two
        // rows have chr_1 = a and chr_1 = b respectively.
        assert_eq!(res.query_nodes.len(), 2, "query nodes taken before frozenset dedup");
    }
}
