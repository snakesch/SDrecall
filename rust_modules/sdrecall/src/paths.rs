//! `Paths` — port of `src/const.py::SDrecallPaths` (the authoritative path layout).
//!
//! This is the pure-derivation half of the Python class: the `__init__`
//! derivations (`assembly`, `sample_id`, `target_tag`, `basename`, `work_dir`)
//! and **all** `*_path()` getters, producing byte-identical path strings to the
//! Python. The side-effecting parts (`os.makedirs`, `samtools faidx`,
//! `get_insert_size_distribution`, the `check_*_validity` freshness probes) are
//! intentionally **not** ported here — they are I/O and belong in `sdrecall-io`
//! per the T0 split (DESIGN §4: "Pure derivation — NO mkdir, NO samtools").
//!
//! TODO(T9-integration): move this into `sdrecall-utils::paths` once the
//! concurrent sd-prep work settles. It lives locally in the `sdrecall` crate for
//! the scaffold pass so we don't touch a crate another agent depends on.

use std::path::{Path, PathBuf};

use sdrecall_utils::{Result, SdError};

/// One input accepted everywhere an RG is named — collapses the int / numeric-str
/// / `RG\d+`-str Python overloads of `_normalize_rg_label`.
///
/// `RgRef::Index(0)` ↔ Python `0`; `RgRef::Label("RG3")` ↔ Python `"RG3"` or
/// `"3"`.
#[derive(Clone, Copy, Debug)]
pub enum RgRef<'a> {
    Index(u32),
    Label(&'a str),
}

impl From<u32> for RgRef<'_> {
    fn from(i: u32) -> Self {
        RgRef::Index(i)
    }
}

impl<'a> From<&'a str> for RgRef<'a> {
    fn from(s: &'a str) -> Self {
        RgRef::Label(s)
    }
}

/// Centralised file-path layout for one SDrecall run.
///
/// Mirrors the fields derived in `SDrecallPaths.__init__` (`src/const.py:69-127`).
/// Frag-size stats are computed in `sdrecall-io` and SET on `Paths` afterwards
/// (keeps the path layer I/O-free).
#[derive(Clone, Debug)]
pub struct Paths {
    // ── absolute input file paths (os.path.abspath in Python) ────────────
    pub ref_genome: PathBuf,
    pub input_bam: PathBuf,
    pub reference_sd_map: PathBuf,
    /// Empty string in Python when no BED is given; `None` here.
    pub target_bed: Option<PathBuf>,
    pub output_dir: PathBuf,
    pub repo_dir: PathBuf,

    // ── derived identifiers ──────────────────────────────────────────────
    pub assembly: String,
    pub sample_id: String,
    pub target_tag: String,
    pub basename: String,
    pub work_dir: PathBuf,

    // ── standard directory tree (work_dir/{recall_results,realign_groups,intermediates}) ─
    pub recall_results_dir: PathBuf,
    pub realign_groups_dir: PathBuf,
    pub tmp_dir: PathBuf,

    // ── frag-size stats: computed in sdrecall-io, set after derivation ──
    // (Python `get_insert_size_distribution` returns (avg, median, std).)
    pub avg_frag_size: Option<f64>,
    pub median_frag_size: Option<f64>,
    pub frag_size_std: Option<f64>,
}

impl Paths {
    /// Pure derivation — **NO** mkdir, **NO** samtools, **NO** insert-size.
    ///
    /// Mirrors `SDrecallPaths.__init__` (`src/const.py:69-127`) minus the
    /// side-effecting steps. The caller (orchestrator) wires the I/O (directory
    /// creation, frag-size stats) from `sdrecall-io` at T9 integration.
    ///
    /// `ref_genome_tag` short-circuits assembly detection exactly like the
    /// Python `ref_genome_tag` parameter (`src/const.py:101`).
    #[allow(clippy::too_many_arguments)]
    pub fn derive(
        ref_genome: &Path,
        input_bam: &Path,
        reference_sd_map: &Path,
        output_dir: &Path,
        target_bed: Option<&Path>,
        sample_id: Option<&str>,
        target_tag: Option<&str>,
        ref_genome_tag: Option<&str>,
        repo_dir: &Path,
    ) -> Result<Self> {
        // os.path.abspath of every input (src/const.py:93-98).
        let ref_genome = abspath(ref_genome);
        let input_bam = abspath(input_bam);
        let reference_sd_map = abspath(reference_sd_map);
        let target_bed = target_bed.map(abspath);
        let output_dir = abspath(output_dir);
        let repo_dir = abspath(repo_dir);

        // assembly: ref_genome_tag short-circuit, else cascade detection.
        let assembly = match ref_genome_tag {
            Some(tag) => tag.to_string(),
            None => extract_assembly_version(&reference_sd_map, &ref_genome)?,
        };

        // sample_id (src/const.py:105 / _extract_sample_id 166-170).
        let sample_id = match sample_id {
            Some(s) => s.to_string(),
            None => extract_sample_id(&input_bam),
        };

        // target_tag (src/const.py:108 / _extract_target_tag 172-181).
        let target_tag = match target_tag {
            Some(t) => t.to_string(),
            None => extract_target_tag(target_bed.as_deref()),
        };

        // basename = "_".join([sample_id, assembly, (target_tag?)]) + "_SDrecall"
        // (src/const.py:111-117). target_tag is always non-empty (defaults to
        // "exome"), so it is always appended — matching the Python `if target_tag`.
        let mut dir_parts = vec![sample_id.clone(), assembly.clone()];
        if !target_tag.is_empty() {
            dir_parts.push(target_tag.clone());
        }
        let basename = format!("{}_SDrecall", dir_parts.join("_"));

        let work_dir = output_dir.join(&basename);
        let recall_results_dir = work_dir.join("recall_results");
        let realign_groups_dir = work_dir.join("realign_groups");
        let tmp_dir = work_dir.join("intermediates");

        Ok(Paths {
            ref_genome,
            input_bam,
            reference_sd_map,
            target_bed,
            output_dir,
            repo_dir,
            assembly,
            sample_id,
            target_tag,
            basename,
            work_dir,
            recall_results_dir,
            realign_groups_dir,
            tmp_dir,
            avg_frag_size: None,
            median_frag_size: None,
            frag_size_std: None,
        })
    }

    // ─────────────────────────────────────────────────────────────────────
    //  recall_results_dir-anchored getters (src/const.py:230-249)
    // ─────────────────────────────────────────────────────────────────────

    /// `recall_results/{sample_id}.pooled.raw.bam` (const.py:230-232).
    pub fn pooled_raw_bam_path(&self) -> PathBuf {
        self.recall_results_dir
            .join(format!("{}.pooled.raw.bam", self.sample_id))
    }

    /// `recall_results/{sample_id}.sdrecall.raw.vcf.gz` (const.py:234-236).
    pub fn recall_raw_vcf_path(&self) -> PathBuf {
        self.recall_results_dir
            .join(format!("{}.sdrecall.raw.vcf.gz", self.sample_id))
    }

    /// `recall_results/{sample_id}.pooled.clean.bam` (const.py:238-240).
    pub fn pooled_filtered_bam_path(&self) -> PathBuf {
        self.recall_results_dir
            .join(format!("{}.pooled.clean.bam", self.sample_id))
    }

    /// `recall_results/{sample_id}.sdrecall.clean.vcf.gz` (const.py:242-244).
    pub fn recall_filtered_vcf_path(&self) -> PathBuf {
        self.recall_results_dir
            .join(format!("{}.sdrecall.clean.vcf.gz", self.sample_id))
    }

    /// `recall_results/{sample_id}.sdrecall.vcf.gz` (const.py:247-249).
    pub fn final_recall_vcf_path(&self) -> PathBuf {
        self.recall_results_dir
            .join(format!("{}.sdrecall.vcf.gz", self.sample_id))
    }

    // ─────────────────────────────────────────────────────────────────────
    //  work_dir-anchored getters (src/const.py:251-363)
    // ─────────────────────────────────────────────────────────────────────

    /// `{sample_id}_realign_meta_table.tsv` (const.py:251-253).
    pub fn realign_meta_table_path(&self) -> PathBuf {
        self.work_dir
            .join(format!("{}_realign_meta_table.tsv", self.sample_id))
    }

    /// `{basename}_qnode_grouping.graphml` (const.py:336-338).
    pub fn qnode_grouping_graph(&self) -> PathBuf {
        self.work_dir
            .join(format!("{}_qnode_grouping.graphml", self.basename))
    }

    /// `{basename}.{target_tag}.multialign.bed` (const.py:341-343).
    pub fn multi_align_bed_path(&self) -> PathBuf {
        self.work_dir.join(format!(
            "{}.{}.multialign.bed",
            self.basename, self.target_tag
        ))
    }

    /// `raw_SD_binary_map.tsv` (const.py:345-347).
    pub fn raw_sd_binary_map_path(&self) -> PathBuf {
        self.work_dir.join("raw_SD_binary_map.tsv")
    }

    /// `filtered_SD_binary_map.tsv` (const.py:349-351).
    pub fn filtered_sd_binary_map_path(&self) -> PathBuf {
        self.work_dir.join("filtered_SD_binary_map.tsv")
    }

    /// `{basename}_multiplexed_SDs.graphml` (const.py:353-355).
    pub fn multiplex_graph_path(&self) -> PathBuf {
        self.work_dir
            .join(format!("{}_multiplexed_SDs.graphml", self.basename))
    }

    /// `multiplex_graph.replace(".graphml", ".trim.annoPC.graphml")`
    /// (const.py:357-359).
    pub fn annotated_graph_path(&self) -> PathBuf {
        replace_graphml_suffix(&self.multiplex_graph_path(), ".trim.annoPC.graphml")
    }

    /// `multiplex_graph.replace(".graphml", ".directed.overlap.{chrom}.graphml")`
    /// (const.py:361-363).
    pub fn directed_graph_path(&self, chrom: &str) -> PathBuf {
        replace_graphml_suffix(
            &self.multiplex_graph_path(),
            &format!(".directed.overlap.{chrom}.graphml"),
        )
    }

    /// `{basename}target_overlapping_query_SD.bed` (const.py:365-367).
    ///
    /// NOTE: the Python has no `.` separator between `basename` and the suffix
    /// — reproduced exactly.
    pub fn target_overlapping_query_sd_bed_path(&self) -> PathBuf {
        self.work_dir
            .join(format!("{}target_overlapping_query_SD.bed", self.basename))
    }

    /// `total_intrinsic_alignments.bam` (const.py:428-430).
    pub fn total_intrinsic_bam_path(&self) -> PathBuf {
        self.work_dir.join("total_intrinsic_alignments.bam")
    }

    /// `data/{assembly}/raw_intrin_align/assembly_intrinsic_align.WGAC.{assembly}.reformat.bam`
    /// under `repo_dir` (const.py:432-440). The Python raises on missing/too-small;
    /// the existence/size checks are I/O and deferred to T9 integration — this
    /// getter returns the pure path only.
    pub fn raw_intrinsic_bam_path(&self) -> PathBuf {
        self.repo_dir
            .join("data")
            .join(&self.assembly)
            .join("raw_intrin_align")
            .join(format!(
                "assembly_intrinsic_align.WGAC.{}.reformat.bam",
                self.assembly
            ))
    }

    /// `realign_groups/all_target_recall_SD_regions.bed` (const.py:452-454).
    pub fn total_recall_sd_region_bed_path(&self) -> PathBuf {
        self.realign_groups_dir
            .join("all_target_recall_SD_regions.bed")
    }

    /// `all_RG_related_homo_regions.bed` (const.py:472-474).
    pub fn total_homo_regions_bed_path(&self) -> PathBuf {
        self.work_dir.join("all_RG_related_homo_regions.bed")
    }

    /// FAI index path `{ref_genome}.fai` (const.py:256-263, path only — the
    /// `samtools faidx` side-effect is I/O, deferred to T9 integration).
    pub fn ref_genome_fai_path(&self) -> PathBuf {
        append_suffix(&self.ref_genome, ".fai")
    }

    // ─────────────────────────────────────────────────────────────────────
    //  RG (realign-group) getters — each normalizes the label internally
    //  (src/const.py:370-446). PURE joins (Python mkdir side-effect dropped).
    // ─────────────────────────────────────────────────────────────────────

    /// `realign_groups/{RGn}` (const.py:370-375). Pure join; mkdir is an
    /// io-crate concern.
    pub fn rg_dir(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.realign_groups_dir.join(label))
    }

    /// `recall_results/{sample_id}.sdrecall.only_{RGn}.raw.bam` (const.py:377-380).
    pub fn rg_raw_masked_bam_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.recall_results_dir.join(format!(
            "{}.sdrecall.only_{}.raw.bam",
            self.sample_id, label
        )))
    }

    /// `(r1, r2)` fastq pair under `recall_results` (const.py:382-386).
    pub fn rg_realign_fastqs_path(&self, rg: RgRef<'_>) -> Result<(PathBuf, PathBuf)> {
        let label = Self::normalize_rg_label(rg)?;
        let r1 = self.recall_results_dir.join(format!(
            "{}.sdrecall.only_{}.r1.fastq",
            self.sample_id, label
        ));
        let r2 = self.recall_results_dir.join(format!(
            "{}.sdrecall.only_{}.r2.fastq",
            self.sample_id, label
        ));
        Ok((r1, r2))
    }

    /// `{RGdir}/{RGn}.bed` (const.py:388-391).
    pub fn rg_query_bed_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.rg_dir(rg)?.join(format!("{label}.bed")))
    }

    /// `{RGdir}/{RGn}_counterparts.bed` (const.py:393-396).
    pub fn rg_counterparts_bed_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.rg_dir(rg)?.join(format!("{label}_counterparts.bed")))
    }

    /// `{RGdir}/{RGn}_related_homo_regions.bed` (const.py:398-401).
    pub fn all_homo_regions_bed_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self
            .rg_dir(rg)?
            .join(format!("{label}_related_homo_regions.bed")))
    }

    /// `{RGdir}/{RGn}_related_homo_regions.raw.fastq` (const.py:403-406).
    pub fn all_homo_regions_fastq_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self
            .rg_dir(rg)?
            .join(format!("{label}_related_homo_regions.raw.fastq")))
    }

    /// `{RGdir}/{RGn}.masked.fasta` (const.py:408-411).
    pub fn masked_genome_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.rg_dir(rg)?.join(format!("{label}.masked.fasta")))
    }

    /// `{masked_genome}.fai` (const.py:413-416).
    pub fn masked_genome_fai_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        Ok(append_suffix(&self.masked_genome_path(rg)?, ".fai"))
    }

    /// `{RGdir}/{RGn}.masked.contigsize.genome` (const.py:418-421).
    pub fn masked_genome_contigsize_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self
            .rg_dir(rg)?
            .join(format!("{label}.masked.contigsize.genome")))
    }

    /// `{RGdir}/{RGn}.masked.mmi` (const.py:423-426).
    pub fn minimap_index_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.rg_dir(rg)?.join(format!("{label}.masked.mmi")))
    }

    /// `{RGdir}/{RGn}.intrinsic.bam` (const.py:442-445).
    pub fn intrinsic_bam_path(&self, rg: RgRef<'_>) -> Result<PathBuf> {
        let label = Self::normalize_rg_label(rg)?;
        Ok(self.rg_dir(rg)?.join(format!("{label}.intrinsic.bam")))
    }

    /// Convert int / numeric-str / `RG\d+`-str → `RG{n}` (const.py:301-330).
    /// The one RG-label normalizer.
    pub fn normalize_rg_label(rg: RgRef<'_>) -> Result<String> {
        match rg {
            RgRef::Index(i) => Ok(format!("RG{i}")),
            RgRef::Label(s) => {
                if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) {
                    // numeric string → treat as index (Python str.isdigit()).
                    Ok(format!("RG{s}"))
                } else if is_rg_label(s) {
                    // already matches ^RG\d+$ → use directly.
                    Ok(s.to_string())
                } else {
                    Err(SdError::InvalidRgLabel(s.to_string()))
                }
            }
        }
    }
}

// ───────────────────────── pure derivation helpers ──────────────────────

/// `os.path.abspath` semantics: join with cwd if relative, then lexically
/// normalize. Rust's `Path` has no built-in abspath; this mirrors the common
/// case (already-absolute inputs pass through; relative inputs are joined to
/// the current dir). No filesystem access (matches `os.path.abspath`, which is
/// purely lexical).
fn abspath(p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|cwd| cwd.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    }
}

/// Append a literal suffix to a path (e.g. `.fai`) without touching extension
/// parsing — mirrors Python's `f"{path}.fai"` string concatenation.
fn append_suffix(p: &Path, suffix: &str) -> PathBuf {
    let mut s = p.as_os_str().to_os_string();
    s.push(suffix);
    PathBuf::from(s)
}

/// Mirror Python `path.replace(".graphml", new_suffix)`. The multiplex graph
/// path always ends in `.graphml`, so this swaps that exact tail.
fn replace_graphml_suffix(p: &Path, new_suffix: &str) -> PathBuf {
    let s = p.to_string_lossy();
    let replaced = s.replace(".graphml", new_suffix);
    PathBuf::from(replaced)
}

/// `_extract_assembly_version` (const.py:130-164). Substring-match cascade
/// over the sd-map basename first, then ref-genome (full path), then ref
/// basename; raise on miss.
fn extract_assembly_version(reference_sd_map: &Path, ref_genome: &Path) -> Result<String> {
    let sd_map_name = basename(reference_sd_map);
    if let Some(a) = match_assembly(&sd_map_name) {
        return Ok(a);
    }

    // Python checks the full ref_genome string here (not just basename).
    let ref_genome_str = ref_genome.to_string_lossy();
    if let Some(a) = match_assembly(&ref_genome_str) {
        return Ok(a);
    }

    // Fallback: ref basename. NOTE: the Python's final cascade has a quirk —
    // its chm13/t2t branches test `sd_map_name`, not `filename`. We reproduce
    // that exactly so behaviour is byte-identical.
    let filename = basename(ref_genome);
    if filename.contains("hg19") || filename.contains("GRCh37") {
        return Ok("hg19".to_string());
    }
    if filename.contains("hg38") || filename.contains("GRCh38") {
        return Ok("hg38".to_string());
    }
    if sd_map_name.contains("chm13") || sd_map_name.contains("CHM13") {
        return Ok("chm13".to_string());
    }
    if sd_map_name.contains("t2t") || sd_map_name.contains("T2T") {
        return Ok("chm13".to_string());
    }

    Err(SdError::UnknownAssembly {
        sd_map: reference_sd_map.to_string_lossy().into_owned(),
        reference: ref_genome.to_string_lossy().into_owned(),
    })
}

/// The hg19/hg38/chm13/t2t substring cascade used twice in
/// `_extract_assembly_version` (const.py:134-141, 144-151).
//
// The chm13 and t2t arms intentionally return the same "chm13" string: the
// Python cascade keeps them as separate `elif` branches (t2t is an alias for
// the CHM13 assembly), and we mirror that branch structure 1:1 for readable
// parity with the source rather than collapsing the conditions.
#[allow(clippy::if_same_then_else)]
fn match_assembly(name: &str) -> Option<String> {
    if name.contains("hg19") || name.contains("GRCh37") {
        Some("hg19".to_string())
    } else if name.contains("hg38") || name.contains("GRCh38") {
        Some("hg38".to_string())
    } else if name.contains("chm13") || name.contains("CHM13") {
        Some("chm13".to_string())
    } else if name.contains("t2t") || name.contains("T2T") {
        Some("chm13".to_string())
    } else {
        None
    }
}

/// `_extract_sample_id` (const.py:166-170): basename split on `.`, take `[0]`.
fn extract_sample_id(input_bam: &Path) -> String {
    let name = basename(input_bam);
    name.split('.').next().unwrap_or("").to_string()
}

/// `_extract_target_tag` (const.py:172-181): "exome" default / `default_target`
/// → "exome" / basename split on `.` take `[0]`.
fn extract_target_tag(target_bed: Option<&Path>) -> String {
    match target_bed {
        None => "exome".to_string(),
        Some(bed) => {
            let bed_str = bed.to_string_lossy();
            // Python checks the substring on the *full* target_bed string.
            if bed_str.is_empty() || bed_str.contains("default_target") {
                return "exome".to_string();
            }
            let name = basename(bed);
            name.split('.').next().unwrap_or("").to_string()
        }
    }
}

/// `os.path.basename` — the final path component as a String (lossy for
/// non-UTF-8, which never occurs for these genomics paths).
fn basename(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// `re.match(r'^RG\d+$', s)` — anchored at the start (Python `re.match`), and
/// `$` forces all trailing chars to be digits.
fn is_rg_label(s: &str) -> bool {
    let Some(rest) = s.strip_prefix("RG") else {
        return false;
    };
    !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build the canonical HG002 example config used as the golden table.
    /// Mirrors the example fixtures referenced in the migration docs.
    fn hg002_paths() -> Paths {
        Paths::derive(
            Path::new("/refs/Homo_sapiens_assembly38.hg38.fasta"),
            Path::new("/data/HG002.hg38.bam"),
            Path::new("/refs/sd_map.hg38.tsv"),
            Path::new("/tmp/out"),
            Some(Path::new("/refs/CMRG.bed")),
            None, // sample_id derived from BAM
            None, // target_tag derived from BED
            None, // assembly detected from sd-map
            Path::new("/repo/SDrecall"),
        )
        .expect("derive must succeed for the golden config")
    }

    #[test]
    fn derives_identifiers_like_python() {
        let p = hg002_paths();
        // _extract_sample_id: "HG002.hg38.bam".split('.')[0]
        assert_eq!(p.sample_id, "HG002");
        // _extract_assembly_version: "sd_map.hg38.tsv" contains "hg38"
        assert_eq!(p.assembly, "hg38");
        // _extract_target_tag: "CMRG.bed".split('.')[0]
        assert_eq!(p.target_tag, "CMRG");
        // basename = "_".join([sample_id, assembly, target_tag]) + "_SDrecall"
        assert_eq!(p.basename, "HG002_hg38_CMRG_SDrecall");
        // work_dir = output_dir / basename
        assert_eq!(
            p.work_dir,
            PathBuf::from("/tmp/out/HG002_hg38_CMRG_SDrecall")
        );
    }

    #[test]
    fn ref_genome_tag_short_circuits_assembly() {
        let p = Paths::derive(
            Path::new("/refs/genome.fasta"),
            Path::new("/data/SAMPLE.bam"),
            Path::new("/refs/ambiguous_map.tsv"),
            Path::new("/tmp/out"),
            None,
            None,
            None,
            Some("chm13"),
            Path::new("/repo"),
        )
        .unwrap();
        assert_eq!(p.assembly, "chm13");
        // No target BED → target_tag "exome".
        assert_eq!(p.target_tag, "exome");
        assert_eq!(p.basename, "SAMPLE_chm13_exome_SDrecall");
    }

    /// Golden-string table: assert each getter equals the exact string the
    /// Python `SDrecallPaths` produces for the HG002 config above.
    #[test]
    fn getters_match_python_strings() {
        let p = hg002_paths();
        let wd = "/tmp/out/HG002_hg38_CMRG_SDrecall";
        let rr = format!("{wd}/recall_results");
        let rg = format!("{wd}/realign_groups");

        // recall_results-anchored
        assert_eq!(
            p.pooled_raw_bam_path().to_str().unwrap(),
            format!("{rr}/HG002.pooled.raw.bam")
        );
        assert_eq!(
            p.recall_raw_vcf_path().to_str().unwrap(),
            format!("{rr}/HG002.sdrecall.raw.vcf.gz")
        );
        assert_eq!(
            p.pooled_filtered_bam_path().to_str().unwrap(),
            format!("{rr}/HG002.pooled.clean.bam")
        );
        assert_eq!(
            p.recall_filtered_vcf_path().to_str().unwrap(),
            format!("{rr}/HG002.sdrecall.clean.vcf.gz")
        );
        assert_eq!(
            p.final_recall_vcf_path().to_str().unwrap(),
            format!("{rr}/HG002.sdrecall.vcf.gz")
        );

        // work_dir-anchored
        assert_eq!(
            p.realign_meta_table_path().to_str().unwrap(),
            format!("{wd}/HG002_realign_meta_table.tsv")
        );
        assert_eq!(
            p.qnode_grouping_graph().to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecall_qnode_grouping.graphml")
        );
        assert_eq!(
            p.multi_align_bed_path().to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecall.CMRG.multialign.bed")
        );
        assert_eq!(
            p.multiplex_graph_path().to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecall_multiplexed_SDs.graphml")
        );
        // .replace(".graphml", ".trim.annoPC.graphml")
        assert_eq!(
            p.annotated_graph_path().to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecall_multiplexed_SDs.trim.annoPC.graphml")
        );
        // .replace(".graphml", ".directed.overlap.chr1.graphml")
        assert_eq!(
            p.directed_graph_path("chr1").to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecall_multiplexed_SDs.directed.overlap.chr1.graphml")
        );
        // NOTE the missing separator before "target_overlapping" — Python quirk.
        assert_eq!(
            p.target_overlapping_query_sd_bed_path().to_str().unwrap(),
            format!("{wd}/HG002_hg38_CMRG_SDrecalltarget_overlapping_query_SD.bed")
        );
        assert_eq!(
            p.total_intrinsic_bam_path().to_str().unwrap(),
            format!("{wd}/total_intrinsic_alignments.bam")
        );
        assert_eq!(
            p.total_recall_sd_region_bed_path().to_str().unwrap(),
            format!("{rg}/all_target_recall_SD_regions.bed")
        );
        assert_eq!(
            p.total_homo_regions_bed_path().to_str().unwrap(),
            format!("{wd}/all_RG_related_homo_regions.bed")
        );
        // raw_intrinsic_bam_path is repo_dir-anchored.
        assert_eq!(
            p.raw_intrinsic_bam_path().to_str().unwrap(),
            "/repo/SDrecall/data/hg38/raw_intrin_align/assembly_intrinsic_align.WGAC.hg38.reformat.bam"
        );
        // ref_genome_fai_path appends ".fai" to the abspath'd ref.
        assert_eq!(
            p.ref_genome_fai_path().to_str().unwrap(),
            "/refs/Homo_sapiens_assembly38.hg38.fasta.fai"
        );
    }

    #[test]
    fn rg_getters_match_python_strings() {
        let p = hg002_paths();
        let wd = "/tmp/out/HG002_hg38_CMRG_SDrecall";
        let rr = format!("{wd}/recall_results");
        let rgd = format!("{wd}/realign_groups/RG3");

        // index and "RG3" label and "3" numeric-string all normalize to RG3.
        for rg in [RgRef::Index(3), RgRef::Label("RG3"), RgRef::Label("3")] {
            assert_eq!(p.rg_dir(rg).unwrap().to_str().unwrap(), rgd);
        }

        assert_eq!(
            p.rg_raw_masked_bam_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rr}/HG002.sdrecall.only_RG3.raw.bam")
        );
        let (r1, r2) = p.rg_realign_fastqs_path(RgRef::Index(3)).unwrap();
        assert_eq!(
            r1.to_str().unwrap(),
            format!("{rr}/HG002.sdrecall.only_RG3.r1.fastq")
        );
        assert_eq!(
            r2.to_str().unwrap(),
            format!("{rr}/HG002.sdrecall.only_RG3.r2.fastq")
        );
        assert_eq!(
            p.rg_query_bed_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.bed")
        );
        assert_eq!(
            p.rg_counterparts_bed_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3_counterparts.bed")
        );
        assert_eq!(
            p.all_homo_regions_bed_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3_related_homo_regions.bed")
        );
        assert_eq!(
            p.all_homo_regions_fastq_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3_related_homo_regions.raw.fastq")
        );
        assert_eq!(
            p.masked_genome_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.masked.fasta")
        );
        assert_eq!(
            p.masked_genome_fai_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.masked.fasta.fai")
        );
        assert_eq!(
            p.masked_genome_contigsize_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.masked.contigsize.genome")
        );
        assert_eq!(
            p.minimap_index_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.masked.mmi")
        );
        assert_eq!(
            p.intrinsic_bam_path(RgRef::Index(3))
                .unwrap()
                .to_str()
                .unwrap(),
            format!("{rgd}/RG3.intrinsic.bam")
        );
    }

    #[test]
    fn normalize_rg_label_matches_python() {
        // int / numeric-str / RG-str → RG{n}
        assert_eq!(Paths::normalize_rg_label(RgRef::Index(0)).unwrap(), "RG0");
        assert_eq!(Paths::normalize_rg_label(RgRef::Label("3")).unwrap(), "RG3");
        assert_eq!(
            Paths::normalize_rg_label(RgRef::Label("RG7")).unwrap(),
            "RG7"
        );
        // invalid → Err (Python raises ValueError).
        assert!(Paths::normalize_rg_label(RgRef::Label("x")).is_err());
        assert!(Paths::normalize_rg_label(RgRef::Label("RGx")).is_err());
        assert!(Paths::normalize_rg_label(RgRef::Label("RG")).is_err());
        assert!(Paths::normalize_rg_label(RgRef::Label("")).is_err());
    }

    #[test]
    fn assembly_detection_cascade() {
        // hg19 from sd-map
        assert_eq!(
            extract_assembly_version(Path::new("/x/GRCh37_map.tsv"), Path::new("/y/g.fasta"))
                .unwrap(),
            "hg19"
        );
        // t2t alias → chm13, detected from ref path when sd-map is silent
        assert_eq!(
            extract_assembly_version(Path::new("/x/map.tsv"), Path::new("/y/t2t_genome.fasta"))
                .unwrap(),
            "chm13"
        );
        // unknown → Err
        assert!(
            extract_assembly_version(Path::new("/x/map.tsv"), Path::new("/y/genome.fasta"))
                .is_err()
        );
    }
}
