//! Subprocess wrappers for external bioinformatics tools.
//!
//! These are the deliberately-external leaf subprocesses per the plan's
//! external-tools policy (design §6): alignment algorithms (minimap2) and
//! variant-calling algorithms (bcftools call) stay out-of-process. BAM
//! merge/markdup also runs via samtools subprocess for the integration pass
//! (the in-process rust-htslib merge_bams needs SQ-line reconciliation,
//! deferred).
//!
//! Each wrapper is a hard-error-on-nonzero-exit function following the
//! `vcf-ops/src/norm.rs` pattern.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use sdrecall_utils::{Result, SdError};

// ─────────────────────────── minimap2 ────────────────────────────────────

/// Align paired FASTQ files against a masked reference using Python's minimap2
/// short-read command line.
///
/// Mirrors `shell_utils.sh::independent_minimap2_masked` with the SDrecall
/// preset (`-ax sr --eqx --MD -F 1000 --end-bonus 10 -R ...`). Output is a
/// coordinate-sorted, indexed BAM.
pub fn minimap2_align(
    r1: &Path,
    r2: &Path,
    reference: &Path,
    sample_id: &str,
    output_bam: &Path,
    threads: usize,
) -> Result<()> {
    let t = threads.to_string();
    let mmi = reference.with_extension("mmi");
    let tmp_mmi = mmi.with_extension(format!("mmi.tmp.{}", std::process::id()));
    let rg = minimap2_read_group(sample_id);
    let script = format!(
        "set -o pipefail; \
         minimap2 -x sr -d {tmp_mmi} {ref_} && \
         mv -f {tmp_mmi} {mmi} && \
         minimap2 -ax sr --eqx --MD -F 1000 --end-bonus 10 -t {t} -R {rg} {mmi} {r1} {r2} \
           | samtools sort -@ {t} -o {out} - && \
         samtools index -@ {t} {out}",
        t = t,
        ref_ = sq(reference),
        mmi = sq(&mmi),
        tmp_mmi = sq(&tmp_mmi),
        rg = shquote(&rg),
        r1 = sq(r1),
        r2 = sq(r2),
        out = sq(output_bam),
    );
    run_bash(&script, "minimap2 align")
}

// ─────────────────────────── bcftools ────────────────────────────────────

/// Call variants with `bcftools mpileup | bcftools call` (multi-allelic
/// caller). Mirrors `shell_utils.sh::bcftools_call_per_RG`.
///
/// Writes a bgzipped, indexed VCF.
pub fn bcftools_call(
    input_bam: &Path,
    ref_genome: &Path,
    output_vcf: &Path,
    filter_tag: &str,
    threads: usize,
) -> Result<()> {
    let t = threads.to_string();
    remove_vcf_indexes(output_vcf)?;
    let script = format!(
        "set -o pipefail; \
         export OPENBLAS_NUM_THREADS={t}; \
         bcftools mpileup --indels-2.0 --threads {t} -A -a FORMAT/AD,FORMAT/DP -q 10 -Q 15 -f {ref_} {bam} | \
         bcftools call --threads {t} -mv -P '4e-2' -f GQ -Ou | \
         bcftools norm --threads {t} -m -both -f {ref_} --multi-overlaps 0 -a -Ou - | \
         bcftools norm --threads {t} -d exact - | \
         bcftools view --threads {t} -i 'ALT!=\"*\"' -Ov | \
         awk 'BEGIN{{FS=OFS=\"\\t\"}} {{printf \"%s\\t%s\\t%s\\t%s\\t%s\", $1, $2, $3, toupper($4), toupper($5); for(i=6;i<=NF;i++) printf \"\\t%s\", $i; printf \"\\n\"}}' | \
         bcftools filter --threads {t} -e 'GT != \"mis\"' -s {tag} - | \
         bcftools sort -Oz -o {out} && \
         bcftools index -f {out}",
        t = t,
        ref_ = sq(ref_genome),
        bam = sq(input_bam),
        out = sq(output_vcf),
        tag = shquote(filter_tag),
    );
    run_bash(&script, "bcftools call")
}

/// Subset a VCF to records overlapping a BED region file.
/// `bcftools view -R bed -Oz -o out input` + index.
pub fn bcftools_view_regions(
    input_vcf: &Path,
    region_bed: &Path,
    output_vcf: &Path,
    threads: usize,
) -> Result<()> {
    let t = threads.to_string();
    remove_vcf_indexes(output_vcf)?;
    let script = format!(
        "bcftools view --threads {t} -R {bed} -Oz -o {out} {inp} && \
         bcftools index -f {out}",
        t = t,
        bed = sq(region_bed),
        out = sq(output_vcf),
        inp = sq(input_vcf),
    );
    run_bash(&script, "bcftools view -R")
}

/// Concatenate VCF files (order-preserving) via `bcftools concat`.
pub fn bcftools_concat(inputs: &[&Path], output_vcf: &Path, threads: usize) -> Result<()> {
    if inputs.is_empty() {
        return Err(SdError::Compute("bcftools_concat: no input VCFs".into()));
    }
    let t = threads.to_string();
    remove_vcf_indexes(output_vcf)?;
    let inp_list: String = inputs.iter().map(|p| sq(p)).collect::<Vec<_>>().join(" ");
    let script = format!(
        "set -o pipefail; \
         bcftools concat --threads {t} -a -Oz -o {out} {inp} && \
         bcftools index -f {out}",
        t = t,
        out = sq(output_vcf),
        inp = inp_list,
    );
    run_bash(&script, "bcftools concat")
}

// ─────────────────────────── samtools ────────────────────────────────────

const SAMTOOLS_MERGE_CHUNK_SIZE: usize = 256;

/// Merge multiple BAMs into one coordinate-sorted, indexed BAM.
///
/// Uses `samtools merge` subprocess (the in-process `sdrecall_io::merge_bams`
/// needs SQ-line reconciliation, deferred). Inputs are passed via `samtools
/// merge -b` list files and chunked to avoid OS argv limits and file-descriptor
/// pressure. Single-input case copies directly.
pub fn samtools_merge(inputs: &[&Path], output_bam: &Path, threads: usize) -> Result<()> {
    if inputs.is_empty() {
        return Err(SdError::Compute("samtools_merge: no input BAMs".into()));
    }
    if inputs.len() == 1 {
        let script = format!(
            "cp {inp} {out} && samtools index -@ {t} {out}",
            inp = sq(inputs[0]),
            out = sq(output_bam),
            t = threads,
        );
        return run_bash(&script, "samtools cp+index (single BAM)");
    }

    let mut temp_paths = Vec::new();
    let result = samtools_merge_chunked(inputs, output_bam, threads, &mut temp_paths)
        .and_then(|_| samtools_index(output_bam, threads));

    for path in temp_paths {
        let _ = std::fs::remove_file(path);
    }

    result
}

/// The collate → fixmate → sort → markdup pipeline that deduplicates a
/// coordinate-sorted BAM. Mirrors `realign_and_recall.py:139-153`.
///
/// Writes `output_bam` (the markdup'd, coordinate-sorted, indexed BAM).
pub fn samtools_markdup_pipeline(
    input_bam: &Path,
    output_bam: &Path,
    threads: usize,
) -> Result<()> {
    let t = threads.to_string();
    let script = format!(
        "set -o pipefail; \
         samtools collate -@ {t} -Ou {inp} \
           | samtools fixmate -@ {t} -mu - - \
           | samtools sort -@ {t} -u - \
           | samtools markdup -@ {t} -r - {out} && \
         samtools index -@ {t} {out}",
        t = t,
        inp = sq(input_bam),
        out = sq(output_bam),
    );
    run_bash(&script, "samtools markdup pipeline")
}

/// Run `samtools depth -a` on a BAM, returning the output file path.
///
/// Output format: `chrom\tpos\tdepth` (1-based positions).
pub fn samtools_depth(input_bam: &Path, output_tsv: &Path, threads: usize) -> Result<()> {
    let t = threads.to_string();
    let script = format!(
        "samtools depth -@ {t} -a {bam} > {out}",
        t = t,
        bam = sq(input_bam),
        out = sq(output_tsv),
    );
    run_bash(&script, "samtools depth")
}

/// Slice a BAM to reads overlapping a BED region file.
/// `samtools view -b -L bed -o out bam` + index.
pub fn samtools_view_region(
    input_bam: &Path,
    region_bed: &Path,
    output_bam: &Path,
    threads: usize,
) -> Result<()> {
    let t = threads.to_string();
    let script = format!(
        "samtools view -@ {t} -b -L {bed} -o {out} {inp} && \
         samtools index -@ {t} {out}",
        t = t,
        bed = sq(region_bed),
        out = sq(output_bam),
        inp = sq(input_bam),
    );
    run_bash(&script, "samtools view -L")
}

/// Build a BAM index (.bai) for the given BAM.
pub fn samtools_index(bam: &Path, threads: usize) -> Result<()> {
    let script = format!("samtools index -@ {t} {bam}", t = threads, bam = sq(bam),);
    run_bash(&script, "samtools index")
}

/// Coordinate-sort a BAM into `output_bam` and build its index.
pub fn samtools_sort_index(input_bam: &Path, output_bam: &Path, threads: usize) -> Result<()> {
    let t = threads.to_string();
    let script = format!(
        "samtools sort -O bam -@ {t} -o {out} {inp} && \
         samtools index -@ {t} {out}",
        t = t,
        out = sq(output_bam),
        inp = sq(input_bam),
    );
    run_bash(&script, "samtools sort+index")
}

/// Sort a VCF with bcftools, bgzip it, and build a tabix index.
pub fn bcftools_sort_index(input_vcf: &Path, output_vcf: &Path, _threads: usize) -> Result<()> {
    remove_vcf_indexes(output_vcf)?;
    let script = format!(
        "bcftools sort -Oz -o {out} {inp} && \
         bcftools index -f {out}",
        out = sq(output_vcf),
        inp = sq(input_vcf),
    );
    run_bash(&script, "bcftools sort+tabix")
}

/// Build a fresh VCF index after deleting both possible stale sidecars.
pub fn bcftools_index_vcf(vcf: &Path, _threads: usize) -> Result<()> {
    remove_vcf_indexes(vcf)?;
    let script = format!("bcftools index -f {vcf}", vcf = sq(vcf));
    run_bash(&script, "bcftools index")
}

// ─────────────────────────── internals ───────────────────────────────────

/// Run a bash script, hard-error on non-zero exit. Mirrors
/// `vcf-ops/src/norm.rs::sort_vcf`.
pub(crate) fn run_bash(script: &str, label: &str) -> Result<()> {
    log::debug!("[tools] {label}: {script}");
    let status = Command::new("bash")
        .arg("-c")
        .arg(script)
        .status()
        .map_err(|e| SdError::Io {
            path: format!("<{label}>"),
            source: e,
        })?;
    if !status.success() {
        return Err(SdError::Compute(format!(
            "{label} failed (exit {:?})",
            status.code()
        )));
    }
    Ok(())
}

/// Minimal single-quote shell escaping for a path.
pub(crate) fn sq(p: &Path) -> String {
    let s = p.to_string_lossy();
    shquote(&s)
}

fn samtools_merge_chunked(
    inputs: &[&Path],
    output_bam: &Path,
    threads: usize,
    temp_paths: &mut Vec<PathBuf>,
) -> Result<()> {
    if inputs.len() <= SAMTOOLS_MERGE_CHUNK_SIZE {
        return samtools_merge_from_list(inputs, output_bam, threads, temp_paths);
    }

    let mut chunk_bams = Vec::new();
    for chunk in inputs.chunks(SAMTOOLS_MERGE_CHUNK_SIZE) {
        let idx = temp_paths.len();
        let chunk_bam = merge_temp_path(output_bam, idx, "bam");
        temp_paths.push(chunk_bam.clone());
        samtools_merge_from_list(chunk, &chunk_bam, threads, temp_paths)?;
        chunk_bams.push(chunk_bam);
    }

    let chunk_refs: Vec<&Path> = chunk_bams.iter().map(|p| p.as_path()).collect();
    samtools_merge_chunked(&chunk_refs, output_bam, threads, temp_paths)
}

fn samtools_merge_from_list(
    inputs: &[&Path],
    output_bam: &Path,
    threads: usize,
    temp_paths: &mut Vec<PathBuf>,
) -> Result<()> {
    let list_path = merge_temp_path(output_bam, temp_paths.len(), "list");
    temp_paths.push(list_path.clone());
    write_path_list(inputs, &list_path)?;

    let status = Command::new("samtools")
        .arg("merge")
        .arg("-@")
        .arg(threads.to_string())
        .arg("-f")
        .arg("-o")
        .arg(output_bam)
        .arg("-b")
        .arg(&list_path)
        .status()
        .map_err(|e| SdError::Io {
            path: "<samtools merge>".to_string(),
            source: e,
        })?;

    if !status.success() {
        return Err(SdError::Compute(format!(
            "samtools merge failed (exit {:?})",
            status.code()
        )));
    }

    Ok(())
}

fn write_path_list(inputs: &[&Path], list_path: &Path) -> Result<()> {
    let mut file =
        std::io::BufWriter::new(std::fs::File::create(list_path).map_err(|e| SdError::Io {
            path: list_path.display().to_string(),
            source: e,
        })?);

    for input in inputs {
        writeln!(file, "{}", input.display()).map_err(|e| SdError::Io {
            path: list_path.display().to_string(),
            source: e,
        })?;
    }
    file.flush().map_err(|e| SdError::Io {
        path: list_path.display().to_string(),
        source: e,
    })?;

    Ok(())
}

fn merge_temp_path(output_bam: &Path, idx: usize, ext: &str) -> PathBuf {
    let parent = output_bam.parent().unwrap_or_else(|| Path::new("."));
    let name = output_bam
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("merged.bam");
    parent.join(format!(
        ".{name}.merge.{}.{}.{}",
        std::process::id(),
        idx,
        ext
    ))
}

fn remove_vcf_indexes(vcf: &Path) -> Result<()> {
    for path in vcf_index_paths(vcf) {
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(source) => {
                return Err(SdError::Io {
                    path: path.display().to_string(),
                    source,
                });
            }
        }
    }
    Ok(())
}

fn vcf_index_paths(vcf: &Path) -> [PathBuf; 2] {
    [append_path_suffix(vcf, ".csi"), append_path_suffix(vcf, ".tbi")]
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(suffix);
    PathBuf::from(s)
}

fn shquote(s: &str) -> String {
    format!("'{}'", s.replace('\'', "'\\''"))
}

fn minimap2_read_group(sample_id: &str) -> String {
    format!("@RG\\tID:{sample_id}\\tLB:SureSelectXT\\tPL:ILLUMINA\\tPU:1064\\tSM:{sample_id}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shell_quoting_handles_special_chars() {
        assert_eq!(sq(Path::new("/tmp/a b.bam")), "'/tmp/a b.bam'");
        assert_eq!(sq(Path::new("/tmp/it's.bam")), "'/tmp/it'\\''s.bam'");
    }

    #[test]
    fn minimap2_read_group_uses_escaped_tabs() {
        let rg = minimap2_read_group("HG002");
        assert_eq!(
            rg,
            r"@RG\tID:HG002\tLB:SureSelectXT\tPL:ILLUMINA\tPU:1064\tSM:HG002"
        );
        assert!(!rg.contains('\t'));
    }

    #[test]
    fn merge_temp_path_stays_next_to_output() {
        let tmp = merge_temp_path(Path::new("/tmp/out.bam"), 7, "list");
        assert_eq!(tmp.parent().unwrap(), Path::new("/tmp"));
        assert!(tmp
            .file_name()
            .unwrap()
            .to_string_lossy()
            .contains(".out.bam.merge."));
        assert!(tmp.to_string_lossy().ends_with(".7.list"));
    }

    #[test]
    fn vcf_index_paths_use_htslib_sidecar_names() {
        let paths = vcf_index_paths(Path::new("/tmp/a.vcf.gz"));
        assert_eq!(paths[0], PathBuf::from("/tmp/a.vcf.gz.csi"));
        assert_eq!(paths[1], PathBuf::from("/tmp/a.vcf.gz.tbi"));
    }
}
