//! Indel left-normalization — the ONE deliberately-external leaf.
//!
//! Python's `sort_vcf` (merge L398-429) runs, BEFORE the co-iteration:
//! ```text
//! bcftools norm -m -both -f REF --multi-overlaps 0 -a   (split multiallelics + LEFT-ALIGN)
//!   | bcftools norm -d exact                            (dedup)
//!   | bcftools filter -i 'ALT[0]!="*" && COUNT(GT="alt")>0'
//!   | bcftools sort -Oz -o OUT
//! ```
//! and the differential pass criterion is "record-identical *after* `bcftools
//! norm`". Left-alignment against a reference FASTA has subtle edge cases that
//! `bcftools` is the canonical implementation of, so per the plan's
//! external-tools policy (design §6, recommendation A) this stays as an isolated
//! leaf subprocess — a deliberate external-ALGORITHM choice, NOT a missing-dep
//! fallback. Everything else in the crate is in-process rust-htslib.

use sdrecall_utils::{Result, SdError};
use std::path::{Path, PathBuf};
use std::process::Command;

// leaf-subprocess: external algorithm, see plan policy (design §6).
/// Run the Python `sort_vcf` pipeline via `bcftools`: split-multiallelic +
/// left-align against `ref_genome`, exact-dedup, drop spanning-`*`/no-alt records,
/// coordinate-sort, write bgzipped + index. `out` is the bgzipped VCF path.
///
/// `bcftools` is a required dependency (it is on the SDrecall env PATH); a missing
/// binary or a non-zero exit is a hard error (dependency-availability rule).
pub fn sort_vcf(input: &Path, ref_genome: &Path, out: &Path, threads: u8) -> Result<()> {
    let threads = threads.to_string();
    remove_vcf_indexes(out)?;
    // One shell pipeline, mirroring the Python command exactly so the output is
    // byte-comparable. `set -o pipefail` makes any stage's failure fail the whole.
    let script = format!(
        "set -o pipefail; \
         bcftools norm --threads {t} -m -both -f {refg} --multi-overlaps 0 -a -Ou {inp} | \
         bcftools norm --threads {t} -d exact -Ou - | \
         bcftools filter --threads {t} -i 'ALT[0] != \"*\" && COUNT(GT=\"alt\") > 0' -Ou - | \
         bcftools sort -Oz -o {out} - && \
         bcftools index -f {out}",
        t = threads,
        refg = shell_quote(ref_genome),
        inp = shell_quote(input),
        out = shell_quote(out),
    );

    let status = Command::new("bash")
        .arg("-c")
        .arg(&script)
        .status()
        .map_err(|e| SdError::Vcf(format!("spawn bcftools norm: {e}")))?;

    if !status.success() {
        return Err(SdError::Vcf(format!(
            "bcftools norm/sort failed (exit {:?}) for {}",
            status.code(),
            input.display()
        )));
    }
    Ok(())
}

// leaf-subprocess: external algorithm, see plan policy.
/// The final post-merge `bcftools norm -d exact | bcftools sort` (merge L671-676):
/// exact-dedup the merged output and coordinate-sort into `out`.
pub fn norm_dedup_sort(input: &Path, out: &Path, threads: u8) -> Result<()> {
    let threads = threads.to_string();
    remove_vcf_indexes(out)?;
    let script = format!(
        "set -o pipefail; \
         bcftools norm --threads {t} -d exact -Ou {inp} | \
         bcftools sort -Oz -o {out} - && \
         bcftools index -f {out}",
        t = threads,
        inp = shell_quote(input),
        out = shell_quote(out),
    );
    let status = Command::new("bash")
        .arg("-c")
        .arg(&script)
        .status()
        .map_err(|e| SdError::Vcf(format!("spawn bcftools norm -d exact: {e}")))?;
    if !status.success() {
        return Err(SdError::Vcf(format!(
            "bcftools norm -d exact | sort failed (exit {:?})",
            status.code()
        )));
    }
    Ok(())
}

/// Minimal single-quote shell escaping for a path (paths here are pipeline-internal
/// temp files / the reference genome, never attacker-controlled, but quote anyway).
fn shell_quote(p: &Path) -> String {
    let s = p.to_string_lossy();
    format!("'{}'", s.replace('\'', "'\\''"))
}

fn remove_vcf_indexes(vcf: &Path) -> Result<()> {
    for path in vcf_index_paths(vcf) {
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => {
                return Err(SdError::Vcf(format!(
                    "remove stale VCF index {}: {e}",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

fn vcf_index_paths(vcf: &Path) -> [PathBuf; 2] {
    [
        append_path_suffix(vcf, ".csi"),
        append_path_suffix(vcf, ".tbi"),
    ]
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(suffix);
    PathBuf::from(s)
}
