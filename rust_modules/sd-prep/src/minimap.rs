//! The ONE minimap2-rs wrapper — `align_similarity`, reused by [`crate::traversal`]
//! (homologous-sequence comparison) and [`crate::intrinsic`] (counterpart→masked
//! alignment). Replaces the Python subprocess `minimap2 -x asm10/asm20 --eqx --cs
//! -c <target> <query>` (graph_traversal.py l.140) with an in-process FFI call.
//!
//! ## What it computes
//!
//! `align_similarity(query, target, preset)` builds a minimap2 index from `target`
//! (the FIRST positional of the Python CLI = the qnode region) and maps `query`
//! (the SECOND positional = the cnode region), returning `match_len / block_len`
//! of the best mapping — exactly the Python `matches / aln_len` from PAF columns 9
//! and 10 (graph_traversal.py l.152-156). Returns `0.0` when there is no mapping
//! (Python returns `(False, 0.0)`), matching the "no alignment found" branch.
//!
//! ## Version-skew note (dependency-availability rule + parity hazard)
//!
//! minimap2-rs `0.1.31` bundles **minimap2 2.30**; the Python oracle shells the
//! system `minimap2` (here `2.28-r1209`). For `asm10`/`asm20` the `matches/aln_len`
//! ratio is deterministic given identical seeds, but a minor-version change to the
//! chaining/extension heuristics can shift the ratio near the `>= 0.95` decision
//! boundary. The threshold is therefore a named const ([`SIMILARITY_THRESHOLD`]) so
//! the parity-sensitive value lives in exactly one place, and the skew is recorded
//! here. Required dep (no hand-rolled aligner) per the dependency-availability rule.

use minimap2::Aligner;
use sdrecall_utils::{Result, SdError};

/// Counterpart-acceptance similarity cutoff: a cnode is kept iff `match_len /
/// block_len >= SIMILARITY_THRESHOLD` (graph_traversal.py l.349 keeps
/// `similarity >= 0.95`). Named so the bundled-2.30-vs-system-2.28 skew has a
/// single place to be re-tuned if a differential surfaces a boundary divergence.
pub const SIMILARITY_THRESHOLD: f64 = 0.95;

/// The minimap2 preset to use — only the two the pipeline needs.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Preset {
    /// `asm10` — traversal homologous-sequence comparison (graph_traversal l.140).
    Asm10,
    /// `asm20` — intrinsic alignment of counterpart seqs vs masked genome
    /// (intrinsic_alignment.py l.61).
    Asm20,
}

/// Align `query` against `target` and return the best mapping's `match_len /
/// block_len` (the Python PAF `matches / aln_len`). Returns `0.0` if no mapping.
///
/// `target` is the index/reference (Python's first positional, the qnode region);
/// `query` is mapped against it (second positional, the cnode region). `&[u8]`
/// borrows the caller's reference slices (zero-copy). Mirrors the Python
/// single-thread CLI (`-t 1`) with `--cs -c` (CIGAR enabled).
pub fn align_similarity(query: &[u8], target: &[u8], preset: Preset) -> Result<f64> {
    if query.is_empty() || target.is_empty() {
        return Ok(0.0);
    }
    let builder = match preset {
        Preset::Asm10 => Aligner::builder().asm10(),
        Preset::Asm20 => Aligner::builder().asm20(),
    };
    let aligner = builder
        .with_cigar()
        .with_index_threads(1)
        .with_seq(target)
        .map_err(|e| SdError::Compute(format!("minimap2 index build failed: {e}")))?;

    // map(seq, cs=true, md=false, max_frag_len=None, extra_flags=None, qname).
    let mappings = aligner
        .map(query, true, false, None, None, None)
        .map_err(|e| SdError::Compute(format!("minimap2 map failed: {e}")))?;

    // Python parses the FIRST PAF line and returns its matches/aln_len. minimap2-rs
    // returns mappings already ordered with the primary/best first, so taking the
    // best mapping reproduces "the first PAF record". Pick the maximal block_len
    // mapping defensively (the chain with the most aligned bases = the representative
    // PAF line for these single-contig asm alignments).
    let best = mappings
        .iter()
        .max_by_key(|m| m.block_len)
        .filter(|m| m.block_len > 0);
    Ok(match best {
        Some(m) => m.match_len as f64 / m.block_len as f64,
        None => 0.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_full_similarity() {
        // A long enough seq for minimap2 to seed/chain; identical → similarity 1.0.
        let seq: Vec<u8> = (0..400).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let sim = align_similarity(&seq, &seq, Preset::Asm10).unwrap();
        assert!(sim > 0.99, "identical seqs similarity {sim} should be ~1.0");
    }

    #[test]
    fn empty_inputs_zero_similarity() {
        assert_eq!(align_similarity(b"", b"ACGT", Preset::Asm10).unwrap(), 0.0);
        assert_eq!(align_similarity(b"ACGT", b"", Preset::Asm20).unwrap(), 0.0);
    }

    #[test]
    fn unmappable_sequences_zero_similarity() {
        // A homopolymer query against an unrelated high-complexity target produces
        // no chainable seed → minimap2 returns no mapping → similarity 0.0 (the
        // Python "no alignment found" → (False, 0.0) branch). Using a homopolymer
        // (not a periodic ACGT) avoids spurious repetitive-seed micro-alignments.
        let target: Vec<u8> = {
            // deterministic pseudo-random high-complexity sequence (LCG).
            let mut s: u64 = 0xC0FFEE;
            (0..600)
                .map(|_| {
                    s = s
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[((s >> 40) % 4) as usize]
                })
                .collect()
        };
        let query = vec![b'A'; 500];
        let sim = align_similarity(&query, &target, Preset::Asm10).unwrap();
        assert!(
            sim < SIMILARITY_THRESHOLD,
            "unmappable similarity {sim} unexpectedly high"
        );
    }

    #[test]
    fn one_snp_high_similarity() {
        // Same seq with one base changed → similarity just below 1.0 but high.
        let seq: Vec<u8> = {
            let mut s: u64 = 0xBEEF;
            (0..600)
                .map(|_| {
                    s = s
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[((s >> 40) % 4) as usize]
                })
                .collect()
        };
        let mut q = seq.clone();
        q[300] = if q[300] == b'A' { b'T' } else { b'A' };
        let sim = align_similarity(&q, &seq, Preset::Asm10).unwrap();
        assert!(sim > 0.99 && sim <= 1.0, "single-SNP similarity {sim}");
    }
}
