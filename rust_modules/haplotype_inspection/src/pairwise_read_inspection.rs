//! Pairwise read inspection: haplotype/error vector extraction and variant counting.
//!
//! Ports Python functions from `pairwise_read_inspection.py`:
//! - `get_read_id` → `read_id`
//! - `get_hapvector_from_cigar` → `extract_hap_vector`
//! - `get_errorvector_from_cigar` → `extract_error_vector`
//! - `extract_read_qseqs` → `extract_read_qseqs`
//! - `count_snv`, `count_continuous_indel_blocks`, `count_var`, `count_continuous_blocks`
//! - `encode_base`
//!
//! Haplotype encoding:
//!    1  = reference match
//!   -4  = SNV
//!   -6  = deletion
//!   >1  = insertion marker (length × 4, placed at first pos of next ref-consuming op)
//! > -10  = padding (NaN / not-a-value)
//!
//! Error probabilities (qual_arrays) are float32 values in [0, 1] where smaller = better quality.
//!
//! # Ownership model
//!
//! - **Records are borrowed** (`&Record`) — we never take ownership of BAM records.
//! - **Caches store owned values** — `HashMap<String, Array1<…>>` needs owned keys.
//! - **Batch returns are owned** — caller takes full ownership.

use rust_htslib::bam::Record;
use rust_htslib::bam::record::Cigar;
use ndarray::Array1;
use std::collections::HashMap;
use std::sync::LazyLock;
use log::debug;


// ─── Helper functions ─────────────────────────────────────────────────────────

/// Generate a unique identifier for a BAM read, matching Python's `get_read_id`.
///
/// Format: `"qname:flag"` where flag is the SAM flag integer.
#[inline]
pub fn read_id(record: &Record) -> String {
    let qname = String::from_utf8_lossy(record.qname());
    format!("{}:{}", qname, record.flags())
}

/// Lookup table mapping every possible Phred score (0..=255) to its error
/// probability `10^(-Q/10)`. Built once on first use, then read directly —
/// avoids a `powf` call for every base in the hot error-vector path.
static PHRED_TO_PROB: LazyLock<[f32; 256]> = LazyLock::new(|| {
    let mut table = [0.0f32; 256];
    for (q, p) in table.iter_mut().enumerate() {
        *p = 10f32.powf(-(q as f32) / 10.0);
    }
    table
});

/// Convert a Phred quality score to an error probability: `10^(-Q/10)`.
/// Reads from a precomputed table; the value is identical to computing
/// `10f32.powf(-(phred as f32) / 10.0)` directly.
#[inline]
fn phred_to_prob(phred: u8) -> f32 {
    PHRED_TO_PROB[phred as usize]
}

// ─── Haplotype vector extraction ──────────────────────────────────────────────

/// Extract a haplotype vector from a BAM record's CIGAR string.
///
/// Faithfully ports Python's `get_hapvector_from_cigar` (pairwise_read_inspection.py:165-264).
/// The vector length equals the reference genome span consumed by this alignment.
///
/// Borrows the Record read-only; only needs CIGAR access. No copy of the Record is made.
///
/// # Encoding
/// -  `1`  = reference match (`=` / RefSkip)
/// - `-4`  = SNV (`X` mismatch)
/// - `-6`  = deletion (`D`)
/// - `>1`  = insertion marker (`length × 4`), placed at the first position of the
///   next reference-consuming operation after the insertion
///
/// # Panics
/// - If the CIGAR contains `M` (op 0) — requires `=`/`X` mode (--eqx).
/// - If the reference-consuming length is 0.
pub fn extract_hap_vector(record: &Record) -> Array1<i16> {
    let cigar = record.cigar();

    // Single pass: compute ref_len and assert no M ops
    let mut ref_len: usize = 0;
    for c in cigar.iter() {
        match c {
            Cigar::Match(_) => {
                panic!(
                    "CIGAR requires =/X mode but contains M (op 0). \
                     Use --eqx when aligning."
                );
            }
            Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_len += *len as usize;
            }
            _ => {}
        }
    }
    assert!(ref_len > 0, "Reference consumption length is 0 for CIGAR");
    debug!("[extract_hap_vector] ref_len={} cigar_ops={}", ref_len, cigar.len());

    let mut hapvector: Vec<i16> = Vec::with_capacity(ref_len);
    let mut query_pos: usize = 0;
    let mut pending_ins: i16 = 0; // 0 = no pending insertion

    for c in cigar.iter() {
        match c {
            // ── RefSkip (N, op 3) ──────────────────────────────────
            // Python checks this with `if operation == 3` BEFORE the match block.
            Cigar::RefSkip(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hapvector.extend(std::iter::repeat_n(1i16, n));
                } else {
                    hapvector.push(pending_ins);
                    hapvector.extend(std::iter::repeat_n(1i16, n.saturating_sub(1)));
                    pending_ins = 0;
                }
            }
            // ── Equal (=, op 7) ────────────────────────────────────
            Cigar::Equal(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hapvector.extend(std::iter::repeat_n(1i16, n));
                } else {
                    hapvector.push(pending_ins);
                    hapvector.extend(std::iter::repeat_n(1i16, n.saturating_sub(1)));
                    pending_ins = 0;
                }
                query_pos += n;
            }
            // ── Diff (X, op 8) ─────────────────────────────────────
            // Python has an N-base check here, but it is effectively dead
            // due to int8 encoding of query_sequence. We match that behavior:
            // all X bases are treated as true mismatches (→ -4).
            Cigar::Diff(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hapvector.extend(std::iter::repeat_n(-4i16, n));
                } else {
                    hapvector.push(pending_ins);
                    if n > 1 {
                        hapvector.extend(std::iter::repeat_n(-4i16, n - 1));
                    }
                    pending_ins = 0;
                }
                query_pos += n;
            }
            // ── Ins (I, op 1) ──────────────────────────────────────
            // Defer insertion marker to next ref-consuming operation.
            Cigar::Ins(len) => {
                let n = *len as usize;
                query_pos += n;
                // Python: if index > 0
                if !hapvector.is_empty() {
                    pending_ins = (n as i16) * 4;
                    debug!("[extract_hap_vector] insertion len={} marker={} at ref_pos={}", n, pending_ins, hapvector.len());
                }
            }
            // ── Del (D, op 2) ──────────────────────────────────────
            Cigar::Del(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hapvector.extend(std::iter::repeat_n(-6i16, n));
                } else {
                    hapvector.push(pending_ins);
                    if n > 1 {
                        hapvector.extend(std::iter::repeat_n(-6i16, n - 1));
                    }
                    pending_ins = 0;
                }
            }
            // ── SoftClip (S, op 4) ────────────────────────────────
            Cigar::SoftClip(len) => {
                query_pos += *len as usize;
            }
            // HardClip, Pad, Match (already checked): no-op
            _ => {}
        }
    }

    // Suppress unused warning — query_pos mirrors Python's tracking
    let _ = query_pos;

    Array1::from_vec(hapvector)
}

// ─── Error vector extraction ──────────────────────────────────────────────────

/// Extract an error probability vector from a BAM record's CIGAR and qualities.
///
/// Faithfully ports Python's `get_errorvector_from_cigar` (pairwise_read_inspection.py:269-362).
/// The vector length equals `reference_end - reference_start`.
///
/// Borrows the Record read-only; accesses CIGAR and quality arrays. No copy of the Record.
///
/// # Encoding
/// - `(0.0, 1.0]` = error probability from Phred: `10^(-Q/10)`
/// - `0.0`        = placeholder for insertions, deletions, and reference skips
///   (Python uses sentinel 99 → converted to 0 at the end)
///
/// # Key insertion behavior
/// When an insertion is encountered, the error probability at the reference position
/// **immediately before** the insertion (`ref_pos - 1`) is set to 0.0.
/// This differs from the hap vector where the insertion marker is placed at the
/// **first position of the next** ref-consuming operation.
///
/// # Panics
/// - If the CIGAR contains `M` (op 0) — requires `=`/`X` mode.
pub fn extract_error_vector(record: &Record) -> Array1<f32> {
    let cigar = record.cigar();
    let qual = record.qual();

    // Single pass: compute ref span and assert no M ops
    let mut ref_span: usize = 0;
    for c in cigar.iter() {
        match c {
            Cigar::Match(_) => {
                panic!(
                    "CIGAR requires =/X mode but contains M (op 0). \
                     Use --eqx when aligning."
                );
            }
            Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_span += *len as usize;
            }
            _ => {}
        }
    }

    debug!("[extract_error_vector] ref_span={} qual_len={}", ref_span, qual.len());

    let mut err_vector: Vec<f32> = Vec::with_capacity(ref_span);
    let mut query_pos: usize = 0;
    let mut ref_pos: usize = 0;

    for c in cigar.iter() {
        match c {
            // ── RefSkip (N, op 3): fill with 0.0 ──────────────────
            // Python: errorvector[ref_consume:...] = 99; later 99 → 0
            Cigar::RefSkip(len) => {
                let n = *len as usize;
                err_vector.extend(std::iter::repeat_n(0.0f32, n));
                ref_pos += n;
            }
            // ── SoftClip (S, op 4) ────────────────────────────────
            Cigar::SoftClip(len) => {
                query_pos += *len as usize;
            }
            // ── Equal (=, op 7): Phred → error probability ────────
            Cigar::Equal(len) => {
                let n = *len as usize;
                for i in 0..n {
                    err_vector.push(phred_to_prob(qual[query_pos + i]));
                }
                query_pos += n;
                ref_pos += n;
            }
            // ── Diff (X, op 8): same as Equal for error vector ────
            Cigar::Diff(len) => {
                let n = *len as usize;
                for i in 0..n {
                    err_vector.push(phred_to_prob(qual[query_pos + i]));
                }
                query_pos += n;
                ref_pos += n;
            }
            // ── Ins (I, op 1): mark position before insertion ─────
            // Python: errorvector[ref_consume - 1] = 99; later 99 → 0
            Cigar::Ins(len) => {
                if ref_pos > 0 {
                    debug!("[extract_error_vector] insertion overwrite at ref_pos={}", ref_pos - 1);
                    err_vector[ref_pos - 1] = 0.0;
                }
                query_pos += *len as usize;
            }
            // ── Del (D, op 2): fill with 0.0 ─────────────────────
            // Python: errorvector[ref_consume:...] = 99; later 99 → 0
            Cigar::Del(len) => {
                let n = *len as usize;
                err_vector.extend(std::iter::repeat_n(0.0f32, n));
                ref_pos += n;
            }
            _ => {}
        }
    }

    debug_assert_eq!(
        err_vector.len(),
        ref_span,
        "Error vector length {} != expected ref span {}",
        err_vector.len(),
        ref_span
    );

    Array1::from_vec(err_vector)
}

// ─── Combined extraction ──────────────────────────────────────────────────────

/// Extract both haplotype and error vectors from a BAM record.
///
/// Convenience wrapper around `extract_hap_vector` and `extract_error_vector`.
///
/// # Returns
/// `(hap_vector, err_vector, reference_start, reference_end)`
pub fn extract_hap_err_vectors(record: &Record) -> (Array1<i16>, Array1<f32>, i64, i64) {
    let hap_vector = extract_hap_vector(record);
    let err_vector = extract_error_vector(record);
    let ref_start = record.pos();
    let ref_end = record.cigar().end_pos();
    debug!(
        "[extract_hap_err_vectors] hap_len={} err_len={} ref_start={} ref_end={}",
        hap_vector.len(), err_vector.len(), ref_start, ref_end
    );
    (hap_vector, err_vector, ref_start, ref_end)
}
// ============================================================================
// Variant Counting Functions
// Ports from pairwise_read_inspection.py
// ============================================================================

/// Count SNV positions (encoded as -4) in the haplotype array.
/// Each SNV is counted separately, even if consecutive.
///
/// Ports from pairwise_read_inspection.py lines 44-51
#[inline]
pub fn count_snv(array: &Array1<i16>) -> i32 {
    array.iter().filter(|&&v| v == -4).count() as i32
}

/// Count continuous indel blocks (deletions -6 or insertions >1).
/// Consecutive indels of any type are counted as a single block.
///
/// Ports from pairwise_read_inspection.py lines 60-67
pub fn count_continuous_indel_blocks(array: &Array1<i16>) -> i32 {
    // Create boolean array: is_var = (array == -6) | (array > 1)
    let is_var = array.mapv(|v| v == -6 || v > 1);
    count_continuous_blocks(&is_var)
}

/// Count total variants (SNVs + indel blocks).
///
/// Ports from pairwise_read_inspection.py lines 70-72
#[inline]
pub fn count_var(array: &Array1<i16>) -> i32 {
    count_snv(array) + count_continuous_indel_blocks(array)
}

/// Helper: Count the number of continuous blocks of True values in a boolean array.
///
/// Ports from pairwise_read_inspection.py lines 24-41
///
/// Algorithm:
/// 1. Pad array with False at beginning and end
/// 2. Detect block starts: extended_arr[i] == True AND extended_arr[i+1] == False
/// 3. Count block starts
pub fn count_continuous_blocks(arr: &Array1<bool>) -> i32 {
    if arr.is_empty() {
        return 0;
    }

    // Pad one False to the beginning and end of arr
    let mut extended_arr = Array1::<bool>::from_elem(arr.len() + 2, false);
    extended_arr.slice_mut(ndarray::s![1..arr.len() + 1]).assign(arr);

    // Count block starts: where extended_arr[i] is True and extended_arr[i+1] is False
    let mut block_count = 0;
    for i in 0..extended_arr.len() - 1 {
        if extended_arr[i] && !extended_arr[i + 1] {
            block_count += 1;
        }
    }

    block_count
}


// ─── Read query sequence extraction ───────────────────────────────────────────

/// Encode a single ASCII base to integer: A=0, T=1, C=2, G=3, anything else=4 (N).
///
/// Matches Python's `base_dict = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}`.
#[inline]
pub fn encode_base(byte: u8) -> i8 {
    match byte {
        b'A' | b'a' => 0,
        b'T' | b't' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 3,
        _ => 4, // N or any ambiguity code
    }
}

/// Result of extracting query sequence data from a BAM record.
///
/// Ports Python's `extract_read_qseqs` (pairwise_read_inspection.py:446-462)
/// combined with `prepare_ref_query_idx_map` (pairwise_read_inspection.py:99-130).
#[derive(Clone, Debug)]
pub struct ReadQseqData {
    /// Mapping from reference offset (relative to `reference_start`) → query index.
    /// Length = `reference_end - reference_start`.
    /// Value of -1 means that reference position is a deletion (no query base).
    /// This is `ref_positions` in Python (`prepare_ref_query_idx_map` output).
    pub ref_to_query: Vec<i32>,

    /// Mapping from query index → absolute reference position.
    /// Length = query sequence length (including soft-clipped bases).
    /// Value of -1 means that query base is an insertion or soft-clip (no ref position).
    /// This is `qseq_ref_positions` in Python (`read.get_reference_positions(full_length=True)`).
    pub query_to_ref: Vec<i32>,

    /// Encoded query sequence: A=0, T=1, C=2, G=3, N=4.
    /// Length = query sequence length.
    /// This is `query_sequence_encoded` in Python.
    pub qseq_encoded: Vec<i8>,

    /// Base quality scores (raw Phred, not ASCII-offset).
    /// Length = query sequence length.
    /// This is `query_sequence_qualities` in Python.
    pub qseq_qualities: Vec<i8>,
}

/// Extract query sequence data from a BAM record by walking the CIGAR.
///
/// Builds both the query→ref and ref→query position maps in a single CIGAR walk,
/// along with the encoded query sequence and quality scores.
///
/// Ports Python's `extract_read_qseqs` + `prepare_ref_query_idx_map`.
///
/// # Arguments
/// * `record` - BAM record (borrowed read-only)
///
/// # Returns
/// `ReadQseqData` containing all four arrays.
///
/// # Panics
/// If the CIGAR contains `M` (op 0) — requires `=`/`X` mode (--eqx).
pub fn extract_read_qseqs(record: &Record) -> ReadQseqData {
    let cigar = record.cigar();
    let seq = record.seq();
    let qual = record.qual();
    let ref_start = record.pos(); // 0-based

    // Compute ref span and query length from CIGAR
    let mut ref_span: usize = 0;
    let mut query_len: usize = 0;
    for c in cigar.iter() {
        match c {
            Cigar::Match(_) => {
                panic!(
                    "CIGAR requires =/X mode but contains M (op 0). \
                     Use --eqx when aligning."
                );
            }
            Cigar::Equal(len) | Cigar::Diff(len) => {
                ref_span += *len as usize;
                query_len += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_span += *len as usize;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                query_len += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    debug!(
        "[extract_read_qseqs] ref_start={} ref_span={} query_len={} cigar_ops={}",
        ref_start, ref_span, query_len, cigar.len()
    );

    // Initialize output arrays
    let mut ref_to_query: Vec<i32> = vec![-1i32; ref_span];     // ref offset → query index
    let mut query_to_ref: Vec<i32> = vec![-1i32; query_len];    // query index → abs ref position
    let mut qseq_encoded: Vec<i8> = Vec::with_capacity(query_len);
    let mut qseq_qualities: Vec<i8> = Vec::with_capacity(query_len);

    // Encode query sequence bases
    for i in 0..seq.len() {
        qseq_encoded.push(encode_base(seq[i]));
    }

    // Copy quality scores (Phred, not ASCII-offset — pysam/htslib already converts)
    for &q in qual {
        qseq_qualities.push(q as i8);
    }

    // Walk CIGAR to build both maps simultaneously
    let mut ref_offset: usize = 0;  // offset from reference_start
    let mut query_pos: usize = 0;   // position in query sequence

    for c in cigar.iter() {
        match c {
            // ── Equal (=) or Diff (X): consumes both ref and query ────
            Cigar::Equal(len) | Cigar::Diff(len) => {
                let n = *len as usize;
                for i in 0..n {
                    ref_to_query[ref_offset + i] = (query_pos + i) as i32;
                    query_to_ref[query_pos + i] = ref_start as i32 + (ref_offset + i) as i32;
                }
                ref_offset += n;
                query_pos += n;
            }
            // ── Deletion (D): consumes ref only → ref_to_query stays -1 ───
            Cigar::Del(len) => {
                ref_offset += *len as usize;
                // ref_to_query already -1 for these positions
            }
            // ── RefSkip (N): consumes ref only → ref_to_query stays -1 ────
            Cigar::RefSkip(len) => {
                ref_offset += *len as usize;
            }
            // ── Insertion (I): consumes query only → query_to_ref stays -1 ─
            Cigar::Ins(len) => {
                query_pos += *len as usize;
                // query_to_ref already -1 for these positions
            }
            // ── SoftClip (S): consumes query only → query_to_ref stays -1 ──
            Cigar::SoftClip(len) => {
                query_pos += *len as usize;
                // query_to_ref already -1 for these positions
            }
            // HardClip, Pad: consume neither
            _ => {}
        }
    }

    debug!(
        "[extract_read_qseqs] built ref_to_query[{}] query_to_ref[{}] qseq[{}] qual[{}]",
        ref_to_query.len(), query_to_ref.len(), qseq_encoded.len(), qseq_qualities.len()
    );

    ReadQseqData {
        ref_to_query,
        query_to_ref,
        qseq_encoded,
        qseq_qualities,
    }
}

/// Extract query sequence data with caching by read ID.
///
/// Mirrors Python's caching pattern in `extract_read_qseqs` where
/// `read_ref_pos_dict[read_id]` is checked before re-computing.
///
/// # Arguments
/// * `record` - BAM record (borrowed read-only)
/// * `cache`  - Mutable reference to the cache (keyed by `"qname:flag"`)
///
/// # Returns
/// A clone of the cached `ReadQseqData`.
pub fn extract_read_qseqs_cached(
    record: &Record,
    cache: &mut HashMap<String, ReadQseqData>,
) -> ReadQseqData {
    let rid = read_id(record);
    cache
        .entry(rid)
        .or_insert_with(|| extract_read_qseqs(record))
        .clone()
}

#[cfg(test)]
mod tests {
    use super::*;
    use log::debug;
    use ndarray::array;

    // ── Tests for phred_to_prob ────────────────────────────────────────────

    #[test]
    fn test_phred_to_prob() {
        // Phred 0 → prob 1.0
        assert!((phred_to_prob(0) - 1.0).abs() < 1e-6);
        // Phred 10 → prob 0.1
        assert!((phred_to_prob(10) - 0.1).abs() < 1e-6);
        // Phred 20 → prob 0.01
        assert!((phred_to_prob(20) - 0.01).abs() < 1e-6);
        // Phred 30 → prob 0.001
        assert!((phred_to_prob(30) - 0.001).abs() < 1e-6);
    }

    // ── Tests for extract_hap_vector / extract_error_vector ────────────────
    // These tests construct BAM Records programmatically using rust-htslib.

    use rust_htslib::bam::record::CigarString;
    use rust_htslib::bam::Record as BamRecord;

    /// Helper: create a BAM record with given CIGAR, sequence, qualities, pos, and name.
    fn make_record(cigar: CigarString, seq: &[u8], qual: &[u8], pos: i64) -> BamRecord {
        make_named_record(b"testread", cigar, seq, qual, pos)
    }

    fn make_named_record(
        name: &[u8],
        cigar: CigarString,
        seq: &[u8],
        qual: &[u8],
        pos: i64,
    ) -> BamRecord {
        let mut record = BamRecord::new();
        record.set(name, Some(&cigar), seq, qual);
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        record
    }

    #[test]
    fn test_hap_vector_all_matches() {
        // CIGAR: 5= (5 sequence matches)
        let cigar = CigarString(vec![Cigar::Equal(5)]);
        let seq = b"ACGTG";
        let qual = &[30, 30, 30, 30, 30];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        assert_eq!(hap.to_vec(), vec![1i16; 5]);
    }

    #[test]
    fn test_hap_vector_with_snv() {
        // CIGAR: 3=1X2= → 6 ref bases
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]);
        let seq = b"ACGAAC";
        let qual = &[30; 6];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        assert_eq!(hap.to_vec(), vec![1, 1, 1, -4, 1, 1]);
    }

    #[test]
    fn test_hap_vector_with_deletion() {
        // CIGAR: 3=2D3= → 8 ref positions
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Del(2),
            Cigar::Equal(3),
        ]);
        let seq = b"ACGACG"; // 6 query bases (3 + 3)
        let qual = &[30; 6];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        assert_eq!(hap.to_vec(), vec![1, 1, 1, -6, -6, 1, 1, 1]);
    }

    #[test]
    fn test_hap_vector_with_insertion() {
        // CIGAR: 3=2I3= → 6 ref positions, 8 query bases
        // Insertion marker: 2 * 4 = 8, placed at first pos of next = (position 3)
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Ins(2),
            Cigar::Equal(3),
        ]);
        let seq = b"ACGTTACG"; // 3 + 2ins + 3
        let qual = &[30; 8];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        // [1, 1, 1, 8, 1, 1] — insertion marker at pos 3
        assert_eq!(hap.to_vec(), vec![1, 1, 1, 8, 1, 1]);
    }

    #[test]
    fn test_hap_vector_insertion_before_deletion() {
        // CIGAR: 3=1I2D3= → 8 ref positions
        // Insertion marker: 1 * 4 = 4, placed at first pos of Del (position 3)
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Ins(1),
            Cigar::Del(2),
            Cigar::Equal(3),
        ]);
        let seq = b"ACGTACG"; // 3 + 1ins + 3
        let qual = &[30; 7];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        // [1, 1, 1, 4, -6, 1, 1, 1] — ins marker at pos 3, then 1 del, then matches
        assert_eq!(hap.to_vec(), vec![1, 1, 1, 4, -6, 1, 1, 1]);
    }

    #[test]
    fn test_hap_vector_insertion_before_mismatch() {
        // CIGAR: 3=1I1X2= → 6 ref positions
        // Insertion marker: 1 * 4 = 4, placed at first pos of X (position 3)
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Ins(1),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]);
        let seq = b"ACGTAAC"; // 3 + 1ins + 1mm + 2
        let qual = &[30; 7];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        // Insertion marker overwrites the single mismatch position
        assert_eq!(hap.to_vec(), vec![1, 1, 1, 4, 1, 1]);
    }

    #[test]
    fn test_hap_vector_softclip_prefix() {
        // CIGAR: 2S3= → 3 ref positions, 5 query bases
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Equal(3),
        ]);
        let seq = b"TTACG";
        let qual = &[30; 5];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        assert_eq!(hap.to_vec(), vec![1, 1, 1]);
    }

    #[test]
    fn test_hap_vector_insertion_at_start_ignored() {
        // CIGAR: 2I3= → insertion at very start is ignored (Python: if index > 0)
        let cigar = CigarString(vec![
            Cigar::Ins(2),
            Cigar::Equal(3),
        ]);
        let seq = b"TTACG"; // 2ins + 3
        let qual = &[30; 5];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        // Insertion at the start is dropped; just 3 matches
        assert_eq!(hap.to_vec(), vec![1, 1, 1]);
    }

    #[test]
    fn test_error_vector_all_matches() {
        // CIGAR: 5= with quality [30, 20, 10, 40, 30]
        let cigar = CigarString(vec![Cigar::Equal(5)]);
        let seq = b"ACGTG";
        let qual = &[30u8, 20, 10, 40, 30];
        let record = make_record(cigar, seq, qual, 100);

        let err = extract_error_vector(&record);
        assert_eq!(err.len(), 5);
        assert!((err[0] - phred_to_prob(30)).abs() < 1e-7);
        assert!((err[1] - phred_to_prob(20)).abs() < 1e-7);
        assert!((err[2] - phred_to_prob(10)).abs() < 1e-7);
        assert!((err[3] - phred_to_prob(40)).abs() < 1e-7);
        assert!((err[4] - phred_to_prob(30)).abs() < 1e-7);
    }

    #[test]
    fn test_error_vector_with_deletion() {
        // CIGAR: 2=2D2= → 6 ref positions, 4 query bases
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Del(2),
            Cigar::Equal(2),
        ]);
        let seq = b"ACGT";
        let qual = &[30, 20, 25, 35];
        let record = make_record(cigar, seq, qual, 100);

        let err = extract_error_vector(&record);
        assert_eq!(err.len(), 6);
        assert!((err[0] - phred_to_prob(30)).abs() < 1e-7);
        assert!((err[1] - phred_to_prob(20)).abs() < 1e-7);
        assert_eq!(err[2], 0.0); // deletion
        assert_eq!(err[3], 0.0); // deletion
        assert!((err[4] - phred_to_prob(25)).abs() < 1e-7);
        assert!((err[5] - phred_to_prob(35)).abs() < 1e-7);
    }

    #[test]
    fn test_error_vector_with_insertion() {
        // CIGAR: 3=2I3= → 6 ref positions, 8 query bases
        // Insertion sets err_vector[ref_pos-1] = 0.0 (position 2)
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Ins(2),
            Cigar::Equal(3),
        ]);
        let seq = b"ACGTTACG"; // 3 + 2ins + 3
        let qual = &[30, 20, 25, 10, 10, 35, 40, 30]; // 8 quality values
        let record = make_record(cigar, seq, qual, 100);

        let err = extract_error_vector(&record);
        assert_eq!(err.len(), 6);
        assert!((err[0] - phred_to_prob(30)).abs() < 1e-7); // pos 0
        assert!((err[1] - phred_to_prob(20)).abs() < 1e-7); // pos 1
        assert_eq!(err[2], 0.0);                             // pos 2: overwritten by insertion
        // After insertion: query_pos = 5, so next = bases use qual[5], qual[6], qual[7]
        assert!((err[3] - phred_to_prob(35)).abs() < 1e-7); // pos 3
        assert!((err[4] - phred_to_prob(40)).abs() < 1e-7); // pos 4
        assert!((err[5] - phred_to_prob(30)).abs() < 1e-7); // pos 5
    }

    #[test]
    fn test_error_vector_with_mismatch() {
        // CIGAR: 2=1X2= → 5 ref positions
        // Mismatch gets same treatment as match in error vector
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]);
        let seq = b"ACAAC";
        let qual = &[30, 20, 15, 25, 35];
        let record = make_record(cigar, seq, qual, 100);

        let err = extract_error_vector(&record);
        assert_eq!(err.len(), 5);
        assert!((err[2] - phred_to_prob(15)).abs() < 1e-7); // mismatch base
    }

    #[test]
    fn test_hap_err_vectors_combined() {
        // CIGAR: 3=1X1I2= → 6 ref positions (3match + 1mm + 2match)
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Diff(1),
            Cigar::Ins(1),
            Cigar::Equal(2),
        ]);
        let seq = b"ACGACAC"; // 3 + 1mm + 1ins + 2
        let qual = &[30, 30, 30, 20, 10, 30, 30]; // 7 query bases
        let record = make_record(cigar, seq, qual, 200);

        let (hap, err, start, end) = extract_hap_err_vectors(&record);

        // Hap: [1, 1, 1, -4, ins_marker=4, 1]
        assert_eq!(hap.to_vec(), vec![1, 1, 1, -4, 4, 1]);

        // Err: insertion overwrites err[ref_pos-1] where ref_pos=4 after the X(1).
        assert_eq!(err.len(), 6);
        assert!((err[0] - phred_to_prob(30)).abs() < 1e-7);
        assert!((err[1] - phred_to_prob(30)).abs() < 1e-7);
        assert!((err[2] - phred_to_prob(30)).abs() < 1e-7);
        assert_eq!(err[3], 0.0); // overwritten by insertion
        assert!((err[4] - phred_to_prob(30)).abs() < 1e-7);
        assert!((err[5] - phred_to_prob(30)).abs() < 1e-7);

        assert_eq!(start, 200);
        assert_eq!(end, 206);
    }

    #[test]
    fn test_hap_vector_softclip_both_sides() {
        // CIGAR: 2S4=3S → 4 ref positions, 9 query bases
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Equal(4),
            Cigar::SoftClip(3),
        ]);
        let seq = b"TTACGTGGG";
        let qual = &[10; 9];
        let record = make_record(cigar, seq, qual, 100);

        let hap = extract_hap_vector(&record);
        assert_eq!(hap.to_vec(), vec![1, 1, 1, 1]);
    }

    // ── Tests for extract_read_qseqs ─────────────────────────────────────────

    /// Call at the start of every test to enable `RUST_LOG=debug` output.
    fn init_log() {
        let _ = env_logger::try_init();
    }

    #[test]
    fn test_encode_base() {
        init_log();
        debug!("encode_base: A={} T={} C={} G={} N={} a={} n={}",
            encode_base(b'A'), encode_base(b'T'), encode_base(b'C'),
            encode_base(b'G'), encode_base(b'N'), encode_base(b'a'),
            encode_base(b'n'));
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'T'), 1);
        assert_eq!(encode_base(b'C'), 2);
        assert_eq!(encode_base(b'G'), 3);
        assert_eq!(encode_base(b'N'), 4);
        // lowercase
        assert_eq!(encode_base(b'a'), 0);
        assert_eq!(encode_base(b't'), 1);
        assert_eq!(encode_base(b'c'), 2);
        assert_eq!(encode_base(b'g'), 3);
        // ambiguity codes → 4
        assert_eq!(encode_base(b'R'), 4);
        assert_eq!(encode_base(b'Y'), 4);
    }

    #[test]
    fn test_extract_qseqs_all_matches() {
        init_log();
        // CIGAR: 5= at pos 100 → ref span 100..105, query len 5
        let cigar = CigarString(vec![Cigar::Equal(5)]);
        let record = make_record(cigar, b"ATCGN", &[30, 25, 20, 15, 10], 100);

        let data = extract_read_qseqs(&record);

        debug!("all_matches: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}, qual={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded, data.qseq_qualities);

        // ref_to_query: offset 0→qi0, 1→qi1, 2→qi2, 3→qi3, 4→qi4
        assert_eq!(data.ref_to_query, vec![0, 1, 2, 3, 4]);
        // query_to_ref: qi0→100, qi1→101, qi2→102, qi3→103, qi4→104
        assert_eq!(data.query_to_ref, vec![100, 101, 102, 103, 104]);
        // encoded: A=0, T=1, C=2, G=3, N=4
        assert_eq!(data.qseq_encoded, vec![0, 1, 2, 3, 4]);
        // qualities
        assert_eq!(data.qseq_qualities, vec![30, 25, 20, 15, 10]);
    }

    #[test]
    fn test_extract_qseqs_with_deletion() {
        init_log();
        // CIGAR: 3= 2D 3= at pos 10 → ref span 8 (3+2+3), query len 6
        // query: A C G T A C
        // ref:   A C G - - T A C
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Del(2),
            Cigar::Equal(3),
        ]);
        let record = make_record(cigar, b"ACGTAC", &[30; 6], 10);

        let data = extract_read_qseqs(&record);

        debug!("deletion: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded);

        // ref_to_query length = 8 (ref span)
        // offsets 0,1,2 → qi 0,1,2; offsets 3,4 → -1 (deletion); offsets 5,6,7 → qi 3,4,5
        assert_eq!(data.ref_to_query, vec![0, 1, 2, -1, -1, 3, 4, 5]);

        // query_to_ref length = 6 (query len)
        // qi 0→10, 1→11, 2→12, 3→15, 4→16, 5→17
        assert_eq!(data.query_to_ref, vec![10, 11, 12, 15, 16, 17]);

        // encoded: A=0, C=2, G=3, T=1, A=0, C=2
        assert_eq!(data.qseq_encoded, vec![0, 2, 3, 1, 0, 2]);
    }

    #[test]
    fn test_extract_qseqs_with_insertion() {
        init_log();
        // CIGAR: 2= 2I 3= at pos 50 → ref span 5 (2+3), query len 7
        // query: A C T G G A T
        //         matched  ins  matched
        // ref:   A C _ _ G A T
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Ins(2),
            Cigar::Equal(3),
        ]);
        let record = make_record(cigar, b"ACTGGAT", &[30; 7], 50);

        let data = extract_read_qseqs(&record);

        debug!("insertion: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded);

        // ref_to_query length = 5 (ref span)
        // offset 0→qi0, 1→qi1, 2→qi4, 3→qi5, 4→qi6
        assert_eq!(data.ref_to_query, vec![0, 1, 4, 5, 6]);

        // query_to_ref length = 7
        // qi0→50, qi1→51, qi2→-1(ins), qi3→-1(ins), qi4→52, qi5→53, qi6→54
        assert_eq!(data.query_to_ref, vec![50, 51, -1, -1, 52, 53, 54]);
    }

    #[test]
    fn test_extract_qseqs_with_softclip() {
        init_log();
        // CIGAR: 2S 3= 1S at pos 100 → ref span 3, query len 6
        // query: [S S] A T G [S]
        let cigar = CigarString(vec![
            Cigar::SoftClip(2),
            Cigar::Equal(3),
            Cigar::SoftClip(1),
        ]);
        let record = make_record(cigar, b"CCATGA", &[10, 10, 30, 30, 30, 10], 100);

        let data = extract_read_qseqs(&record);

        debug!("softclip: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded);

        // ref_to_query length = 3
        // offset 0→qi2, 1→qi3, 2→qi4
        assert_eq!(data.ref_to_query, vec![2, 3, 4]);

        // query_to_ref length = 6
        // qi0→-1(SC), qi1→-1(SC), qi2→100, qi3→101, qi4→102, qi5→-1(SC)
        assert_eq!(data.query_to_ref, vec![-1, -1, 100, 101, 102, -1]);

        // encoded: C=2, C=2, A=0, T=1, G=3, A=0
        assert_eq!(data.qseq_encoded, vec![2, 2, 0, 1, 3, 0]);
    }

    #[test]
    fn test_extract_qseqs_with_snv() {
        init_log();
        // CIGAR: 2= 1X 2= at pos 200 → ref span 5, query len 5
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]);
        let record = make_record(cigar, b"ACTAG", &[30, 30, 25, 30, 30], 200);

        let data = extract_read_qseqs(&record);

        debug!("snv: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded);

        // Both maps are 1:1 (SNV doesn't affect positional mapping)
        assert_eq!(data.ref_to_query, vec![0, 1, 2, 3, 4]);
        assert_eq!(data.query_to_ref, vec![200, 201, 202, 203, 204]);
        // A=0, C=2, T=1, A=0, G=3
        assert_eq!(data.qseq_encoded, vec![0, 2, 1, 0, 3]);
    }

    #[test]
    fn test_extract_qseqs_complex_cigar() {
        init_log();
        // CIGAR: 1S 2= 1I 2= 1D 1= 1S at pos 0
        // query: [S] A C [I] G T _ A [S]
        //  qi:    0  1 2  3  4 5   6  7
        // ref:       A C     G T - A
        //  ro:       0 1     2 3 4 5
        let cigar = CigarString(vec![
            Cigar::SoftClip(1),
            Cigar::Equal(2),
            Cigar::Ins(1),
            Cigar::Equal(2),
            Cigar::Del(1),
            Cigar::Equal(1),
            Cigar::SoftClip(1),
        ]);
        let record = make_record(cigar, b"NACGTAAN", &[5, 30, 30, 30, 30, 30, 30, 5], 0);

        let data = extract_read_qseqs(&record);

        debug!("complex: ref_to_query={:?}, query_to_ref={:?}, qseq={:?}, qual={:?}",
            data.ref_to_query, data.query_to_ref, data.qseq_encoded, data.qseq_qualities);

        // ref span = 2(=) + 2(=) + 1(D) + 1(=) = 6
        assert_eq!(data.ref_to_query.len(), 6);
        // ro 0→qi1, ro 1→qi2, ro 2→qi4, ro 3→qi5, ro 4→-1(D), ro 5→qi6
        assert_eq!(data.ref_to_query, vec![1, 2, 4, 5, -1, 6]);

        // query len = 1(S) + 2(=) + 1(I) + 2(=) + 1(=) + 1(S) = 8
        assert_eq!(data.query_to_ref.len(), 8);
        // qi0→-1(SC), qi1→0, qi2→1, qi3→-1(I), qi4→2, qi5→3, qi6→5, qi7→-1(SC)
        assert_eq!(data.query_to_ref, vec![-1, 0, 1, -1, 2, 3, 5, -1]);

        // encoded: N=4, A=0, C=2, G=3, T=1, A=0, A=0, N=4
        assert_eq!(data.qseq_encoded, vec![4, 0, 2, 3, 1, 0, 0, 4]);
    }

    #[test]
    fn test_extract_qseqs_caching() {
        init_log();
        let cigar = CigarString(vec![Cigar::Equal(3)]);
        let record = make_record(cigar, b"ACG", &[30, 25, 20], 100);

        let mut cache: HashMap<String, ReadQseqData> = HashMap::new();

        // First call populates cache
        let data1 = extract_read_qseqs_cached(&record, &mut cache);
        assert_eq!(cache.len(), 1);

        // Second call returns cached value
        let data2 = extract_read_qseqs_cached(&record, &mut cache);
        assert_eq!(cache.len(), 1); // no new entry

        debug!("caching: data1.qseq={:?}, data2.qseq={:?}", data1.qseq_encoded, data2.qseq_encoded);
        assert_eq!(data1.qseq_encoded, data2.qseq_encoded);
        assert_eq!(data1.ref_to_query, data2.ref_to_query);
        assert_eq!(data1.query_to_ref, data2.query_to_ref);
    }

    #[test]
    fn test_extract_qseqs_matches_python_deletion_example() {
        init_log();
        // Python docstring example from prepare_ref_query_idx_map:
        // query: A,C,G,T,A,C,G,T (Between 3rd and 4th base: deletion of 2 ref positions)
        // qseq_ref_positions: [10, 11, 12, 15, 16, 17, 18, 19]
        // ref_positions output: {10:0, 11:1, 12:2, (13:-1, 14:-1), 15:3, 16:4, 17:5, 18:6, 19:7}
        //
        // CIGAR: 3= 2D 5= at pos 10
        let cigar = CigarString(vec![
            Cigar::Equal(3),
            Cigar::Del(2),
            Cigar::Equal(5),
        ]);
        let record = make_record(cigar, b"ACGTACGT", &[30; 8], 10);

        let data = extract_read_qseqs(&record);

        debug!("python_del_example: ref_to_query={:?}, query_to_ref={:?}",
            data.ref_to_query, data.query_to_ref);

        // ref span = 3 + 2 + 5 = 10 (ref offsets 0..9, abs positions 10..19)
        assert_eq!(data.ref_to_query.len(), 10);
        // Python: {10:0, 11:1, 12:2, 13:-1, 14:-1, 15:3, 16:4, 17:5, 18:6, 19:7}
        assert_eq!(data.ref_to_query, vec![0, 1, 2, -1, -1, 3, 4, 5, 6, 7]);

        // query_to_ref = [10, 11, 12, 15, 16, 17, 18, 19]
        assert_eq!(data.query_to_ref, vec![10, 11, 12, 15, 16, 17, 18, 19]);
    }

    #[test]
    fn test_extract_qseqs_matches_python_insertion_example() {
        init_log();
        // Python docstring example from prepare_ref_query_idx_map:
        // query: A,C,G,T,A,C,G,T (3rd and 4th are inserted)
        // qseq_ref_positions: [10, 11, -1, -1, 12, 13, 14, 15]
        // ref_positions output: [0, 1, 4, 5, 6, 7]
        //
        // CIGAR: 2= 2I 4= at pos 10
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Ins(2),
            Cigar::Equal(4),
        ]);
        let record = make_record(cigar, b"ACGTACGT", &[30; 8], 10);

        let data = extract_read_qseqs(&record);

        debug!("python_ins_example: ref_to_query={:?}, query_to_ref={:?}",
            data.ref_to_query, data.query_to_ref);

        // ref span = 2 + 4 = 6 (ref offsets 0..5, abs positions 10..15)
        assert_eq!(data.ref_to_query.len(), 6);
        // Python: [0, 1, 4, 5, 6, 7]
        assert_eq!(data.ref_to_query, vec![0, 1, 4, 5, 6, 7]);

        // query_to_ref = [10, 11, -1, -1, 12, 13, 14, 15]
        assert_eq!(data.query_to_ref, vec![10, 11, -1, -1, 12, 13, 14, 15]);
    }

    // ========================================================================
    // Variant Counting Tests (from pairwise_read_inspection.py)
    // ========================================================================

    #[test]
    fn test_count_snv_basic() {
        // Test basic SNV counting (SNV encoded as -4)
        let hap_vec = array![1, -4, 1, -4, -4, 1];
        assert_eq!(count_snv(&hap_vec), 3);
    }

    #[test]
    fn test_count_snv_empty() {
        let hap_vec = array![1, 1, 1];
        assert_eq!(count_snv(&hap_vec), 0);
    }

    #[test]
    fn test_count_snv_all_snv() {
        let hap_vec = array![-4, -4, -4, -4];
        assert_eq!(count_snv(&hap_vec), 4);
    }

    #[test]
    fn test_count_continuous_indel_blocks_single_deletion() {
        // Single deletion block: -6
        let hap_vec = array![1, 1, -6, -6, -6, 1, 1];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 1);
    }

    #[test]
    fn test_count_continuous_indel_blocks_single_insertion() {
        // Single insertion block: >1 (e.g., 8 means 2bp insertion)
        let hap_vec = array![1, 1, 8, 1, 1];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 1);
    }

    #[test]
    fn test_count_continuous_indel_blocks_multiple() {
        // Multiple indel blocks: deletion, insertion, deletion
        let hap_vec = array![1, -6, -6, 1, 8, 1, -6, 1];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 3);
    }

    #[test]
    fn test_count_continuous_indel_blocks_consecutive_different_types() {
        // Consecutive deletion and insertion should be ONE block
        let hap_vec = array![1, -6, 8, 1];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 1);
    }

    #[test]
    fn test_count_continuous_indel_blocks_empty() {
        let hap_vec = array![1, 1, 1, -4, -4];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 0);
    }

    #[test]
    fn test_count_continuous_indel_blocks_at_boundaries() {
        // Indel at start and end
        let hap_vec = array![-6, 1, 1, 8];
        assert_eq!(count_continuous_indel_blocks(&hap_vec), 2);
    }

    #[test]
    fn test_count_var_combined() {
        // 2 SNVs + 2 indel blocks = 4 variants
        let hap_vec = array![1, -4, -6, -6, 1, -4, 8, 1];
        assert_eq!(count_var(&hap_vec), 4);
    }

    #[test]
    fn test_count_var_only_snvs() {
        let hap_vec = array![-4, 1, -4, 1, -4];
        assert_eq!(count_var(&hap_vec), 3);
    }

    #[test]
    fn test_count_var_only_indels() {
        let hap_vec = array![1, -6, 1, 8, 1];
        assert_eq!(count_var(&hap_vec), 2);
    }

    #[test]
    fn test_count_var_empty() {
        let hap_vec = array![1, 1, 1];
        assert_eq!(count_var(&hap_vec), 0);
    }

    #[test]
    fn test_count_continuous_blocks_helper_empty() {
        let bool_arr = array![false, false, false];
        assert_eq!(count_continuous_blocks(&bool_arr), 0);
    }

    #[test]
    fn test_count_continuous_blocks_helper_single_block() {
        let bool_arr = array![false, true, true, true, false];
        assert_eq!(count_continuous_blocks(&bool_arr), 1);
    }

    #[test]
    fn test_count_continuous_blocks_helper_multiple_blocks() {
        let bool_arr = array![true, true, false, true, false, true, true, true];
        assert_eq!(count_continuous_blocks(&bool_arr), 3);
    }

    #[test]
    fn test_count_continuous_blocks_helper_all_true() {
        let bool_arr = array![true, true, true, true];
        assert_eq!(count_continuous_blocks(&bool_arr), 1);
    }

}
