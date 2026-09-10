//! CIGAR-walk counts oracle (T1 tier-3, T4 golden-switch safety net).
//!
//! A Python-INDEPENDENT ground truth for the encoding-dependent variant counters. The
//! oracle reads the raw CIGAR, never the hap-vector encoding, so the *same* harness
//! validates the current (overwrite) encoding and the golden (summation) encoding —
//! only the asserted `count_snv` relationship changes. This crate now uses GOLDEN.
//!
//! Truths (encoding-agnostic):
//!   * `true #SNV          = Σ len(X ops)`
//!   * `true #indel-blocks = number of maximal runs of consecutive I/D ops`
//!
//! Findings this harness pins down for the GOLDEN (summation) encoding:
//!   * `count_continuous_indel_blocks(hap)` == ref-position #indel-blocks — ALWAYS (even compounds).
//!   * `count_snv(hap)` == true #SNV                                       — a compound mismatch
//!     under an insertion (`I` then `X`, whose value ends in digit 6) is now RECOVERED: summation
//!     keeps the `-4` signal (`-4 + 10*L`) instead of overwriting it. This is the DivB fix the old
//!     overwrite encoding lacked (where the compound `-4` was lost, undercounting by one per compound).

use haplotype_inspection::pairwise_read_inspection::{
    count_continuous_indel_blocks, count_snv, count_var, extract_hap_vector,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::Record;

/// Build a BAM record from CIGAR ops. The hap vector is derived purely from the CIGAR, so the
/// sequence content is irrelevant; we only size it to the query-consuming op length.
fn make_record(ops: &[Cigar]) -> Record {
    let query_len: usize = ops
        .iter()
        .map(|c| match c {
            Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => *n as usize,
            _ => 0,
        })
        .sum();
    let seq = vec![b'A'; query_len];
    let qual = vec![30u8; query_len];
    let cigar = CigarString(ops.to_vec());
    let mut rec = Record::new();
    rec.set(b"oracle", Some(&cigar), &seq, &qual);
    rec.set_pos(1000);
    rec.set_tid(0);
    rec.set_mapq(60);
    rec
}

// ── Encoding-agnostic CIGAR truths ────────────────────────────────────────────

fn true_snv(ops: &[Cigar]) -> i32 {
    ops.iter()
        .map(|c| if let Cigar::Diff(n) = c { *n as i32 } else { 0 })
        .sum()
}

/// Number of maximal runs of consecutive I/D ops in the op sequence.
fn true_indel_blocks(ops: &[Cigar]) -> i32 {
    let mut blocks = 0;
    let mut prev_indel = false;
    for c in ops {
        let is_indel = matches!(c, Cigar::Ins(_) | Cigar::Del(_));
        if is_indel && !prev_indel {
            blocks += 1;
        }
        prev_indel = is_indel;
    }
    blocks
}

/// Compound bases: an `I` op immediately followed by an `X` op, so the first mismatch position is
/// summed with the insertion marker (`-4 + 10*L`, ending in digit 6). The OLD overwrite encoding
/// lost the `-4` here (undercount); GOLDEN recovers it, so this is now an exact "compound count".
fn compound_count(ops: &[Cigar]) -> i32 {
    ops.windows(2)
        .filter(|w| matches!(w[0], Cigar::Ins(_)) && matches!(w[1], Cigar::Diff(_)))
        .count() as i32
}

/// Ref-position indel-block truth: mark each reference position as indel (a deletion, or the
/// position an insertion marker attaches to) directly from CIGAR structure — NOT from the i16
/// values — then count maximal runs. This is the space `count_continuous_indel_blocks` actually
/// operates in, so it stays exact even when an insertion marker consumes a lone match position
/// between two deletions and merges them (a SECOND overwrite-encoding loss mode, distinct from the
/// compound-SNV undercount, and NOT fixed by golden — insertions still attach to a ref position).
/// Validating against this catches encoding bugs where a deletion/insertion fails to register as a
/// variant position (e.g. a golden change to the deletion sentinel that forgets the count predicate).
fn ref_position_indel_blocks(ops: &[Cigar]) -> i32 {
    let mut is_indel: Vec<bool> = Vec::new();
    let mut pending_ins = false;
    for op in ops {
        match op {
            Cigar::Equal(n) | Cigar::Diff(n) | Cigar::RefSkip(n) => {
                for i in 0..*n {
                    is_indel.push(i == 0 && pending_ins);
                    if i == 0 {
                        pending_ins = false;
                    }
                }
            }
            Cigar::Del(n) => {
                is_indel.resize(is_indel.len() + *n as usize, true);
                pending_ins = false;
            }
            Cigar::Ins(_) => {
                if !is_indel.is_empty() {
                    pending_ins = true; // attaches to the next ref position
                }
            }
            _ => {}
        }
    }
    let mut blocks = 0;
    let mut prev = false;
    for &b in &is_indel {
        if b && !prev {
            blocks += 1;
        }
        prev = b;
    }
    blocks
}

// ── Hand-built cases ──────────────────────────────────────────────────────────

#[test]
fn snv_only_matches_cigar_truth() {
    // 3= 1X 2=  → one SNV, no indels
    let ops = vec![Cigar::Equal(3), Cigar::Diff(1), Cigar::Equal(2)];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(count_snv(&hap), true_snv(&ops));
    assert_eq!(count_snv(&hap), 1);
    assert_eq!(count_continuous_indel_blocks(&hap), 0);
}

#[test]
fn deletion_is_one_indel_block() {
    let ops = vec![Cigar::Equal(3), Cigar::Del(2), Cigar::Equal(3)];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(count_continuous_indel_blocks(&hap), true_indel_blocks(&ops));
    assert_eq!(count_continuous_indel_blocks(&hap), 1);
    assert_eq!(count_snv(&hap), 0);
}

#[test]
fn insertion_is_one_indel_block() {
    let ops = vec![Cigar::Equal(3), Cigar::Ins(2), Cigar::Equal(3)];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(count_continuous_indel_blocks(&hap), true_indel_blocks(&ops));
    assert_eq!(count_continuous_indel_blocks(&hap), 1);
    assert_eq!(count_snv(&hap), 0);
}

#[test]
fn adjacent_ins_del_is_a_single_block() {
    // I then D are consecutive indel ops → one block (vector: marker then -6 run)
    let ops = vec![
        Cigar::Equal(3),
        Cigar::Ins(2),
        Cigar::Del(2),
        Cigar::Equal(3),
    ];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(count_continuous_indel_blocks(&hap), 1);
    assert_eq!(true_indel_blocks(&ops), 1);
}

#[test]
fn compound_ins_then_snv_recovered() {
    // 5= 3I 1X 5=  → golden sums the X's -4 with the 3bp insertion: -4 + 30 = 26 (ends in 6).
    // TRUE: 1 SNV + 1 indel block = 2 variant events — both recovered under golden.
    let ops = vec![
        Cigar::Equal(5),
        Cigar::Ins(3),
        Cigar::Diff(1),
        Cigar::Equal(5),
    ];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();

    assert_eq!(true_snv(&ops), 1);
    assert_eq!(compound_count(&ops), 1);
    // Golden: count_snv == true #SNV (the compound mismatch is recovered, not lost).
    assert_eq!(count_snv(&hap), true_snv(&ops));
    assert_eq!(
        count_snv(&hap),
        1,
        "the compound SNV is recovered under golden summation"
    );
    assert_eq!(count_continuous_indel_blocks(&hap), 1);
    // The true variant count is 2, and golden now reports 2 (the old encoding reported 1).
    assert_eq!(count_var(&hap), 2);
    assert_eq!(true_snv(&ops) + true_indel_blocks(&ops), 2);
}

#[test]
fn trailing_insertion_is_dropped() {
    // 5X 1I → the trailing insertion has no following ref-consuming op to anchor to, so its marker
    // is never flushed and the indel is lost. This mirrors the leading-insertion rule and is
    // faithful to Python's get_hapvector_from_cigar (defer-to-next-ref-op). Real minimap2 --eqx
    // reads end aligned (or soft-clipped), so this never fires in the pipeline.
    let ops = vec![Cigar::Diff(5), Cigar::Ins(1)];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(
        count_continuous_indel_blocks(&hap),
        0,
        "trailing insertion has no anchor → dropped"
    );
    assert_eq!(
        true_indel_blocks(&ops),
        1,
        "op-level truth would be 1; the encoding cannot represent it"
    );
}

#[test]
fn insertion_on_isolated_match_merges_indel_blocks() {
    // 3= 1D 1I 1= 1D 3= : op-level there are TWO indel runs (1D-1I, then 1D) separated by the lone
    // 1= match. But the insertion marker lands on that match, turning it into an indel position, so
    // the two deletions merge into ONE ref-position block. Second information-loss mode of the
    // overwrite encoding; golden does NOT fix it (insertions still attach to a ref position).
    let ops = vec![
        Cigar::Equal(3),
        Cigar::Del(1),
        Cigar::Ins(1),
        Cigar::Equal(1),
        Cigar::Del(1),
        Cigar::Equal(3),
    ];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(
        count_continuous_indel_blocks(&hap),
        1,
        "insertion marker consumes the separating match"
    );
    assert_eq!(ref_position_indel_blocks(&ops), 1);
    assert_eq!(
        true_indel_blocks(&ops),
        2,
        "op-level (biological) count is 2 distinct indel events"
    );
}

#[test]
fn compound_ins_then_longer_snv_recovers_all() {
    // 5= 3I 2X 5= → golden sums only the first mismatch with the insertion (-4+30=26);
    // the second mismatch stays -4. Both are SNVs, so all are recovered.
    let ops = vec![
        Cigar::Equal(5),
        Cigar::Ins(3),
        Cigar::Diff(2),
        Cigar::Equal(5),
    ];
    let hap = extract_hap_vector(&make_record(&ops)).unwrap();
    assert_eq!(true_snv(&ops), 2);
    assert_eq!(compound_count(&ops), 1);
    assert_eq!(count_snv(&hap), 2);
}

// ── Randomized property check (deterministic, dependency-free) ─────────────────

/// Tiny deterministic LCG so the property check is reproducible without a proptest dependency.
struct Lcg(u64);
impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// Generate a "clean" structurally-valid CIGAR: starts with a ref-consuming op (no leading
/// insertion), never stacks two `I` ops back-to-back (which the overwrite encoding would collapse),
/// and always consumes some reference. These constraints match real minimap2 `--eqx` reads.
fn random_clean_cigar(rng: &mut Lcg) -> Vec<Cigar> {
    let n_ops = 1 + rng.below(9) as usize; // 1..=9 ops after the leader
    let mut ops: Vec<Cigar> = Vec::with_capacity(n_ops + 1);
    // leader: Equal or Diff
    let lead_len = 1 + rng.below(6) as u32;
    ops.push(if rng.below(2) == 0 {
        Cigar::Equal(lead_len)
    } else {
        Cigar::Diff(lead_len)
    });
    let mut prev_was_ins = false;
    for _ in 0..n_ops {
        let len = 1 + rng.below(4) as u32;
        // pick a kind, forbidding I directly after I
        let kind = rng.below(4);
        let op = match kind {
            0 => Cigar::Equal(len),
            1 => Cigar::Diff(len),
            2 if !prev_was_ins => Cigar::Ins(len),
            2 => Cigar::Equal(len), // would-be I after I → substitute a match
            _ => Cigar::Del(len),
        };
        prev_was_ins = matches!(op, Cigar::Ins(_));
        ops.push(op);
    }
    // Anchor the end with a match so no insertion is left trailing/unflushed (real --eqx reads end
    // aligned). Without this, a trailing `I` would be silently dropped — see trailing_insertion_is_dropped.
    ops.push(Cigar::Equal(1 + rng.below(3) as u32));
    ops
}

#[test]
fn property_counts_match_cigar_truth_over_random_cigars() {
    let mut rng = Lcg(0x5DEECE66D);
    let iters = 5000;
    let mut total_compounds = 0i32;
    let mut cases_with_compound = 0usize;

    for _ in 0..iters {
        let ops = random_clean_cigar(&mut rng);
        let hap = extract_hap_vector(&make_record(&ops)).unwrap();

        // Indel-block count matches the ref-position structural truth (the space the function
        // operates in), exact even when an insertion merges op-level-distinct runs.
        assert_eq!(
            count_continuous_indel_blocks(&hap),
            ref_position_indel_blocks(&ops),
            "indel-block count diverged from ref-position truth for ops {ops:?}"
        );

        // GOLDEN: SNV count equals the true #SNV exactly — compound (ins-then-mismatch)
        // bases are recovered (`-4 + 10*L` still classifies as an SNV), no undercount.
        let compounds = compound_count(&ops);
        assert_eq!(
            count_snv(&hap),
            true_snv(&ops),
            "snv count != true_snv for ops {ops:?}"
        );

        // count_var is just the sum of its two parts.
        assert_eq!(
            count_var(&hap),
            count_snv(&hap) + count_continuous_indel_blocks(&hap)
        );

        if compounds > 0 {
            cases_with_compound += 1;
            total_compounds += compounds;
        }
    }

    // The generator must actually exercise the compound path, else the recovery assertion is vacuous.
    assert!(
        cases_with_compound > 100,
        "expected many compound cases, got {cases_with_compound} ({total_compounds} compound bases)"
    );
    eprintln!(
        "[cigar_oracle] {iters} random CIGARs: {cases_with_compound} contained compound ins+SNV bases, \
         {total_compounds} compound SNVs recovered by the golden summation encoding."
    );
}
