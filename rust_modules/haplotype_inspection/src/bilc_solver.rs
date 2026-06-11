/// BILC (Binary Integer Linear Constraint) solver for haplotype selection.
///
/// Port of `fp_control/bilc.py::lp_solve_remained_haplotypes`.
///
/// Uses the HiGHS solver (via the `highs` crate) to maximize the sum of
/// haplotype coefficients while satisfying per-region constraints on how many
/// haplotypes may be flagged as misaligned.
use std::fmt;

use highs::{HighsModelStatus, RowProblem, Sense};
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::structs::{BilcRecord, RegionKey};

/// Status returned by the BILC solver.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BilcStatus {
    Optimal,
    Infeasible,
    Other(String),
}

impl fmt::Display for BilcStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BilcStatus::Optimal => write!(f, "Optimal"),
            BilcStatus::Infeasible => write!(f, "Infeasible"),
            BilcStatus::Other(s) => write!(f, "{s}"),
        }
    }
}

impl From<HighsModelStatus> for BilcStatus {
    fn from(status: HighsModelStatus) -> Self {
        match status {
            HighsModelStatus::Optimal => BilcStatus::Optimal,
            HighsModelStatus::Infeasible => BilcStatus::Infeasible,
            other => BilcStatus::Other(format!("{other:?}")),
        }
    }
}

/// Solve the BILC problem: determine which haplotypes to keep vs drop.
///
/// Each unique hap_id becomes a binary variable (0 = keep, 1 = drop).
/// The objective minimizes the negative coefficient sum (i.e. maximizes
/// the total coefficient of dropped haplotypes — haplotypes with higher
/// misalignment likelihood get dropped first).
///
/// Per-region constraints limit how many haplotypes can be dropped
/// in each genomic region, based on the total variant count and
/// rank distribution.
///
/// # Arguments
/// * `records` — one row per (region, haplotype) pair, may repeat hap_id
///   across multiple regions.
///
/// # Returns
/// * `(select_hap_ids, drop_hap_ids, status)` where select = kept (solution 0)
///   and drop = flagged misaligned (solution 1).
pub fn lp_solve_remained_haplotypes(
    records: &[BilcRecord],
) -> (FxHashSet<i32>, FxHashSet<i32>, BilcStatus) {
    let (sel, drop, status, _obj) = lp_solve_remained_haplotypes_with_obj(records);
    (sel, drop, status)
}

/// Same as [`lp_solve_remained_haplotypes`] but also returns the objective value
/// (sum of coefficients of dropped haplotypes). Useful for cross-validation.
pub fn lp_solve_remained_haplotypes_with_obj(
    records: &[BilcRecord],
) -> (FxHashSet<i32>, FxHashSet<i32>, BilcStatus, f64) {
    // --- Step 1: Extract unique (hap_id, coefficient) pairs ---
    // Use FxHashMap to preserve first-seen coefficient per hap_id,
    // matching the Python drop_duplicates() which keeps the first occurrence.
    let mut seen = FxHashMap::default();
    for rec in records {
        seen.entry(rec.hap_id).or_insert(rec.coefficient);
    }
    // Deterministic ordering: sort by hap_id for reproducibility.
    let mut hap_coefficients: Vec<(i32, f64)> = seen.into_iter().collect();
    hap_coefficients.sort_by_key(|&(hid, _)| hid);

    let hap_no = hap_coefficients.len();
    if hap_no == 0 {
        return (FxHashSet::default(), FxHashSet::default(), BilcStatus::Optimal, 0.0);
    }

    // Build bidirectional index <-> hap_id mappings.
    let index_to_hapid: Vec<i32> = hap_coefficients.iter().map(|&(hid, _)| hid).collect();
    let hapid_to_index: FxHashMap<i32, usize> = index_to_hapid
        .iter()
        .enumerate()
        .map(|(idx, &hid)| (hid, idx))
        .collect();

    info!(
        "BILC: {} unique haplotypes, index_to_hapid: {:?}",
        hap_no,
        &index_to_hapid[..std::cmp::min(10, hap_no)]
    );

    // --- Step 2: Group records by (chrom, start, end) ---
    let mut region_groups: Vec<(RegionKey, Vec<&BilcRecord>)> = Vec::new();
    {
        // Use an index map to keep insertion order.
        let mut key_to_idx: FxHashMap<RegionKey, usize> = FxHashMap::default();
        for rec in records {
            let key = RegionKey {
                chrom: rec.chrom.clone(),
                start: rec.start,
                end: rec.end,
            };
            if let Some(&idx) = key_to_idx.get(&key) {
                region_groups[idx].1.push(rec);
            } else {
                let idx = region_groups.len();
                key_to_idx.insert(key.clone(), idx);
                region_groups.push((key, vec![rec]));
            }
        }
    }
    let inspect_region_no = region_groups.len();

    info!(
        "BILC: constraint matrix shape ({inspect_region_no}, {hap_no})"
    );

    // --- Step 3: Build the HiGHS model ---
    let mut pb = RowProblem::default();

    // Add binary variables: one per unique haplotype.
    // Objective = minimize(-coefficient), i.e. maximize coefficient for selected=1 haps.
    // Col handles are stored in order matching index_to_hapid.
    let cols: Vec<_> = hap_coefficients
        .iter()
        .enumerate()
        .map(|(i, &(_hid, coeff))| {
            let cost = -coeff; // minimize negative = maximize
            let col = pb.add_integer_column(cost, 0.0..=1.0);
            if i < 10 {
                info!(
                    "BILC: var {} hap_id={} cost={}",
                    i, index_to_hapid[i], cost
                );
            }
            col
        })
        .collect();

    // --- Step 4: Add per-region constraints ---
    for (key, group) in &region_groups {
        let region_str = format!("{}:{}-{}", key.chrom, key.start, key.end);

        // Unique hap_ids in this region.
        let mut included_hids: FxHashSet<i32> = FxHashSet::default();
        let mut total_var_count: i32 = 0;
        let mut rank_1_hids: FxHashSet<i32> = FxHashSet::default();
        let mut rank_2_hids: FxHashSet<i32> = FxHashSet::default();

        for rec in group {
            included_hids.insert(rec.hap_id);
            total_var_count += rec.var_count;
            if rec.varc_rank <= 1 {
                rank_1_hids.insert(rec.hap_id);
            }
            if rec.varc_rank <= 2 {
                rank_2_hids.insert(rec.hap_id);
            }
        }

        let n_haps = included_hids.len();
        let rank_1_count = rank_1_hids.len();
        let _rank_2_count = rank_2_hids.len();

        // Upper bound: how many haplotypes can be flagged as misaligned.
        let upper_bound = if n_haps <= 1 {
            n_haps
        } else if n_haps > 4 && total_var_count >= 1 {
            n_haps - 2
        } else {
            n_haps - rank_1_count
        };

        // Build sparse row: (col_handle, 1.0) for each hap in this region.
        let row_factors: Vec<_> = included_hids
            .iter()
            .map(|hid| {
                let idx = hapid_to_index[hid];
                (cols[idx], 1.0)
            })
            .collect();

        // Constraint: 0 <= sum(x_i for i in region_haps) <= upper_bound
        pb.add_row(0.0..=(upper_bound as f64), &row_factors);

        info!(
            "BILC: region {}, {} haps {:?}, indices {:?}, upper_bound={}",
            region_str,
            n_haps,
            included_hids.iter().collect::<Vec<_>>(),
            included_hids
                .iter()
                .map(|hid| hapid_to_index[hid])
                .collect::<Vec<_>>(),
            upper_bound,
        );
    }

    // --- Step 5: Solve ---
    let mut model = pb.optimise(Sense::Minimise);
    model.make_quiet();
    let solved = model.solve();

    let status = solved.status();
    info!("BILC: model status = {status:?}");

    let bilc_status = BilcStatus::from(status);

    // --- Step 6: Extract solution ---
    let solution = solved.get_solution();
    let col_values = solution.columns();

    let mut select_hap_ids = FxHashSet::default();
    let mut drop_hap_ids = FxHashSet::default();
    let mut objective_value: f64 = 0.0;

    for (i, &val) in col_values.iter().enumerate() {
        let hid = index_to_hapid[i];
        let coeff = hap_coefficients[i].1;
        // Binary variable: 0 = keep, 1 = drop.
        // Use threshold 0.5 for robustness against floating-point noise.
        if val < 0.5 {
            select_hap_ids.insert(hid);
        } else {
            drop_hap_ids.insert(hid);
            objective_value += coeff;
        }
    }

    info!(
        "BILC: select_hap_ids={select_hap_ids:?}, drop_hap_ids={drop_hap_ids:?}, obj_value={objective_value}"
    );

    (select_hap_ids, drop_hap_ids, bilc_status, objective_value)
}

// =============================================================================
// Tests
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a BilcRecord with defaults for fields not under test.
    fn rec(chrom: &str, start: i32, end: i32, hap_id: i32, coeff: f64, var_count: i32, varc_rank: i32) -> BilcRecord {
        BilcRecord {
            chrom: chrom.to_string(),
            start,
            end,
            hap_id,
            coefficient: coeff,
            var_count,
            varc_rank,
        }
    }

    #[test]
    fn test_empty_records() {
        let (sel, drop, status) = lp_solve_remained_haplotypes(&[]);
        assert_eq!(status, BilcStatus::Optimal);
        assert!(sel.is_empty());
        assert!(drop.is_empty());
    }

    #[test]
    fn test_single_haplotype_single_region() {
        // One hap in one region — upper_bound = n_haps = 1 (since n_haps <= 1).
        // The solver can assign 0 or 1 freely; with negative cost it prefers 1 (drop).
        let records = vec![rec("chr1", 1000, 2000, 10, 5.0, 3, 1)];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        // cost = -5.0; minimizing => wants x=1 (drop it)
        assert!(drop.contains(&10));
        assert!(!sel.contains(&10));
    }

    #[test]
    fn test_two_haps_high_vs_low_coefficient() {
        // Two haps in one region. n_haps=2, var_count sum > 0 but n_haps <= 4,
        // so upper_bound = n_haps - rank_1_count.
        // hap 1: coeff=10, varc_rank=2 (high misalignment score → should be dropped)
        // hap 2: coeff=1,  varc_rank=1 (low score → should be kept)
        // upper_bound = 2 - 1 (rank_1_count=1) = 1 → at most 1 can be dropped
        let records = vec![
            rec("chr1", 1000, 2000, 1, 10.0, 3, 2),
            rec("chr1", 1000, 2000, 2, 1.0, 2, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        // Solver prefers dropping hap 1 (coeff=10 → cost=-10, bigger payoff)
        assert!(drop.contains(&1));
        assert!(sel.contains(&2));
    }

    #[test]
    fn test_five_haps_with_var_count_constraint() {
        // 5 haps, one region. n_haps=5 > 4 and total_var_count >= 1,
        // so upper_bound = 5 - 2 = 3. At most 3 can be dropped.
        let records = vec![
            rec("chr1", 100, 200, 0, 20.0, 1, 3),
            rec("chr1", 100, 200, 1, 15.0, 1, 3),
            rec("chr1", 100, 200, 2, 10.0, 0, 2),
            rec("chr1", 100, 200, 3, 2.0, 0, 1),
            rec("chr1", 100, 200, 4, 1.0, 0, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        assert_eq!(drop.len(), 3, "should drop exactly 3 (upper_bound = 3)");
        // The top-3 by coefficient (0, 1, 2) should be dropped.
        assert!(drop.contains(&0));
        assert!(drop.contains(&1));
        assert!(drop.contains(&2));
        assert!(sel.contains(&3));
        assert!(sel.contains(&4));
    }

    #[test]
    fn test_hap_appears_in_multiple_regions() {
        // Hap 1 appears in region A and B. Hap 2 only in A, hap 3 only in B.
        // Region A: haps {1, 2}, n=2 <= 4, rank_1_count=1(hap2), upper_bound = 2-1 = 1
        // Region B: haps {1, 3}, n=2 <= 4, rank_1_count=1(hap3), upper_bound = 2-1 = 1
        // Hap 1 has highest coefficient → solver wants to drop it.
        // Dropping hap 1 uses 1 slot in both regions, so it's feasible.
        let records = vec![
            rec("chr1", 100, 200, 1, 8.0, 1, 2),
            rec("chr1", 100, 200, 2, 1.0, 1, 1),
            rec("chr1", 300, 400, 1, 8.0, 1, 2),
            rec("chr1", 300, 400, 3, 1.0, 1, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        assert!(drop.contains(&1));
        assert!(sel.contains(&2));
        assert!(sel.contains(&3));
    }

    #[test]
    fn test_multi_region_constraint_interaction() {
        // Haps 10, 20, 30 share region A (5 haps total with 40, 50).
        // Haps 10, 20 are also in region B (just 2 haps).
        // Region A: 5 haps, var_count=5 >= 1, upper_bound = 5-2 = 3
        // Region B: 2 haps, var_count=0 < 1, rank_1_count = 1(hap20), upper_bound = 2-1 = 1
        //   → at most 1 of {10, 20} can be dropped in region B
        // Hap 10 has the highest coefficient in B, so it gets dropped there.
        let records = vec![
            // Region A
            rec("chr1", 100, 200, 10, 9.0, 1, 3),
            rec("chr1", 100, 200, 20, 7.0, 1, 1),
            rec("chr1", 100, 200, 30, 6.0, 1, 2),
            rec("chr1", 100, 200, 40, 3.0, 1, 2),
            rec("chr1", 100, 200, 50, 1.0, 1, 1),
            // Region B
            rec("chr2", 500, 600, 10, 9.0, 0, 2),
            rec("chr2", 500, 600, 20, 7.0, 0, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        // In region B, at most 1 drop → hap 10 (coeff=9) gets dropped
        assert!(drop.contains(&10));
        // The solver can also drop up to 2 more from region A.
        // Total drops should respect both constraints.
        assert!(drop.len() <= 4); // at most 3 from A, but hap 10 is already counted
    }

    #[test]
    fn test_all_rank_1_no_drops() {
        // Two haps, both rank 1. upper_bound = 2 - 2 = 0.
        // No drops allowed → both must be kept.
        let records = vec![
            rec("chr1", 100, 200, 1, 10.0, 1, 1),
            rec("chr1", 100, 200, 2, 8.0, 1, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        assert!(drop.is_empty());
        assert!(sel.contains(&1));
        assert!(sel.contains(&2));
    }

    #[test]
    fn test_zero_var_count_rank_based_bound() {
        // 3 haps, all var_count=0. n_haps=3 <= 4.
        // rank_1_count: hap 3 has rank 1. So upper_bound = 3 - 1 = 2.
        let records = vec![
            rec("chr1", 100, 200, 1, 10.0, 0, 3),
            rec("chr1", 100, 200, 2, 5.0, 0, 2),
            rec("chr1", 100, 200, 3, 1.0, 0, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        // Drops top 2 by coefficient (hap 1 and hap 2).
        assert_eq!(drop.len(), 2);
        assert!(drop.contains(&1));
        assert!(drop.contains(&2));
        assert!(sel.contains(&3));
    }

    #[test]
    fn test_bilc_status_display() {
        assert_eq!(format!("{}", BilcStatus::Optimal), "Optimal");
        assert_eq!(format!("{}", BilcStatus::Infeasible), "Infeasible");
        assert_eq!(
            format!("{}", BilcStatus::Other("Foo".into())),
            "Foo"
        );
    }

    #[test]
    fn test_duplicate_hap_id_uses_first_coefficient() {
        // hap_id 1 appears with coeff 10 first and 99 second.
        // FxHashMap::or_insert keeps the first ⇒ coeff = 10.
        let records = vec![
            rec("chr1", 100, 200, 1, 10.0, 1, 2),
            rec("chr1", 100, 200, 1, 99.0, 1, 2),
            rec("chr1", 100, 200, 2, 1.0, 1, 1),
        ];
        let (sel, drop, status) = lp_solve_remained_haplotypes(&records);
        assert_eq!(status, BilcStatus::Optimal);
        // With coeff=10, hap 1 should still be preferentially dropped.
        assert!(drop.contains(&1));
        assert!(sel.contains(&2));
    }
}
