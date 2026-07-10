//! Shared data structures for haplotype inspection module

/// One row of total_record_df that enters the BILC solver.
/// Contains only the columns used by lp_solve_remained_haplotypes.
#[derive(Clone, Debug)]
pub struct BilcRecord {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub hap_id: i32,
    pub coefficient: f64,
    pub var_count: i32,
    pub varc_rank: i32,
}

/// Region key for grouping BilcRecords.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct RegionKey {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
}

/// Enriched record for BILC pre-processing.
/// Represents one (region, haplotype) row with all metadata columns.
/// Corresponds to one row of Python's total_record_df after enrichment.
#[derive(Clone, Debug)]
pub struct EnrichedRecord {
    // Core 8 columns from IdentifyMisalignmentResult.record_2d_arr
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub total_depth: i32,
    pub hap_id: i32,
    pub hap_depth: i32,
    pub var_count: i32,
    pub indel_count: i32,
    pub psv_count: i32,
    // Enrichment columns from per-haplotype accumulators
    pub extreme_vard: bool,
    pub scatter_hap: bool,
    pub hap_var_count: i32,
    pub hap_max_sim_scores: f64, // rounded to 0.1
    pub hap_max_psvs: i32,
    // Computed during coefficient calculation
    pub varc_rank: i32, // sim-score rank (Python "varc_rank"); read by the BILC solver
    pub rank: i32,      // coefficient rank (Python "rank" column); TSV/debug only
    pub interval_coefficient: f64,
    pub coefficient: f64,
}
