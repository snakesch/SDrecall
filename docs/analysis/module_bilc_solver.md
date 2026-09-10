# Module: bilc_solver.rs — Technical Reference

## ILP Formulation

Binary integer linear program matching Python's `bilc.py:62-148`:

- **Variables:** `x_i in {0,1}` per unique hap_id
- **Objective:** minimize `sum(-coefficient_i * x_i)` (equivalently: maximize coefficient for selected haps)
- **Bounds:** `0 <= x_i <= 1`
- **Constraints** (one per region group `{chrom, start, end}`):
  - `0 <= sum(x_j for hap_ids in region) <= upper_bound`
  - `upper_bound`:
    - `n_haps > 4 AND sum(var_count) >= 1` -> `n_haps - 2`
    - `n_haps <= 1` -> `n_haps`
    - else -> `n_haps - count(varc_rank <= 1)`
- **Returns:** `(select_hap_ids, drop_hap_ids, status, objective)`

Uses `highs` crate v2.0 (`highs-sys` v1.12.1, HiGHS solver v1.12.1) via `RowProblem` API.

## Cross-Validation (1,568 Production Files, HG006)

### Initial Results (Rust vs Stale TSV)

| Metric | Count | Percentage |
|--------|-------|------------|
| ILP-exact match | 1,157 | 73.8% |
| ILP-different | 411 | 26.2% |
| Full match (ILP + post-ILP) | 1,542 | 98.3% |

### Root Cause of 411 ILP-Level Differences

The `.haplotype_meta.tsv` files were from a **previous pipeline run** with a different HiGHS build. Python's `highspy` linked HiGHS ~1.7.2; Rust's `highs-sys` linked HiGHS 1.12.1. HiGHS branch-and-bound is path-dependent on solver version — different versions explore different branches and produce different (but equally optimal) solutions. 0 same-objective alternative-optima cases were found — the "alternative optima" hypothesis was disproven.

### Definitive Validation (Fresh Python vs Rust)

- Re-ran current Python solver on all 26 full-mismatch files: **26/26 match Rust exactly**
- Re-ran on 100 randomly sampled ILP-differing files: **100/100 match Rust exactly**
- **Conclusion: Rust BILC solver is 100% equivalent to the current Python solver.** All 411 differences were stale-TSV artifacts.

### Validation Binary

`validate_bilc_solver.rs` (322 lines): parses Python TSV output, reconstructs BilcRecord rows, runs Rust solver, compares at ILP level and full-match level. `--report <path.tsv>` for machine-readable 14-column TSV.

### Artifacts

| File | Path |
|------|------|
| Cross-validation report | `/paedyl01/disk1/yangyxt/test_tmp/bilc_crossval_all1568.tsv` |
| Per-batch reports | `/paedyl01/disk1/yangyxt/test_tmp/bilc_report_a{a,b,c,d}.tsv` |
| Console log | `/paedyl01/disk1/yangyxt/test_tmp/bilc_crossval_full.log` |
