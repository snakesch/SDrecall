//! Parallelism budget math — port of `src/utils.py::configure_parallelism`.
//!
//! Splits a total thread budget into `(num_jobs, threads_per_job)` so the
//! pipeline can run `num_jobs` workers each given `threads_per_job` threads.

/// Compute `(num_jobs, threads_per_job)` from a total thread budget.
///
/// For valid inputs this mirrors the Python:
/// ```python
/// num_jobs = np.ceil(total_threads / threads_per_job)
/// return int(num_jobs), int(threads_per_job)
/// ```
/// `num_jobs` is the ceiling of the division; `threads_per_job` is truncated
/// toward zero (matching Python's `int()`).
///
/// Degenerate budgets are clamped **before** the arithmetic so the result can be
/// handed straight to a rayon pool: a `0` total or a non-finite / `<= 0`
/// `threads_per_job` would otherwise produce `0` jobs (rayon reads
/// `num_threads(0)` as "use every CPU") or `inf as usize == usize::MAX` jobs.
/// Both returned values are therefore guaranteed `>= 1`.
pub fn configure_parallelism(total_threads: usize, threads_per_job: f64) -> (usize, usize) {
    let total_threads = total_threads.max(1);
    let threads_per_job = if threads_per_job.is_finite() && threads_per_job > 0.0 {
        threads_per_job
    } else {
        1.0
    };

    let num_jobs = (total_threads as f64 / threads_per_job).ceil() as usize;
    (num_jobs.max(1), (threads_per_job as usize).max(1))
}

/// Clamp a `usize` thread budget into the `u8` range that the external tools
/// (`samtools -@`, htslib `set_threads`) accept, saturating at `255`.
///
/// A plain `threads as u8` silently wraps — `256 as u8 == 0` would *disable*
/// threading on a high-core machine. Saturation is correct here: a per-tool
/// thread count above 255 is already past any useful point.
pub fn clamp_threads_u8(threads: usize) -> u8 {
    threads.min(u8::MAX as usize) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_python_ceiling() {
        // 16 / 4 = 4.0 → ceil 4 jobs, 4 threads each
        assert_eq!(configure_parallelism(16, 4.0), (4, 4));
        // 17 / 4 = 4.25 → ceil 5 jobs, 4 threads each
        assert_eq!(configure_parallelism(17, 4.0), (5, 4));
        // 1 / 4 = 0.25 → ceil 1 job
        assert_eq!(configure_parallelism(1, 4.0), (1, 4));
        // fractional threads_per_job truncates like Python int()
        assert_eq!(configure_parallelism(10, 2.5), (4, 2));
    }

    #[test]
    fn clamps_degenerate_budgets() {
        // 0 total threads is clamped to 1 (was (0, 4) → rayon num_threads(0) = all CPUs).
        assert_eq!(configure_parallelism(0, 4.0), (1, 4));
        // threads_per_job == 0 must not divide-by-zero into usize::MAX jobs.
        assert_eq!(configure_parallelism(8, 0.0), (8, 1));
        // Negative / non-finite threads_per_job fall back to 1.0.
        assert_eq!(configure_parallelism(8, -3.0), (8, 1));
        assert_eq!(configure_parallelism(8, f64::NAN), (8, 1));
        assert_eq!(configure_parallelism(8, f64::INFINITY), (8, 1));
        // Both outputs are always >= 1.
        let (jobs, tpj) = configure_parallelism(0, 0.0);
        assert!(jobs >= 1 && tpj >= 1);
    }

    #[test]
    fn clamp_threads_saturates_at_255() {
        assert_eq!(clamp_threads_u8(0), 0);
        assert_eq!(clamp_threads_u8(8), 8);
        assert_eq!(clamp_threads_u8(255), 255);
        // The bug this guards: 256 as u8 == 0 (threading silently disabled).
        assert_eq!(clamp_threads_u8(256), 255);
        assert_eq!(clamp_threads_u8(100_000), 255);
    }
}
