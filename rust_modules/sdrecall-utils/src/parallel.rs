//! Parallelism budget math — port of `src/utils.py::configure_parallelism`.
//!
//! Splits a total thread budget into `(num_jobs, threads_per_job)` so the
//! pipeline can run `num_jobs` workers each given `threads_per_job` threads.

/// Compute `(num_jobs, threads_per_job)` from a total thread budget.
///
/// Mirrors the Python exactly:
/// ```python
/// num_jobs = np.ceil(total_threads / threads_per_job)
/// return int(num_jobs), int(threads_per_job)
/// ```
/// `num_jobs` is the ceiling of the division; `threads_per_job` is truncated
/// toward zero (matching Python's `int()`).
pub fn configure_parallelism(total_threads: usize, threads_per_job: f64) -> (usize, usize) {
    let num_jobs = (total_threads as f64 / threads_per_job).ceil() as usize;
    (num_jobs, threads_per_job as usize)
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
        // 0 total → 0 jobs
        assert_eq!(configure_parallelism(0, 4.0), (0, 4));
    }
}
