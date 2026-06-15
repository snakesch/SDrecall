//! Fatal-invariant enforcement — hard stop with full context when a structural
//! invariant is violated.
//!
//! This is NOT an error-recovery path. A `fatal_invariant!` invocation means:
//! > "The program reached a state that is logically impossible given the code's
//! > own invariants. Continuing would produce silent corruption."
//!
//! **Contract:**
//! 1. Logs the violation at `error!` level with file/line/message.
//! 2. Captures a full `std::backtrace::Backtrace` and logs it (visible when the
//!    log level is set to `trace` or higher — use `RUST_LOG=trace` for post-mortem).
//! 3. Calls `std::process::abort()` — NOT panic. This:
//!    - Does not unwind (so PyO3 cannot catch it and swallow it as a Python exception).
//!    - Produces a core dump (if `ulimit -c` allows) for gdb-level post-mortem.
//!    - Guarantees the process stops *immediately*.
//!
//! ## When to use
//!
//! Use `fatal_invariant!` only for conditions that are **provably impossible**
//! under the module's own logic — where a violation means a bug in this crate
//! (not in user input, not in a file on disk). For recoverable domain errors, use
//! `SdError` + `Result`.

/// Hard-stop the process on a structural invariant violation.
///
/// Logs the violation + backtrace, then aborts (no unwind, no PyO3 catch).
/// Use for conditions that are logically impossible given the code's invariants.
///
/// # Example
///
/// ```ignore
/// let (a, b) = g.edge_endpoints(e)
///     .unwrap_or_else(|| fatal_invariant!(
///         "dijkstra route references edge {:?} that doesn't exist in the graph", e
///     ));
/// ```
#[macro_export]
macro_rules! fatal_invariant {
    ($($arg:tt)*) => {{
        let bt = std::backtrace::Backtrace::force_capture();
        log::error!(
            "FATAL invariant violation at {}:{}:{} — {}",
            file!(),
            line!(),
            column!(),
            format_args!($($arg)*)
        );
        log::error!("Backtrace:\n{bt}");
        std::process::abort();
    }};
}
