//! `sdrecall-utils` — shared pure types, errors, logging and parallelism math
//! for the SDrecall Rust pipeline.
//!
//! This is one of the two T0 foundation crates (the other is `sdrecall-io`).
//! It holds plumbing every stage crate reuses so that no stage re-implements it
//! (the no-duplicate-helpers rule). It has **no file I/O** and **no heavy deps**
//! — a math-only stage or a quick test can depend on it without compiling
//! `rust-htslib`.
//!
//! Scope today (grown task-by-task): the crate-wide [`SdError`]/[`Result`], the
//! geometry value types ([`GenomicInterval`], [`RegionKey`], [`Strand`],
//! [`HapId`], [`QnameIdx`]), the parallelism budget split
//! [`configure_parallelism`], and a minimal console logger
//! [`init_console_logger`]. The remaining interface contract (`Paths`) is
//! specified in `docs/analysis/tasks/T0_foundation_crates.DESIGN.md` and lands as
//! the stages that need it arrive.

mod error;
#[macro_use]
mod fatal;
mod geometry;
mod logging;
mod parallel;
mod resources;

pub use error::{Result, SdError};
pub use geometry::{GenomicInterval, HapId, QnameIdx, RegionKey, Strand};
pub use logging::init_console_logger;
pub use parallel::{clamp_threads_u8, configure_parallelism};
pub use resources::{
    graph_memory_units, CpuLease, CpuPhase, IslandResourceGuard, MemoryLease, PhaseResources,
};
