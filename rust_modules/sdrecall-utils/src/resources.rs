//! Shared CPU and memory leases for concurrent island processing.

use std::fmt;
use std::sync::{Arc, Condvar, Mutex, MutexGuard};

/// A pipeline phase that may borrow CPU capacity beyond its island's base token.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CpuPhase {
    SamtoolsCollate,
    RustPairing,
    GraphBuild,
    Gce,
    Calling,
}

/// Shared phase-aware resource pool for all islands in one pipeline run.
#[derive(Clone)]
pub struct PhaseResources {
    inner: Arc<ResourceInner>,
}

#[derive(Debug)]
struct ResourceConfig {
    total_cpu: usize,
    max_active_islands: usize,
    total_memory_units: usize,
    collate_limit: usize,
}

#[derive(Debug)]
struct ResourceState {
    remaining_islands: usize,
    extra_cpu_in_use: usize,
    memory_in_use: usize,
    collates_in_use: usize,
}

struct ResourceInner {
    config: ResourceConfig,
    state: Mutex<ResourceState>,
    changed: Condvar,
}

impl fmt::Debug for PhaseResources {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PhaseResources")
            .field("config", &self.inner.config)
            .field("state", &*self.lock_state())
            .finish()
    }
}

impl PhaseResources {
    /// Creates a resource pool with independent CPU, graph-memory, and collate limits.
    pub fn new(
        total_cpu: usize,
        max_active_islands: usize,
        remaining_islands: usize,
        total_memory_units: usize,
        collate_limit: usize,
    ) -> Self {
        Self {
            inner: Arc::new(ResourceInner {
                config: ResourceConfig {
                    total_cpu: total_cpu.max(1),
                    max_active_islands: max_active_islands.max(1),
                    total_memory_units: total_memory_units.max(1),
                    collate_limit: collate_limit.max(1),
                },
                state: Mutex::new(ResourceState {
                    remaining_islands,
                    extra_cpu_in_use: 0,
                    memory_in_use: 0,
                    collates_in_use: 0,
                }),
                changed: Condvar::new(),
            }),
        }
    }

    /// Borrows currently idle CPU capacity for one phase without blocking on extra CPUs.
    pub fn acquire_cpu(&self, phase: CpuPhase) -> CpuLease {
        let needs_collate_slot = phase == CpuPhase::SamtoolsCollate;
        let mut state = self.lock_state();
        while needs_collate_slot && state.collates_in_use >= self.inner.config.collate_limit {
            state = self.wait(state);
        }

        let active = state
            .remaining_islands
            .min(self.inner.config.max_active_islands)
            .max(1);
        let fair_threads = self.inner.config.total_cpu.div_ceil(active).max(1);
        let phase_cap = match phase {
            CpuPhase::SamtoolsCollate | CpuPhase::RustPairing | CpuPhase::Calling => 8,
            CpuPhase::GraphBuild | CpuPhase::Gce => self.inner.config.total_cpu,
        };
        let requested_threads = fair_threads.min(phase_cap).max(1);
        let reserved_base = active.min(self.inner.config.total_cpu);
        let extra_capacity = self
            .inner
            .config
            .total_cpu
            .saturating_sub(reserved_base)
            .saturating_sub(state.extra_cpu_in_use);
        let extra_cpu = requested_threads.saturating_sub(1).min(extra_capacity);

        state.extra_cpu_in_use += extra_cpu;
        if needs_collate_slot {
            state.collates_in_use += 1;
        }
        drop(state);

        CpuLease {
            resources: self.clone(),
            phase,
            threads: 1 + extra_cpu,
            extra_cpu,
            collate_slot: needs_collate_slot,
        }
    }

    /// Blocks until the requested graph-memory units can be reserved.
    pub fn acquire_memory(&self, requested_units: usize) -> MemoryLease {
        let units = requested_units
            .max(1)
            .min(self.inner.config.total_memory_units);
        let mut state = self.lock_state();
        while state.memory_in_use + units > self.inner.config.total_memory_units {
            state = self.wait(state);
        }
        state.memory_in_use += units;
        drop(state);

        MemoryLease {
            resources: self.clone(),
            units,
        }
    }

    /// Returns a guard that marks one island complete when dropped.
    pub fn island_guard(&self) -> IslandResourceGuard {
        IslandResourceGuard {
            resources: self.clone(),
            finished: false,
        }
    }

    /// Returns the number of islands that have not reached a terminal outcome.
    pub fn remaining_islands(&self) -> usize {
        self.lock_state().remaining_islands
    }

    fn finish_island(&self) {
        let mut state = self.lock_state();
        state.remaining_islands = state.remaining_islands.saturating_sub(1);
        drop(state);
        self.inner.changed.notify_all();
    }

    fn lock_state(&self) -> MutexGuard<'_, ResourceState> {
        self.inner
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    fn wait<'a>(&self, state: MutexGuard<'a, ResourceState>) -> MutexGuard<'a, ResourceState> {
        self.inner
            .changed
            .wait(state)
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }
}

/// RAII lease for a phase's total thread budget.
#[derive(Debug)]
pub struct CpuLease {
    resources: PhaseResources,
    phase: CpuPhase,
    threads: usize,
    extra_cpu: usize,
    collate_slot: bool,
}

impl CpuLease {
    /// Returns the total threads granted to this phase, including its base token.
    pub fn threads(&self) -> usize {
        self.threads
    }

    /// Returns the phase associated with this lease.
    pub fn phase(&self) -> CpuPhase {
        self.phase
    }
}

impl Drop for CpuLease {
    fn drop(&mut self) {
        let mut state = self.resources.lock_state();
        state.extra_cpu_in_use = state.extra_cpu_in_use.saturating_sub(self.extra_cpu);
        if self.collate_slot {
            state.collates_in_use = state.collates_in_use.saturating_sub(1);
        }
        drop(state);
        self.resources.inner.changed.notify_all();
    }
}

/// RAII lease for weighted graph-memory capacity.
#[derive(Debug)]
pub struct MemoryLease {
    resources: PhaseResources,
    units: usize,
}

impl MemoryLease {
    /// Returns the number of weighted memory units held by this lease.
    pub fn units(&self) -> usize {
        self.units
    }
}

impl Drop for MemoryLease {
    fn drop(&mut self) {
        let mut state = self.resources.lock_state();
        state.memory_in_use = state.memory_in_use.saturating_sub(self.units);
        drop(state);
        self.resources.inner.changed.notify_all();
    }
}

/// Marks an island complete exactly once on explicit finish or scope exit.
#[derive(Debug)]
pub struct IslandResourceGuard {
    resources: PhaseResources,
    finished: bool,
}

impl IslandResourceGuard {
    /// Marks the island complete immediately instead of waiting for scope exit.
    pub fn finish(mut self) {
        self.resources.finish_island();
        self.finished = true;
    }
}

impl Drop for IslandResourceGuard {
    fn drop(&mut self) {
        if !self.finished {
            self.resources.finish_island();
        }
    }
}

/// Estimates weighted graph-memory units from the number of paired entities.
pub fn graph_memory_units(read_pairs: usize) -> usize {
    1 + read_pairs / 12_000
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;
    use std::time::Duration;

    #[test]
    fn cpu_leases_expand_when_the_island_tail_shrinks() {
        let resources = PhaseResources::new(8, 4, 4, 8, 2);
        let leases: Vec<_> = (0..4)
            .map(|_| resources.acquire_cpu(CpuPhase::GraphBuild))
            .collect();
        assert_eq!(leases.iter().map(CpuLease::threads).sum::<usize>(), 8);
        assert!(leases.iter().all(|lease| lease.threads() == 2));
        drop(leases);

        resources.island_guard().finish();
        let tail = resources.acquire_cpu(CpuPhase::GraphBuild);
        assert_eq!(resources.remaining_islands(), 3);
        assert_eq!(tail.threads(), 3);
    }

    #[test]
    fn concurrent_cpu_leases_never_exceed_the_global_budget() {
        let resources = PhaseResources::new(25, 13, 13, 25, 12);
        let leases: Vec<_> = (0..13)
            .map(|_| resources.acquire_cpu(CpuPhase::GraphBuild))
            .collect();

        assert_eq!(leases.iter().map(CpuLease::threads).sum::<usize>(), 25);
        assert_eq!(
            leases.iter().filter(|lease| lease.threads() == 2).count(),
            12
        );
        assert_eq!(
            leases.iter().filter(|lease| lease.threads() == 1).count(),
            1
        );
    }

    #[test]
    fn collate_limit_blocks_excess_processes() {
        let resources = PhaseResources::new(4, 4, 4, 4, 1);
        let first = resources.acquire_cpu(CpuPhase::SamtoolsCollate);
        let worker_resources = resources.clone();
        let (sender, receiver) = mpsc::channel();
        let worker = std::thread::spawn(move || {
            let lease = worker_resources.acquire_cpu(CpuPhase::SamtoolsCollate);
            sender.send(lease.threads()).expect("send lease size");
        });

        assert!(receiver.recv_timeout(Duration::from_millis(50)).is_err());
        drop(first);
        assert_eq!(receiver.recv_timeout(Duration::from_secs(1)).unwrap(), 1);
        worker.join().unwrap();
    }

    #[test]
    fn memory_units_are_weighted_and_released() {
        assert_eq!(graph_memory_units(0), 1);
        assert_eq!(graph_memory_units(11_999), 1);
        assert_eq!(graph_memory_units(12_000), 2);
        assert_eq!(graph_memory_units(179_653), 15);

        let resources = PhaseResources::new(4, 2, 2, 4, 1);
        let first = resources.acquire_memory(4);
        let worker_resources = resources.clone();
        let (sender, receiver) = mpsc::channel();
        let worker = std::thread::spawn(move || {
            let lease = worker_resources.acquire_memory(1);
            sender.send(lease.units()).expect("send memory units");
        });

        assert!(receiver.recv_timeout(Duration::from_millis(50)).is_err());
        drop(first);
        assert_eq!(receiver.recv_timeout(Duration::from_secs(1)).unwrap(), 1);
        worker.join().unwrap();
    }
}
