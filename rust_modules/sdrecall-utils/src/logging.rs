//! Minimal console logger — the transitional replacement for `src/log.py`.
//!
//! Installs a tiny `log::Log` implementation that writes `[LEVEL] message` to
//! stderr. Deliberately dependency-free (no `env_logger`, no color crate) so the
//! foundation library stays light; a richer colored formatter can replace the
//! body later without changing this public entry point.

use log::{LevelFilter, Metadata, Record};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

struct ConsoleLogger {
    level: LevelFilter,
    started_at: Instant,
}

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            let unix_ms = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map(|duration| duration.as_millis())
                .unwrap_or(0);
            eprintln!(
                "[ts_unix_ms={unix_ms} elapsed_s={:.6}] [{}] {}",
                self.started_at.elapsed().as_secs_f64(),
                record.level(),
                record.args()
            );
        }
    }

    fn flush(&self) {}
}

/// Install the console logger at the given level. Idempotent and best-effort:
/// if a global logger is already set (e.g. a host process installed one), this
/// silently keeps the existing one rather than panicking.
pub fn init_console_logger(level: LevelFilter) {
    let _ = log::set_boxed_logger(Box::new(ConsoleLogger {
        level,
        started_at: Instant::now(),
    }))
    .map(|()| log::set_max_level(level));
}
