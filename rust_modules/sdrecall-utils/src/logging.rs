//! Minimal console logger — the transitional replacement for `src/log.py`.
//!
//! Installs a tiny `log::Log` implementation that writes `[LEVEL] message` to
//! stderr. Deliberately dependency-free (no `env_logger`, no color crate) so the
//! foundation library stays light; a richer colored formatter can replace the
//! body later without changing this public entry point.

use log::{LevelFilter, Metadata, Record};

struct ConsoleLogger {
    level: LevelFilter,
}

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            eprintln!("[{}] {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

/// Install the console logger at the given level. Idempotent and best-effort:
/// if a global logger is already set (e.g. a host process installed one), this
/// silently keeps the existing one rather than panicking.
pub fn init_console_logger(level: LevelFilter) {
    let _ = log::set_boxed_logger(Box::new(ConsoleLogger { level }))
        .map(|()| log::set_max_level(level));
}
