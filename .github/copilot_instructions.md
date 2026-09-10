# SDrecall Rust Migration — Copilot Instructions

SDrecall's production pipeline is implemented as a Rust workspace and top-level
`sdrecall` orchestrator. PyO3 and Maturin bindings were removed on 2026-07-17;
do not reintroduce Python extension boundaries.

For full development instructions see [CLAUDE.md](../CLAUDE.md).
For migration architecture and status see [RUST_MIGRATION_PLAN.md](../docs/analysis/RUST_MIGRATION_PLAN.md).

## Build Environment (SDrecall conda env)

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall

export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

cd rust_modules
CXX=/usr/bin/c++ cargo check --workspace
CXX=/usr/bin/c++ cargo test --workspace
CXX=/usr/bin/c++ cargo clippy --workspace --all-targets -- -D warnings
CXX=/usr/bin/c++ cargo build --release -p sdrecall
```

## Code Quality

- `rustc-hash::FxHashMap` for integer keys, `ahash` for string keys
- Logging via the `log` crate
- Unit tests for all functions
- Preserve differential parity and deterministic VCF/BAM output contracts
