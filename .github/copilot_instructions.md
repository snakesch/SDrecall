# SDrecall Rust Migration — Copilot Instructions

SDrecall migrates performance-critical Python (Phase 2c FP control, 4595 lines) to Rust for 2.5-3.7x speedup. 46/46 functions ported, 182 tests passing.

For full development instructions see [CLAUDE.md](../CLAUDE.md).
For migration architecture and status see [RUST_MIGRATION_PLAN.md](../docs/analysis/RUST_MIGRATION_PLAN.md).

## Build Environment (SDrecall conda env)

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall

export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

cd rust_modules/haplotype_inspection && maturin build --release
pip install target/wheels/haplotype_inspection-*.whl
```

## Code Quality

- `rustc-hash::FxHashMap` for integer keys, `ahash` for string keys
- Logging via `log` crate (forwarded to Python via `pyo3-log`)
- Unit tests for all functions
- Validate against Python implementation
