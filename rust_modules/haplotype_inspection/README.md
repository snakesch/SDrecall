# haplotype_inspection

Rust implementation of SDrecall's haplotype inspection pipeline (Phase 2c FP control). Replaces Python NCLS-based haplotype inspection with Rust using `rust-lapper`, `rust-htslib`, and `highs`.

**Status:** Production Rust library (187 passing unit tests, one ignored timing
test, nine separate CIGAR-oracle tests, and five validation harnesses in
`examples/`).

## Building

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

cd rust_modules
CXX=/usr/bin/c++ cargo build --release -p haplotype_inspection
```

## Usage

The production entry point is
`identify_misaligned_haps::inspect_haplotypes`, called in-process by the
`fp-control` crate. The former Python extension API and its dual-path Python
benchmark scripts were removed on 2026-07-17; differential validation is kept
in Rust tests and `examples/` harnesses.

For the data-flow diagram and validation history see
[T1 — validate haplotype-inspection](../../docs/analysis/tasks/T1_validate_haplotype_inspection.md).
