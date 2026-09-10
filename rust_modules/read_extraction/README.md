# rust-read-extraction

Fast Rust-based BAM to FASTQ converter with region filtering for the SDrecall workflow.

**Status:** Production (273 lines, ~5 functions).

## Building

```bash
cd rust_modules
cargo build -p rust_read_extraction
```

## Usage

```rust
let (r1, r2) = rust_read_extraction::bam_to_fastq(
    "aligned.bam",
    "regions.bed",
    "output_R1.fastq",
    "output_R2.fastq",
    false,
    4,
)?;
```

The production `sdrecall` orchestrator calls this library directly. The former
PyO3 wrapper and wheel distribution were removed on 2026-07-17.

For the data-flow diagram + interface + filtering logic see the `read_extraction` appendix in [RUST_MIGRATION_PLAN.md](../../docs/analysis/RUST_MIGRATION_PLAN.md); the function inventory lives in the source.
