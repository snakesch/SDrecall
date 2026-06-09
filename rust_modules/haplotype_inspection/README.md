# Haplotype Inspection Rust Module

This module provides high-performance Rust implementations for haplotype inspection operations in SDrecall.

## Purpose

Replaces Python NCLS-based haplotype inspection with Rust implementations using:
- `rust-lapper` for fast interval queries (replaces Python NCLS)
- `rust-htslib` for efficient BAM file I/O
- Native Rust for compute-intensive consensus assembly and similarity calculations

## Expected Performance

- BAM → Lapper: 5-8s (vs 19s Python NCLS) - 2.4-3.8× speedup
- Haplotype inspection: 10-15s (vs 60s Python) - 4-6× speedup
- Total per subprocess: 22-30s (vs 82s) - 2.7-3.7× speedup

## Building

```bash
./build.sh
```

Or manually:
```bash
maturin build --release
pip install target/wheels/haplotype_inspection-*.whl
```

## Usage

```python
from haplotype_inspection import inspect_haplotypes_rust

correct_qnames, mismap_qnames = inspect_haplotypes_rust(
    bam_path=bam,
    intrinsic_bam_path=intrinsic_bam,
    hap_qname_info=dict(hap_qname_info),
    qname_hap_info=dict(qname_hap_info),
    # ... other parameters
)
```

## Module Structure

- `lib.rs` - Module entry point and PyO3 bindings
- `bam_lapper.rs` - BAM reading and Lapper construction
- `consensus.rs` - Consensus sequence assembly
- `similarity.rs` - Reference sequence similarity calculations
- `inspection.rs` - Main haplotype inspection logic
- `python_bindings.rs` - Python interface functions
- `structs.rs` - Shared data structures

## Development Status

🚧 Under active development - migrating functions one at a time
