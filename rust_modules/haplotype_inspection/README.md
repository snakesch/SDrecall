# haplotype_inspection

Rust implementation of SDrecall's haplotype inspection pipeline (Phase 2c FP control). Replaces Python NCLS-based haplotype inspection with Rust using `rust-lapper`, `rust-htslib`, and `highs`.

**Status:** Code complete (54 functions, 182 unit tests — 181 active + 1 ignored bench, 5 example harnesses in `examples/`).

## Building

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

maturin build --release
pip install target/wheels/haplotype_inspection-*.whl
```

## Usage

```python
from haplotype_inspection import inspect_haplotypes_rust

correct_qnames, mismap_qnames = inspect_haplotypes_rust(
    bam_path, intrinsic_bam_path, hap_qname_info, qname_hap_info,
    qname_to_node, total_lowqual_qnames, compare_haplotype_meta_tab,
    mean_read_length, recall_mq_cutoff, basequal_median_cutoff,
)
```

For the data-flow diagram + Python interface contract see [T1 — validate haplotype-inspection](../../docs/analysis/tasks/T1_validate_haplotype_inspection.md), and the validated submodule deep-dives in [docs/analysis/](../../docs/analysis/); the function inventory + struct layout live in the source.
