# build_phasing_graph

Rust implementation of SDrecall's phasing graph construction, replacing `fp_control/graph_build.py`. Uses `petgraph` for graph algorithms, `rust-htslib` for BAM I/O, and `rayon` for parallel overlap detection.

**Status:** Production (41 functions, 3,226 lines).

## Building

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
export LIBCLANG_PATH=$CONDA_PREFIX/lib

maturin develop --release
```

## Usage

```python
from build_phasing_graph import build_phasing_graph_rust

result = build_phasing_graph_rust(
    bam_file_path, reference_genome, mean_read_length,
    edge_weight_cutoff=0.201, mapq_filter=10, basequal_median_filter=10,
    filter_noisy=True, use_collate=True, threads=4,
)
```

For the data-flow diagram + Python interface contract see [T2 — phasing-graph parity](../../docs/analysis/tasks/T2_phasing_graph_parity.md); the function inventory + struct layout live in the source.
