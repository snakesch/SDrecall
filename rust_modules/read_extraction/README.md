# rust-read-extraction

Fast Rust-based BAM to FASTQ converter with region filtering for the SDrecall workflow.

**Status:** Production (273 lines, ~5 functions).

## Installation

```bash
pip install rust_read_extraction
```

## Usage

```python
from read_extraction import bam_to_fastq_biobambam

r1_path, r2_path = bam_to_fastq_biobambam(
    input_bam="aligned.bam",
    region_bed="regions.bed",
    output_freads="output_R1.fastq",
    output_rreads="output_R2.fastq",
    multi_aligned=False, threads=4, tmp_dir="/tmp",
)
```

For the data-flow diagram + interface + filtering logic see the `read_extraction` appendix in [RUST_MIGRATION_PLAN.md](../../docs/analysis/RUST_MIGRATION_PLAN.md); the function inventory lives in the source.
