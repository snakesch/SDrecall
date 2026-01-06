# rust-read-extraction

A fast Rust-based BAM to FASTQ converter with region filtering, designed for the [SDrecall](https://github.com/snakesch/SDrecall) workflow.

## Features

- **Fast BAM to FASTQ conversion**: Efficiently extracts paired-end reads from BAM files
- **Region-based filtering**: Extracts only reads overlapping specified BED regions
- **Multi-aligned read filtering**: Optional filtering for multi-aligned reads based on alignment tags (SA, XA, AS, XS)
- **Paired-end aware**: Ensures both mates are extracted when one overlaps the target region

## Installation

```bash
pip install rust_read_extraction
```

## Usage

```python
from rust_read_extraction import bam_to_fastq_biobambam

# Extract reads from specific regions
r1_path, r2_path = bam_to_fastq_biobambam(
    input_bam="aligned.bam",
    region_bed="regions.bed",
    output_freads="output_R1.fastq",
    output_rreads="output_R2.fastq",
    multi_aligned=False,  # Set True to filter multi-aligned reads
    threads=4,
    tmp_dir="/tmp"
)
```

## Parameters

- `input_bam`: Path to indexed BAM file
- `region_bed`: BED file specifying regions of interest
- `output_freads`: Output path for R1 FASTQ
- `output_rreads`: Output path for R2 FASTQ
- `multi_aligned`: Whether to apply multi-aligned read filtering (default: False)
- `threads`: Number of threads for BAM reading (default: 1)
- `tmp_dir`: Temporary directory (default: "/tmp")

## Multi-aligned Read Filtering

When `multi_aligned=True`, the module filters reads based on:
- MAPQ < 60 (mandatory)
- No SA tag (no supplementary alignments)
- Has XA tag (has secondary alignments)
- |AS - XS| <= 10 (alignment score difference)

## License

MIT License

## Author

Xingtian Yang (yangyxt@hku.hk)

