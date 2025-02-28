# SDrecall

SDrecall is a specialized variant caller designed to improve variant detection in segmental duplication (SD) regions where conventional callers often struggle due to mapping ambiguity.

## Overview

SDrecall works by:
1. Identifying problematic SD regions with mapping ambiguity
2. Creating masked reference genomes for these regions
3. Realigning reads and performing targeted variant calling
4. Merging the results with conventional caller output (optional)
5. Annotating common variants using a cohort VCF (optional)

This approach significantly improves variant detection in complex genomic regions where conventional callers typically miss variants or produce false negatives.

## Installation

### Using conda/mamba
Users should first clone this repository to a local directory.

For mamba users, create an environment from YAML:
```bash
mamba env create -f ./env/SDrecall.yml
mamba activate SDrecall
```

### Using docker/singularity
Given the long list of dependencies of SDrecall, we are still working on a docker file / singularity recipe. Any contributions are most welcome.

## Usage

SDrecall provides three main execution modes:

### Complete Pipeline

```bash
# Run the complete SDrecall pipeline
python SDrecall.py run \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -o /path/to/output_dir \
  -t 16
```

### Preparation Only

```bash
# Run only the preparation phase (identifies SD regions, creates masked references)
python SDrecall.py prepare \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -o /path/to/output_dir \
  -t 16 \
  --high_quality_depth 10 \
  --minimum_depth 3
```

### Realignment and Recall Only

```bash
# Run only realignment and recall (requires preparation output)
python SDrecall.py realign \
  -i input.bam \
  -r /path/to/reference.fa \
  -o /path/to/output_dir \
  -t 16 \
  --numba_threads 4
```

### With Supplementary VCF and Cohort Annotation

```bash
# Run with conventional caller integration and cohort annotation
python SDrecall.py run \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -o /path/to/output_dir \
  -t 16 \
  --supplementary_vcf /path/to/deep_variant.vcf \
  --caller_name DeepVariant \
  --cohort_vcf /path/to/control_cohort.vcf \
  --inhouse_common_cutoff 0.05 \
  --conf_level 0.999
```

## Workflow Stages

### 1. Preparation (`prepare_recall_regions.py`)

This stage identifies SD regions with mapping issues:
- Extracts multi-aligned regions based on mapping quality (`pick_multialigned_regions()`)
- Compares to a reference SD map
- Creates a multiplex graph with SD pairs (`build_SD_graph()`)
- Builds masked reference genomes for each SD group (`build_beds_and_masked_genomes()`)

### 2. Realignment and Recall (`realign_and_recall.py`)

This stage performs targeted variant calling:
- Extracts reads from identified SD regions (`imap_prepare_masked_align_region_per_RG()`)
- Realigns to masked references (`imap_process_masked_bam()`)
- Eliminates misalignments (`eliminate_misalignments()`)
- Performs variant calling on filtered alignments
- Tags variants for provenance

### 3. Post-processing (`post_process_vcf()` in `SDrecall.py`)

Final steps may include:
- Annotating variants with cohort data (`identify_common_vars.py`)
- Merging with conventional caller output (`merge_with_priority()` in `src/merge_variants_with_priority.py`)
- Prioritizing variants based on quality metrics

## Inputs

- **BAM file**: Aligned sequencing reads (must be indexed)
- **Reference genome**: FASTA format (hg19 or hg38 supported)
- **Reference SD map**: BED file with segmental duplication coordinates
- **Target BED** (optional): Specific regions to analyze
- **Supplementary VCF** (optional): Conventional caller results to merge with
- **Cohort VCF** (optional): Population data for identifying common variants

## Outputs

The main outputs include:
- **Filtered BAM files**: Realigned reads in SD regions (in `<output_dir>/realigned/`)
- **Variant calls**: VCF files with variants in SD regions (in `<output_dir>/vcf/`)
- **Final merged VCF**: Combined results with appropriate filters/tags (in `<output_dir>/final_vcf/`)
- **Log files**: Detailed processing information

## Advanced Features

- **Mapping quality filtering**: Adjust thresholds with `--mq_cutoff` (default: 41)
- **Depth filtering**: Control with `--high_quality_depth` (default: 10) and `--minimum_depth` (default: 3)
- **Confidence levels**: Set statistical confidence with `--conf_level`
- **Variant filtering**: Filter with `--inhouse_common_cutoff` (default: 0.05) when using cohort data
- **Performance tuning**: Adjust `--threads` for overall parallelism and `--numba_threads` for computational acceleration

## Common Arguments

```
-i, --input_bam        Input BAM file path (must be indexed)
-r, --ref_genome       Reference genome path
-m, --reference_sd_map Reference segmental duplication map
-o, --outdir           Output directory
-b, --target_bed       Optional target regions for analysis (default: whole genome)
-s, --sample_id        Sample ID (default: extracted from BAM filename)
-t, --threads          Number of threads to use (default: 10)
-v, --verbose          Verbosity level (INFO, DEBUG, etc.)
```

## Code development and feature requests
SDrecall is under active development. We welcome all kinds of suggestions and collaborations.

## Contact and correspondence
Xingtian Yang (yangyxt@hku.hk), Louis She (snakesch@connect.hku.hk)
