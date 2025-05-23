# SDrecall

SDrecall is a specialized variant caller designed to improve variant detection in segmental duplication (SD) regions where conventional callers often struggle due to mapping ambiguity.

## Overview of the workflow

Preprocessing:
1. Refining the SD pairs by extracting the subsegments in paired SDs where the sequence similarity is sufficiently high to cause mapping ambiguity. (This step is already done and the refined SD pairs are stored in `data/hg19(38)/ref_SD/WGAC.hg19(38).cigar.trim.homo.expanded.highsim.bed.gz`)
2. Given the input BAM file (For now, we recommend using BWA-MEM or BWA-MEM2, because some analysis step relies on the BAM tags XA, and AS/XS. If you are concerned about the computation speed of BWA-MEM, you can try CORA or GPU supported BWA-MEM, which can greatly accelerate BWA-MEM by dozens of folds. Please let us know if you want support of other aligners, the availability of that aligner mostly depend on the BAM tags it has to mark multi-alignment), identify which regions are poorly mapped due to mapping ambiguity. Intersect the region with the user-defined targeting regions (which is by default the whole exome specified in `data/hg19(38)/default_target/hg19(38).func.coding.pad20.bed.gz`), and then identify which SDs are overlapping with the intersecting regions. The overlapping SDs are the final target SDs where SDrecall tries to recover variants from. 
3. Extract all the SD pairs involving the target SDs, and expand the binary pairs to multiplexed network of SDs (composed by 2 types of edges, one marks the homologous relationship and the other marks the physical overlapping relationship between two SD regions). Sensitively identify all the homologous counterparts of the target SDs.

P.S. There are thousands or even tens of thousands of target_SD:homologous_counterparts clusters, if we do realignment for each group, it will be time-consuming. Therefore, we need to identify which clusters can be merged into a single group. In practice, we merge the clusters where the target SDs do not share similar sequences by graph coloring. 

Realignment and variant calling:
1. Recruiting reads from the homologous counterparts and perform realignment. Use BCFtools to call the variants from the realigned reads.
2. Compose the realigned reads into a graph, where each vertex represents a fragment (a pair of reads) and edges connecting the fragments sharing identical sequences within their overlaps. Identify disjoint maximal cliques to phase the realigned reads and assembled into micro-haplotypes
3. Then eliminate less optimal haplotypes from the realignments with Binary Integer Programming.
4. Variant calling with BCFtools based on the filtered realignments.

Post-processing:
1. Merging the result variants with conventional caller variants (Suggested follow up)
2. Annotating common variants using a cohort VCF (optional but suggested)

SDrecall significantly improves small variant (SNVs and small indels) detection in SDs where conventional callers typically miss variants or produce false negatives. 

The result callset is not with high precision rate like the callsets generated by GATK/DeepVariant. SDrecall is primarily designed for molecular diagnosis of Mendelian diseases patients. Despite the limited precision rate, the false positive control measures conducted in SDrecall managed to control the amount of FP noises survived to be causal variant candidates. Upon systematic evaluation, when targeting the entire exome, SDrecall only left 1-3 rare and deleterious FPs to cloud the final selection of the causal variants among candidates while compensated the detection sensitvity within SDs to approximately 95%. 

For molecular diagnosis of Mendelian disease patients, SDrecall provides comprehensive inspection of SD regions that would otherwise be missed by traditional NGS analysis pipelines, while introducing marginal noise that could interfere with causal variant identification downstream. This is the first generalized variant caller able to rescue variants within any regions suffered from mapping ambiguity. Please cite the following paper if you use SDrecall in your research: https://www.researchgate.net/publication/386155530_SDrecall_A_Scalable_Approach_for_Sensitive_Variant_Detection_in_Segmental_Duplications (under review)

## Installation

### Using conda/mamba
Users should first clone this repository to a local directory.

For mamba/conda users, create an environment from YAML:
```bash
mamba env create -f ./env/SDrecall.yml
mamba activate SDrecall
```

### Using docker/singularity
Given the long list of dependencies of SDrecall, we are still working on a docker file / singularity recipe. Any contributions are most welcome.

## Usage

SDrecall provides three main execution modes:

### Complete Pipeline

### With Supplementary VCF and Cohort Annotation (Recommended way to run SDrecall)

```bash
# Run with conventional caller integration and cohort annotation
SDrecall run \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -b /path/to/target.bed \
  -o /path/to/output_dir \
  -t 16 \
  -s <sample_id> \
  --target_tag <label_of_target_region> \
  --conventional_vcf /path/to/deep_variant.vcf \
  --caller_name DeepVariant \
  --cohort_vcf /path/to/control_cohort.vcf \
  --inhouse_common_cutoff 0.01 \
  --cohort_conf_level 0.99
```

### Without Supplementary VCF and Cohort Annotation

```bash
# Run the complete SDrecall pipeline
SDrecall run \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -o /path/to/output_dir \
  -b /path/to/target.bed \
  -t 16 \
  -s <sample_id> \
  --target_tag <label_of_target_region> \
```

### Preparation Only

```bash
# Run only the preparation phase (identifies SD regions, creates masked references)
SDrecall prepare \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -o /path/to/output_dir \
  -b /path/to/target.bed \
  -t 16 \
  -s <sample_id> \
  --target_tag <label_of_target_region> \
  --high_quality_depth 10 \
  --minimum_depth 3
```

### Realignment and Recall Only

```bash
# Run only realignment and recall (requires preparation output)
SDrecall realign \
  -i input.bam \
  -r /path/to/reference.fa \
  -m /path/to/sd_map.bed \
  -b /path/to/target.bed \
  -o /path/to/output_dir \
  -s <sample_id> \
  -t 16 \
  --target_tag <label_of_target_region> \
  --numba_threads 4
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

- **BAM file**: Aligned sequencing reads (must be sorted by coordinates and indexed)
- **Reference genome**: FASTA format (hg19 or hg38 supported)
- **Reference SD map**: BED file with segmental duplication coordinates (Two gzipped bed files are offered in data/hg19(hg38)/ref_SD)
- **Target BED** : Specific regions to analyze (the targeting regions you want to ensure detection sensitivity. For molecular diagnosis of Mendelian diseases, this can be the whole exome, or the coding regions of functionally relevant genes.)
- **Supplementary VCF** (optional but recommended): Conventional caller results to merge with (The VCF file of the same sample, called by other conventional callers like GATK and DeepVariant. If provided, SDrecall will try to merge its own output with this VCF file to offer a final output VCF for downstream analysis)
- **Cohort VCF** (optional but recommended): Population data for identifying common variants ( It is recommended to perform SDrecall on dozens of control samples with the similar coverage profile. Then merge them with bcftools and have AC and AN INFO tags calculated in the final merged VCF. This way, the AN, AC info for each variant called by SDrecall within your inhouse control cohort can be exploited to estimate whether it is truly a common variant in the general population. This is important because traditional population databases like gnomAD and 1000g is based on NGS data, therefore having gaps on the regions like segmental duplications due to the mapping ambiguity)

## Outputs

The main outputs include:
- **Filtered BAM files**: Realigned reads in SD regions (in `<output_dir>/<sample_id>_<assembly>_<target_tag>_SDrecall/recall_results`)
- **Variant calls**: VCF files with variants in SD regions (in `<output_dir>/<sample_id>_<assembly>_<target_tag>_SDrecall/recall_results/<sample_id>.sdrecall.vcf.gz`)
- **Final merged VCF**: Combined results with appropriate filters/tags (in `<output_dir>/final_vcf/`)
- **Log files**: Detailed processing information

## Advanced Features

- **Mapping quality filtering**: Adjust thresholds with `--mq_cutoff` (default: 41)
- **Depth filtering**: Control with `--high_quality_depth` (default: 10, used for pickup multialigned regions, specifies maximal depth of high MAPQ reads to be considered as insufficient coverage for downstream variant calling) and `--minimum_depth` (default: 3, used for pick up the region suffering multialignments, this specifies the minimal required depth regardless of MAPQs)
- **Confidence levels**: Set statistical confidence with `--conf_level` (default: 0.999, used for common variant estimation)
- **Variant filtering**: Filter with `--inhouse_common_cutoff` (default: 0.01) when using cohort data
- **Performance tuning**: Adjust `--threads` for overall parallelism and `--numba_threads` for computational acceleration

## Common Arguments Use SDrecall --help and SDrecall run/prepare/realign --help to see the full argument list

```
-i, --input_bam        Input BAM file path (must be indexed)
-r, --ref_genome       Reference genome path
-m, --reference_sd_map Reference segmental duplication map
-o, --outdir           Output directory
-b, --target_bed       Required target regions for analysis (default: whole exome, offered in data/hg19(38)/default_target)
-s, --sample_id        Sample ID (default: extracted from BAM filename)
-t, --threads          Number of threads to use (default: 10)
-v, --verbose          Verbosity level (INFO, DEBUG, etc.)
--target_tag           Label of the target region (recommended to specify)
```

## Code development and feature requests
SDrecall is under active development. We welcome all kinds of suggestions and collaborations.

## Contact and correspondence
Xingtian Yang (yangyxt@hku.hk), Louis She (snakesch@connect.hku.hk)

