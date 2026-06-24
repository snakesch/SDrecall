# SDrecall Rust Existing-Output Benchmark

This folder contains a lightweight benchmark for completed Rust SDrecall runs.
It does not rerun SDrecall and does not call variants from BAM. It evaluates the
Rust VCFs already produced by the e2e jobs.

## Files

- `calculate_precision_recall.py`: pure VCF comparison. It reports both
  site-level and genotype-level recall/precision in one run.
- `evaluate_existing_rust_outputs.sh`: evaluates the existing hg38, hg19, and
  T2T/CHM13 Rust output folders.

The old wrappers that rebuilt SDrecall, reran the full pipeline, and called
GATK HaplotypeCaller were removed because they no longer match the current
orchestrated Rust pipeline. The Rust pipeline already emits raw, clean, merged,
and final VCFs:

```text
<run_dir>/recall_results/HG002.sdrecall.raw.vcf.gz
<run_dir>/recall_results/HG002.sdrecall.clean.vcf.gz
<run_dir>/recall_results/HG002.sdrecall.merged.vcf.gz
<run_dir>/recall_results/HG002.sdrecall.vcf.gz
```

## Default Inputs

### Rust Outputs

| Assembly | Rust run directory |
| --- | --- |
| hg38 | `/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38/HG002_hg38_exome_SDrecall` |
| hg19 | `/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg19/HG002_hg19_exome_SDrecall` |
| T2T/CHM13 | `/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_t2t/HG002_chm13_exome_SDrecall` |

For each run, the evaluator uses:

```text
recall_results/HG002.sdrecall.<stage>.vcf.gz
realign_groups/all_target_recall_SD_regions.bed
```

`--vcf-stage final` maps to `HG002.sdrecall.vcf.gz`.

### Truth Sets

| Assembly | Gold VCF | Benchmark BED |
| --- | --- | --- |
| hg38 | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.vcf.gz` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.bed` |
| hg19 | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_hg19_v4.2.1_benchmark.vcf.gz` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_hg19_v4.2.1_benchmark.bed` |
| T2T/CHM13 | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.norm.vcf.gz` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed` |

## Method

For each assembly:

1. Intersect the Rust recall regions with the GIAB benchmark/confident BED:

   ```text
   all_target_recall_SD_regions.bed intersect benchmark.bed
   ```

2. Restrict both gold and called VCFs to that evaluation BED.
3. Normalize both VCFs with:

   ```bash
   bcftools norm -m -both -f <reference.fa> --multi-overlaps 0 -a
   bcftools norm -d exact
   bcftools view -i 'ALT!="*"'
   ```

4. Compare normalized VCFs with `calculate_precision_recall.py`.

The calculator skips `MISALIGNED` calls by default and emits two rows per
assembly:

- `site`: match on `CHROM, POS, REF, ALT`
- `genotype`: match on `CHROM, POS, REF, ALT, GT`

It also reports genotype concordance among site-level true positives.

## Run

Activate the SDrecall environment first so `bcftools`, `bedtools`, and `pysam`
are available:

```bash
source ~/.bashrc
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate SDrecall
```

Run all three completed outputs:

```bash
bash tests/e2e/evaluate_existing_rust_outputs.sh
```

Run a subset:

```bash
bash tests/e2e/evaluate_existing_rust_outputs.sh --assemblies hg38,t2t
```

Evaluate raw calls and keep every FILTER value:

```bash
bash tests/e2e/evaluate_existing_rust_outputs.sh \
  --vcf-stage raw \
  --exclude-filters none
```

Write results elsewhere:

```bash
bash tests/e2e/evaluate_existing_rust_outputs.sh \
  --output-dir /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_eval_$(date +%Y%m%d_%H%M%S)
```

Disable FILTER exclusion:

```bash
bash tests/e2e/evaluate_existing_rust_outputs.sh --exclude-filters none
```

## Outputs

Default output root:

```text
/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_eval
```

Per assembly:

```text
<output>/<stage>/<assembly>/rust_recall_x_benchmark.bed
<output>/<stage>/<assembly>/HG002.gold.targeted.norm.vcf.gz
<output>/<stage>/<assembly>/HG002.sdrecall.<stage>.targeted.norm.vcf.gz
<output>/<stage>/<assembly>/benchmark_results_<assembly>_<stage>.tsv
```

Aggregate:

```text
<output>/all_results.tsv
```
