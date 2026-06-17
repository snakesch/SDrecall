# SDrecall Rust — End-to-End Benchmark Plan

> **Goal:** The Rust-only `sdrecall` binary must independently achieve **≥ 93% Recall** and **≥ 45% Precision** on HG002 (GRCh38), comparing recalled variants within target SD regions against the GIAB v4.2.1 gold-standard VCF. Python is **no longer the gold truth**; the Rust pipeline is evaluated against the external benchmark VCF directly.

---

## 1. Test Data Inventory

### 1.1 Input BAM (30× downsampled, SD-region only)

| File | Path | Size | Description |
|------|------|------|-------------|
| HG002.SD.deduped.bam | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam` | ~4.4 GB | 30× WGS reads mapped to SD regions (GRCh38) |

**Preparation method (documented in `benchmarks/README.md`):**
```bash
sambamba slice -L ${REF_SDs} ${RAW_300x_BAM} | \
samtools view --subsample-seed 0 --subsample 0.1 | cut -f 1 | sort | uniq > ${QNAME_LIST}
samtools view -h -u -N ${QNAME_LIST} ${RAW_BAM} | samtools fastq -1 R1.fq -2 R2.fq -0 /dev/null -n
bwa mem -M ... ucsc.hg38.fasta R1.fq R2.fq | samtools sort | gatk MarkDuplicates ...
```

Original 300× BAM: `ftp://ftp-trace.ncbi.nlm.nih.gov/.../HG002.GRCh38.300x.bam`

### 1.2 Gold-Standard VCF

| File | Path |
|------|------|
| GIAB v4.2.1 benchmark (normalized, SD-intersected, lost-ALT removed) | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.norm.bench.no_lost_alt.vcf.gz` |
| GIAB v4.2.1 benchmark (raw normalized) | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.norm.vcf.gz` |
| Benchmark confidence region | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.bed` |

### 1.3 PacBio HiFi Truth BAM (for visual validation)

| File | Path |
|------|------|
| PacBio HiFi Revio 48× | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam` |

### 1.4 SDrecall Intermediate Outputs (Python baseline reference)

| File | Path | Description |
|------|------|-------------|
| Pooled raw deduped BAM | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.raw.deduped.merged.bam` | Merged realigned reads before filtering |
| Pooled clean BAM | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.clean.bam` | Reads passing haplotype inspection |
| Total intrinsic alignments | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/total_intrinsic_alignments.bam` | Identified misaligned reads |
| Called VCF (annotated) | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.gatk.sdrecall.exome.bench.final.vcf.gz` | SDrecall final callset |
| Multi-aligned BED | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.multi_aligned.bench.bed` | Poor coverage / multi-aligned regions |
| Target SD regions | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed` | All SD regions targeted by recall |

### 1.5 Reference Genome

| File | Path |
|------|------|
| GRCh38 FASTA | `/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta` |

---

## 2. Performance Targets

### 2.1 Primary Targets (site-level — ignore zygosity)

| Metric | Minimum | Stretch Goal | Notes |
|--------|---------|--------------|-------|
| **Recall** | ≥ 93% | ≥ 95% | TP / (TP + FN) at CHROM+POS+REF+ALT level |
| **Precision** | ≥ 45% | ≥ 55% | TP / (TP + FP) at CHROM+POS+REF+ALT level |
| **Runtime** | ≤ 30 min | ≤ 15 min | Full pipeline on HG002 SD BAM (single node, 10 CPUs) |
| **Memory** | ≤ 16 GB | ≤ 8 GB | Peak RSS |

### 2.2 Secondary Targets (genotype-level — zygosity must match)

| Metric | Minimum | Stretch Goal | Notes |
|--------|---------|--------------|-------|
| **GT-Recall** | ≥ 88% | ≥ 92% | Same as recall, but het/hom must match gold |
| **GT-Precision** | ≥ 40% | ≥ 50% | Same as precision, but het/hom must match gold |
| **GT-Accuracy** | ≥ 85% | ≥ 90% | Among site-level TPs, fraction with correct zygosity |

### 2.3 Comparison Modes (both run automatically)

The benchmark script always computes **both** modes in one pass:

1. **Site-level** (`mode=site`): match on CHROM + POS + REF + ALT only. This measures whether the pipeline can detect the variant site regardless of dosage.
2. **Genotype-level** (`mode=genotype`): match on CHROM + POS + REF + ALT + GT (het=1 vs hom-alt=2). This measures whether the pipeline correctly distinguishes heterozygous vs homozygous-alt.

Additionally, **GT concordance** reports among site-level TPs:
- `het_to_hom`: gold says het, called says hom-alt (over-calling)
- `hom_to_het`: gold says hom-alt, called says het (under-calling)

---

## 3. IGV Visual Validation

The following IGV command loads all relevant tracks for manual inspection of specific loci:

```bash
view_igv -i "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.raw.deduped.merged.bam,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.clean.bam,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/total_intrinsic_alignments.bam,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.gatk.sdrecall.exome.bench.final.vcf.gz,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.norm.bench.no_lost_alt.vcf.gz,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.multi_aligned.bench.bed,/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed" -r hg38 -t "chr12:31281001-31282000"
```

**Tracks loaded (top → bottom):**
1. `HG002.pooled.raw.deduped.merged.bam` — all recalled reads before haplotype inspection
2. `HG002.pooled.clean.bam` — reads that passed haplotype filtering
3. `HG002.SD.deduped.bam` — original 30× input BAM
4. `HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam` — PacBio truth for visual ground truth
5. `total_intrinsic_alignments.bam` — reads identified as misaligned (should NOT be in clean BAM)
6. `HG002.gatk.sdrecall.exome.bench.final.vcf.gz` — SDrecall final called variants
7. `HG002_GRCh38_1_22_v4.2.1_benchmark.norm.bench.no_lost_alt.vcf.gz` — gold-standard variants
8. `HG002.SD.deduped.multi_aligned.bench.bed` — multi-aligned / poor-coverage regions
9. `all_target_recall_SD_regions.bed` — targeted recall SD regions

**Example locus:** `chr12:31281001-31282000` (useful for spot-checking recall)

---

## 4. Benchmark Workflow

### 4.1 Step 1: Run Rust SDrecall Pipeline

```bash
# Build the full Rust binary
cd rust_modules && cargo build --release --bin sdrecall

# Run end-to-end
./target/release/sdrecall \
  --bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam \
  --ref /paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta \
  --sd-map <SD_MAP_BED> \
  --output-dir /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_HG002/ \
  --threads 8
```

### 4.2 Step 2: Variant Calling on Recalled BAM

```bash
# Call variants with GATK HaplotypeCaller on the clean BAM output
gatk --java-options "-Xmx16G" HaplotypeCaller \
  -R /paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta \
  -I <RUST_CLEAN_BAM> \
  -L <TARGET_SD_REGIONS> \
  -O <RUST_CALLED_VCF> \
  --native-pair-hmm-threads 4
```

### 4.3 Step 3: Intersect with Benchmark Region

```bash
# Restrict both gold and called VCFs to target SD regions overlapping benchmark confidence
bedtools intersect \
  -a /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed \
  -b /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.bed \
  > /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_HG002/benchmark_SD_region.bed

bcftools view -R benchmark_SD_region.bed <GOLD_VCF> | \
bcftools norm -m -both -f <REF> --multi-overlaps 0 -a | \
bcftools norm -d exact | bcftools view -i 'ALT!="*"' -Oz -o gold_targeted.vcf.gz

bcftools view -R benchmark_SD_region.bed <RUST_CALLED_VCF> | \
bcftools norm -m -both -f <REF> --multi-overlaps 0 -a | \
bcftools norm -d exact | bcftools view -i 'ALT!="*"' -Oz -o called_targeted.vcf.gz
```

### 4.4 Step 4: Calculate Precision & Recall

The script runs BOTH comparison modes (site + genotype) in one invocation:

```bash
python3 tests/e2e/calculate_precision_recall.py \
  --gold gold_targeted.vcf.gz \
  --called called_targeted.vcf.gz \
  --output /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_HG002/benchmark_results.tsv \
  --assembly hg38 \
  --caller SDrecall_rust
```

**Output TSV format (two rows per run):**

| sample | assembly | caller | mode | recall | precision | f1 | tp | fn | fp | gold_total | called_total | gt_accuracy | het_to_hom | hom_to_het |
|--------|----------|--------|------|--------|-----------|----|----|----|----|------------|--------------|-------------|------------|------------|
| HG002 | hg38 | SDrecall_rust | site | 0.9345 | 0.4721 | ... | ... | ... | ... | ... | ... | 0.8891 | 12 | 8 |
| HG002 | hg38 | SDrecall_rust | genotype | 0.8912 | 0.4512 | ... | ... | ... | ... | ... | ... | 0.8891 | 12 | 8 |

**Pass criteria (site-level):**
- `recall >= 0.93`
- `precision >= 0.45`

---

## 5. PBS Parallel Benchmark (3 assemblies simultaneously)

### 5.1 PBS Job Configuration

| Parameter | Value | Notes |
|-----------|-------|-------|
| Queue | `medium` | Shared compute nodes |
| CPUs per job | 10 | 10 of 128 per node |
| RAM per job | 200 GB | Max allowed on medium |
| Walltime | 72:00:00 | Max allowed on medium |
| Jobs | 3 | hg19, hg38, T2T (CHM13) |

### 5.2 Submit All 3 Jobs

```bash
bash tests/e2e/submit_e2e_benchmark_pbs.sh
```

This generates 3 PBS scripts and submits them. Each job:
1. Builds the Rust binary (if needed)
2. Runs the full Rust SDrecall pipeline
3. Calls variants (GATK HaplotypeCaller)
4. Intersects with benchmark confidence region
5. Calculates precision/recall (both site + genotype modes)

### 5.3 Monitor Jobs

```bash
# Check job status
qstat -u $(whoami)

# Real-time log streaming (each assembly)
tail -f /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_logs/e2e_hg19_*.log
tail -f /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_logs/e2e_hg38_*.log
tail -f /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_logs/e2e_t2t_*.log
```

### 5.4 Aggregate Results

After all 3 jobs complete:

```bash
head -1 /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38/benchmark_results_hg38.tsv > all_results.tsv
tail -n +2 -q \
  /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg19/benchmark_results_hg19.tsv \
  /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38/benchmark_results_hg38.tsv \
  /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_t2t/benchmark_results_t2t.tsv \
  >> all_results.tsv
column -t -s $'\t' all_results.tsv
```

---

## 6. Existing Benchmark Scripts (Reference)

These scripts represent the Python-era benchmarking and serve as methodology reference:

| Script | Location | Purpose |
|--------|----------|---------|
| `cal_prec_recall_vcf.sh` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/cal_prec_recall_vcf.sh` | Full benchmark pipeline (bash wrapper) |
| `calculate_precision_recall_vcf.py` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/calculate_precision_recall_vcf.py` | Core precision/recall calculator |
| `cal_prec_recall_vcf.py` | `/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/cal_prec_recall_vcf.py` | Legacy SQLite-backed benchmark recorder |
| `benchmark_variants.sh` | `/paedyl01/disk1/yangyxt/SDrecall-rust-migration/benchmarks/benchmark_variants.sh` | In-repo benchmark (annotation + comparison) |
| `cal_prec_recall_vcf.py` | `/paedyl01/disk1/yangyxt/SDrecall-rust-migration/benchmarks/cal_prec_recall_vcf.py` | In-repo legacy precision/recall |

**Key functions in `benchmark_variants.sh`:**
- `bench_callset_per_sample` — full per-sample benchmark (normalize → intersect → annotate → compare)
- `annotate_gnomAD_common` — gnomAD AF annotation
- `annotate_cadd` — CADD deleteriousness annotation
- `annotate_golden_record` — tag TP variants matching gold standard

---

## 7. Test Execution Plan

### Phase A: Module-Level Validation (current — T1)
- Validate `haplotype_inspection` correct/mismap sets match Python on identical inputs
- Use per-island dump hooks in `fp_control/realign_filter_per_cov.py`

### Phase B: Integration Testing (T4 complete)
- Fused Rust Phase 2c pipeline processes islands end-to-end
- Compare clean BAM qnames against Python pipeline output

### Phase C: Full Pipeline Benchmark (T9 — this document)
- Run entire `sdrecall` binary on HG002 30× SD BAM
- Variant-call the output and benchmark against gold VCF
- **Pass/fail gate: recall ≥ 93%, precision ≥ 45%**

### Phase D: Multi-Sample Validation (stretch)
- Repeat on HG003–HG007 once HG002 passes
- Ensure no sample-specific regressions

---

## 8. Directory Layout

```
tests/
├── e2e/
│   ├── END_TO_END_BENCHMARK.md        ← this file
│   ├── run_benchmark_hg002.sh         ← one-shot benchmark runner (single assembly)
│   ├── submit_e2e_benchmark_pbs.sh    ← PBS submission for all 3 assemblies in parallel
│   └── calculate_precision_recall.py  ← precision/recall (site + genotype modes)
└── (existing legacy)
    ├── phase1_python_reference.py
    └── test_e2e_rust_vs_python.py
```

---

## 9. Progress Log

| Date | Event | Result | Notes |
|------|-------|--------|-------|
| 2026-06-15 | Document created | — | Initial benchmark plan for Rust-only evaluation |
| 2026-06-16 | Baseline benchmark collected | DeepVariant: recall=93.6%, prec=34.3% (site); GATK: recall=87.0%, prec=23.3% (site) | Official Python pipeline results on HG002 hg38. See `BASELINE_SUMMARY.txt` |
| | | | |

---

## 10. Notes

- **No Python gold truth:** The Rust pipeline is benchmarked purely against the GIAB gold-standard VCF. Python SDrecall results serve only as a sanity-check reference, not as the target.
- **Lost-ALT filtering:** Gold VCF variants where the ALT allele has zero supporting reads in the 30× input BAM are excluded from recall calculation (already handled in `no_lost_alt` VCF).
- **Normalization:** Both gold and called VCFs are split at multi-allelic sites and left-aligned before comparison.
- **Runtime measurement:** Use `/usr/bin/time -v` to capture wall-clock and peak RSS.
