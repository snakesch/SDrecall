# Raw Recall Investigation

## Diagnostic Boundary

`HG002.sdrecall.raw.vcf.gz` is the VCF called directly from the pooled
realigned BAM, before fp-control, phasing, HPSUP annotation, clean-BAM
generation, and final priority merge.

Therefore a raw site-recall gap is a first-half pipeline issue. It should be
debugged in these stages before blaming fp-control or phasing:

- target/RG discovery and final target SD-region construction
- homologous counterpart discovery and grouping
- masked-genome and per-RG BED generation
- read extraction into per-RG FASTQs
- minimap2 realignment against masked RG references
- masked-to-genomic BAM remapping and raw BAM merge/dedup
- bcftools raw variant calling

The current working hypothesis for raw recall below the 98% gate is that Rust
is missing some true variant evidence before fp-control, either because
homologous counterparts are not captured sensitively enough, extracted reads are
not parity-equivalent to Python, realignment remains too ambiguous, or raw
bcftools calling is not seeing support that exists in the realigned BAM.

## 2026-06-24 hg38 Findings

The reindexed hg38 raw benchmark had 717 gold sites and 5,812 raw Rust sites:
622 TP, 95 FN, 5,190 FP, for raw site recall 0.867503.

The 95 raw site FNs were mostly SNVs:

- 90 SNVs
- 3 insertions
- 2 deletions

For the 90 SNV FNs, direct BAM inspection showed two classes:

- 50 loci had alternate support in `HG002.pooled.raw.deduped.bam` that passed the
  same `bcftools mpileup -q 10 -Q 15` thresholds.
- 4 loci had alternate support in the original BAM but no alternate support after
  Rust realignment.
- 36 loci had no alternate support in the original sliced/input BAM at the
  normalized gold allele, so they need separate benchmark/input-region review.

The 50 alt-supported raw FNs are still first-half failures, but not simply
"reads absent from the pooled raw BAM". They are usually low allele-fraction
evidence in the per-RG BAM. For example, `chr11:4288101 A>G` in `RG2` had
`AD=95,17` in `bcftools mpileup`, but `PL=0,2,...`, so `bcftools call -mv`
kept the genotype reference and emitted no variant. Similar examples:
`chr1:120468872 G>A` had `AD=93,11`, and `chr10:47745248 A>G` had `AD=81,9`.

This means the remaining raw-recall gap should be debugged as a first-half
read-recruitment and realignment-balance problem: the alleles often exist, but
the realigned per-RG evidence does not reach a heterozygous call. That is
consistent with either under-recruiting true counterpart-supporting reads or
over-recruiting reference-supporting reads relative to Python.

One concrete first-half parity bug was confirmed and fixed in
`rust_read_extraction`: Rust's extraction predicate did not match the Python
shell fallback.

- Python FC/query extraction applies no MAPQ/tag filter after `samtools view -P
  -L`; Rust was dropping `MAPQ >= 60`.
- Python NFC/counterpart extraction uses `![SA] && ([XA] || mapq < 50)`; Rust
  required `XA && abs(AS - XS) <= 10` and also dropped `MAPQ >= 60`.

On implicated hg38 RGs, this predicate difference is large enough to affect raw
recall. Estimated qname additions with the Python-matching predicate:

- `RG0` FC +71,566 qnames; NFC +258,728 qnames
- `RG1` FC +23,482 qnames; NFC +169,673 qnames
- `RG2` FC +5,153 qnames; NFC +120,349 qnames
- `RG5` FC +3,156 qnames; NFC +106,519 qnames
- `RG11` FC +215 qnames; NFC +30,779 qnames
- `RG66` FC +50 qnames; NFC +16,136 qnames

After this fix, the next benchmark should re-run at least hg38 raw recall before
deeper changes to region discovery or caller settings.
