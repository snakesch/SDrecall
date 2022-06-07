# SD analysis
This repository implements segmental duplication (SD) pipeline from Xingtian. This is supposed to be implemented mainly in python3.

## 0. Input (TODO)
* BAM file ready for variant calling / paired-end FASTQ files
* SD region data fetched from BISER

## 1. Major features
- [x] Fetch SD data from BISER
- [x] Trim CIGAR strings of SD regions
- [ ] Extract homologous regions
- [ ] Gene annotation to BED
- [ ] Filter PID causal gene
- [ ] Prepare map file per region (in GATK4.1)
- [ ] Call polyploidy per priority region (in realign_masked_and_HC_multiploidy.sh)
- [ ] bcftools_concatvcfs (in common_bash_utils)
- [ ] Compare with intrinsic VCFs (py)

## 2. Helper functions
- [x] Check BAM validity
- [x] Check VCF validity
- [ ] Check FASTQ validity
- [ ] Extract XA tags from input BAM
- [ ] Extract XA tags from input FASTQ
- [ ] Set environmental variables (\_THREADS)
- [ ] Prepare masked genome

## 3. Others (TODO)
* Number of threads should be user-defined (for systems without PBS)
* Replace GNU parallel with xargs
