# shortVariantVCF
This repository implements short variant calling pipeline from Xingtian.

## 0. Input
* BAM file ready for variant calling / paired-end FASTQ files
* Coordinate data of homologous regions stored in BED format (provided by default)
* tmp_map_file (?)

## 1. Major features
* Prepare map file per region (in GATK4.1)
* Call polyploidy per priority region (in realign_masked_and_HC_multiploidy.sh)
* bcftools_concatvcfs (in common_bash_utils)
* Compare with intrinsic VCFs (py)

## 2. Helper functions / TODO
- [x] Check BAM validity
- [x] Check VCF validity
- [ ] Check FASTQ validity
- [ ] Extract XA tags from input BAM
- [ ] Extract XA tags from input FASTQ
- [ ] Set environmental variables (\_THREADS)
- [ ] Fetch data from SEDEF
- [ ] CIGAR string processing
- [ ] Prepare masked genome

## 3. Others
* Number of threads should be user-defined (for systems without PBS)
* Replace GNU parallel with xargs
