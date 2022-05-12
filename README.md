# shortVariantVCF
This repository implements short variant calling pipeline from Xingtian.

## 0. Input
* BAM file ready for variant calling
* Coordinate data of homologous regions stored in BED format (provided by default)
* tmp_map_file (?)

## 1. Major features
* Prepare map file per region (in GATK4.1)
* Call polyploidy per priority region (in realign_masked_and_HC_multiploidy.sh)
* bcftools_concatvcfs (in common_bash_utils)
* Compare with intrinsic VCFs (py)

## 2. Helper functions
- [x] Check BAM validity
- [x] Check VCF validity

## 3. Others
* Number of threads should be user-defined (for systems without PBS)
* Replace GNU parallel with xargs
