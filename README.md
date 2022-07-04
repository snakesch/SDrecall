# Short Variant Calling

This is the prototype of a tool for calling short variants that may be missed by GATK best practice by ascertaining segmental duplication regions. It only supports hg19 build now.

## Prerequisites
* [samtools](http://www.htslib.org/) >=v1.15
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) >=v2.30.0
* [BISER](https://github.com/0xTCG/biser) >=v1.1
* [seqkit](https://github.com/shenwei356/seqkit) >=v2.2.0
* [seqtk](https://github.com/lh3/seqtk) >=v1.3
* [GATK](https://gatk.broadinstitute.org/hc/en-us) >=v4.2.6.1
* [bwa](https://github.com/lh3/bwa) >=v0.7.17
* [GCC](https://gcc.gnu.org/) >=v9.1.0
* [Python](https://www.python.org/downloads/) >=v3.9.2
* [HTSlib](http://www.htslib.org/download/) >=v1.14
* [VCFPy](https://github.com/bihealth/vcfpy) >=v0.13.4
* [mosdepth](https://github.com/brentp/mosdepth) >=v0.3.3
* [bcftools](http://www.htslib.org/download/) >=v1.14

## Installation
```{bash}
# - Clone from a temporary repository until actual release - #
# Install by conda
conda env create -f ./setup/environment.yml

# Install via Singularity
cd setup && singularity build svsd.sif svsd.def # This step may take some time
singularity shell svsd.sif # Create a Singularity shell
cd .. # Users should find main scripts in parent directory of setup/
```

## Input files
### Required
* Base Quality Score Recalibrated (BQSR) BAM file
* Gene annotation file (use NCBI RefSeq data if not specified)
### Optional
* A gene panel in BEDPE format (See [part 0.4](https://github.com/snakesch/shortVariantVCF#04-annotate-and-extract-regions-of-interest))

## Quick run
:crystal_ball: TODO: Quick run wrapper script :crystal_ball:

## Custom run
As users may allocate different number of threads in each step, custom run allows users to execute each step separately with a defined number of threads. Users are advised to reserve more threads for WGS data.

## Contact and correspondance

Xingtian Yang (Email), Louis She (louisshe@hku.hk)
