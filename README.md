# SDrecall

SDrecall is the prototype of a tool for calling short variants that may be missed by GATK best practice by ascertaining segmental duplication regions. This pipeline is implemented in python3. 

## Prerequisites
* [samtools](http://www.htslib.org/) v1.15
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) v2.30.0
* [BISER](https://github.com/0xTCG/biser) v1.1
* [seqkit](https://github.com/shenwei356/seqkit) v2.2.0
* [seqtk](https://github.com/lh3/seqtk) v1.3
* [GATK](https://gatk.broadinstitute.org/hc/en-us) v4.2.6.1
* [bwa](https://github.com/lh3/bwa) v0.7.17
* [Python](https://www.python.org/downloads/) v3.9.2
* [HTSlib](http://www.htslib.org/download/) v1.14
* [mosdepth](https://github.com/brentp/mosdepth) v0.3.3
* [bcftools](http://www.htslib.org/download/) v1.14
* [pandarallel](https://pypi.org/project/pandarallel/) v1.6.3 (required in step 4)

Note: We do not guarantee compatibility for softwares of more updated versions.

## Installation
For conda users, create an environment from YAML.
```{bash}
conda env create -f ./setup/environment.yml
conda activate SDrecall
```
For [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html) users, please build a container with the provided image file and set up accordingly.

```{bash}
cd setup 
singularity build svsd.sif svsd.def # This step may take some time
singularity shell svsd.sif # Create a Singularity shell
conda init && source ~/.bashrc
conda activate SDrecall
```

## Input files
### Required
* Base Quality Score Recalibrated (BQSR) BAM file
* A panel of genes of interest ([part 0.4](https://github.com/snakesch/SDrecall/blob/main/doc/customRun.md#04-annotate-and-extract-regions-of-interest))
* A prioritized VCF (v4.2) (optional)

## Quick run
```{bash}
usage: SDrecall.py [-h] -i INPUT_BAM [-p PVCF] -r REF_GENOME -o OUTPUT -b BUILD -l LIST -a TABLE [-t THREAD] [--length LENGTH]
                   [--lower LOWER] [-f FRAGLEN] [-g GAPLEN] [--keep_trimmed] [-v VERBOSE]

SDrecall wrapper.

Options:
  -h, --help            show this help message and exit
  -i INPUT_BAM, --input_bam INPUT_BAM
                        Input BAM file
  -p PVCF, --pvcf PVCF  Path of prioritized VCF
  -r REF_GENOME, --ref_genome REF_GENOME
                        reference genome to extract SD regions
  -o OUTPUT, --output OUTPUT
                        output directory of resulting BED files
  -b BUILD, --build BUILD
                        reference genome assembly
  -l LIST, --list LIST  customized gene list
  -a TABLE, --table TABLE
                        gene annotation table
  -t THREAD, --thread THREAD
                        number of threads used with BISER (default = 8)
  --length LENGTH       BED window size for extracting reads from homologous regions (default: 250)
  --lower LOWER         lower bound for unlikely intrinsic variants
  -f FRAGLEN, --fraglen FRAGLEN
                        expected fragment length (default: 300)
  -g GAPLEN, --gaplen GAPLEN
                        small gap cutoff value (default: 10)
  --keep_trimmed        keep trimmed SD BED file for debugging
  -v VERBOSE, --verbose VERBOSE
                        verbosity level (default: INFO)

```
## Custom run
Users may run individual steps separately with script indexed 0-4. Although quick run is recommended for most use cases, custom run is helpful for debugging. Users can also allocate different number of threads for each step. Details of individual steps are described [here](doc/customRun.md).

## Workflow
<p align="center">
  <img src="doc/SDrecall.png" />
</p>

## Contact and correspondance
Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

## Simulation and benchmarking

To be implemented.

## Future improvements
- [ ] Change the source of annotation table
- [ ] Keep all homologous coordinates in a single file (avoid creating separate directories)
- [ ] Replace `samtools depth` with `mosdepth` or other faster algorithms
- [ ] Uplift dependencies
- [ ] Simulation and benchmarking


