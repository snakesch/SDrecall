# SD analysis

:crystal_ball: TODO: Repo description :crystal_ball:
This pipeline extracts segmental duplication (SD) regions from a given genome, annotates the regions and extracts regions that intersect known causal genes of certain diseases. 

## Prerequisites
* [samtools](http://www.htslib.org/) >=v1.15.1
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) >=v2.27.1
* [BISER](https://github.com/0xTCG/biser) >=v1.1
* [seqkit](https://github.com/shenwei356/seqkit) >=v2.2.0
* [seqtk](https://github.com/lh3/seqtk) >=v1.3
* [GATK](https://gatk.broadinstitute.org/hc/en-us) >=v4.2.5
* [bwa](https://github.com/lh3/bwa) >=v0.7.17
* [GCC](https://gcc.gnu.org/) >=v9.1.0
* [Python](https://www.python.org/downloads/) >=v3.9
* [HTSlib](http://www.htslib.org/download/) >=v1.7
* [VCFPy](https://github.com/bihealth/vcfpy) >=v0.13.4
* [mosdepth](https://github.com/brentp/mosdepth)

## Input files
### Required
* Base Quality Score Recalibrated (BQSR) BAM file
* Gene annotation file (use NCBI RefSeq data if not specified)
### Optional
* A gene panel in BEDPE format (See [part 0.4](#### 0.4. Annotate and extract regions of interest))

## Usage
## Quick run
:crystal_ball: TODO: Quick run wrapper script :crystal_ball:

## Customized run
As users may allocate different number of threads in each step, customized run allows users to execute each step separately with a defined number of threads. Users are advised to use more threads for WGS data in some steps discussed below.

### 0. Prepare SD regions (BED format)
#### 0.1. Download reference files
##### 0.1.1. Reference genome (hg19)
Current algorithm only supports hg19 build. 

Users can acquire UCSC reference genome (hg19) [here](https://github.com/creggian/ucsc-hg19-fasta).

##### 0.1.2. Gene annotation file
Users are expected to provide a gene annotation file with the following format.
```
| Chrom | cdsStart | cdsEnd |   gene   | feature | strand | length |
| ----- | -------- | ------ |   ----   | ------- | ------ | ------ |
| chr1  |  11868   | 12227  | DDX11L17 | exon_1  |    +   |   359  |
```
Alternatively, follow the steps below to download from [NCBI RefSeq FTP server](ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt).
```{bash}
./geneAnnotation.py -om <merged output path> -oe <exon output path> -oi <intron output path>
```
Please specify the **full paths** (directory + file name) for all 3 arguments. `-om` writes the merged output (exon & intron); `-oe` writes only exon output; `-oi` writes only intron output. All 3 outputs will be formatted accordingly to the desired output.  

#### 0.2. Extract SD regions by BISER
```{bash}
./1_biserFetch.sh [-h] --ref-genome REF_GENOME --out OUTPUT_PATH [--thread THREAD=8]
```
This script writes extracted regions to <OUTPUT_PATH>/SD_hg19.bed (for hg19 build). By default, 8 threads are used.

Note: Users are advised to use 6-10 threads for this step. (<6 threads will lead to unnecessarily long execution time)

#### 0.3. Trim CIGAR strings of BISER output
```{bash}
./2_trimCIGAR.py [-h] -i INPUT_FN -o OUTPUT_FN [--fraglen | -f  FRAGLEN=300] [--gaplen | -g GAPLEN=10] [--verbose | -v VERBOSE=INFO]
```
Input BED file is trimmed by its CIGAR strings record by record. Resulting BED is written to OUTPUT_FN. By default, FRAGLEN = 300 and GAPLEN = 10. VERBOSE can be any logging level stipulated in [python's logging module](https://docs.python.org/3/library/logging.html#logging-levels).

##### 0.3.1. Trimming logic
For the CIGAR string of each record, we first format the string as a list of tuples (INT, LABEL) where INT is the fragment length and LABEL belongs to the set {D, I, N, S, M}. We aim to extract 1) contiguous block of match/mismatch (M) with fragment length >= FRAGLEN, and 2) discontinuous blocks with total gap length <= GAPLEN. In the case where no M fragment with fragment length >= FRAGLEN is found, we only consider the latter discontinuous blocks. In case overlapping blocks are
reported, the one with the longest fragment length (sum of INT when LABEL == M) is reported. For each block in a record, we record separately on the output BED file.

Example: (FRAGLEN = 300; GAPLEN = 10)
```{python3}
CIGAR = [ (2, M), (1, D), (35, M), (6, I), (7, M), (1, S), (2, I), (292, M), (5, D), (30, M) ]
extracted_block = [ (35, M), (6, I), (7, M), (1, S), (2, D), (292, M) ]
```

```{python3}
CIGAR = [ (2, M), (1, D), (5, M), (6, I), (7, M), (1, S), (2, D), (292, M), (5, D), (30, M) ]
extracted_block = [ (7, M), (1, S), (2, D), (292, M), (5, D), (30, M) ]
```

#### 0.4. Annotate and extract regions of interest 
```{bash}
./3_annotateExtract.py [-h] -i INPUT -r REF -o OUTPATH [-l LIST] [-c|--genecol GENECOL=16] [-v VERBOSE=INFO]
```
This step requires trimmed BED file from part 2 as INPUT. Gene annotation is done with reference to REF specified by `-r`. Users can also provide a list of genes of interest. Current alogorithm will extract specified genes for analysis (`-l`); it takes the whole annotated BED file for analysis if otherwise. GENECOL is a 0-based column index of "Genetic defect" in the given table/list (default = 16). 

The gene/region list should contain the following column:

| Genetic defect |
| -------------- |
| gene_1 |
| gene_2 |
| gene_3 |
| ... |

The following files are written to the path of INPUT:
* `*.homo.expanded.bed`: expanded two-way map of trimmed BED file
* `*.homo.expanded.geneanno.bed`: two-way map with annotations
* `*.homo.expanded.geneanno.region.bed`: two-way map of regions intersecting the provided LIST
The following files are written to OUTPATH:
* `all_homo_regions.bed`: BED containing the coordinates of all genes of interest (equivalent to the annotated BED output if no list is given)
* `<region>_related_homo_region.bed`: extracted BED with only <region> data

Format:

| chr | start pos | end pos | gene |
| --- | --------- | ------- | ---- |
| ... |  ... | ... | ... |

### 1. Preparation
We first extract fragments that align to extracted SD regions in the input BAM file, then build a masked genome against the region.

For single region analysis,
```{bash}
./4_preparation.sh [-h] --input-bam BAM --region-list REGION_BED --ref-genome REF_GENOME --ref-bed REF_BED --out OUTPATH [-mq INT=30] [--thread THREAD=8]
```
For multiple regions analysis,
```{bash}
find DIR -name "*.bed" | xargs -I{} -P XARGS_THREADS -t bash ./4_preparation.sh --input-bam BAM --region-list {} --ref-genome REF_GENOME --ref-bed REF_BED --out OUTPATH [-mq INT=30] [--thread THREAD=8]
```
BAM is the path of base quality score recalibrated (BQSR) BAM file, REF_GENOME is the indexed reference genome FASTA file, REF_BED contains 3 columns (in order): chromosome name, 0 (start pos) and length of the chromosome (in bp). `-mq` is the mapping quality filtering threshold. OUTPATH is the desired output directory. By default, THREAD = 8.

For multiple regions analysis, DIR is the directory containing all BED files created in step 0. Total number of threads used = `XARGS_THREADS x THREAD`.

In this step, we extract regions in REGION_BED from BAM and prepare masked genomes for each of the regions. A table of all files related to a certain region is also created. In the extraction step, we extract reads tagged with XA (multi-aligned) or with mapping quality MQ < _&alpha;_ where _&alpha;_ is the mapping quality threshold specified by `-mq`, which holds a default value of 30. Masked genomes are written to `OUTPATH/masked_genome/` and indexed.

**WARNING: This step may create extremely large FASTA files. Please ensure disk availability before running 4_preparation.sh.**

Note: More threads (>10) should be used for WGS data.

### 2. Masked alignment and polyploid variant calling
```{bash}
./5_maskedAlignPolyVarCall.sh [-h] --input-bam BAM_PATH --data DATAPATH --region-bed BED --ref-genome REF_GENOME [--thread INT=10]
```
Available options:
* `--input-bam`: input BQSR BAM file
* `--data`: parent directory of `masked_genome` and `fastq`
* `--region-bed`: BED file of selected SD regions

Note: The script only handles ONE BED file each time. To analyse multiple regions, users can iterate over each BED file:
```{bash}
for region in $(find BED_DIR -name "*.bed" -type f | sort)
do
    ${ROOT}/5_maskedAlignPolyVarCall.sh --input-bam BAM_PATH --data DATAPATH --region-bed $region --ref-genome ${REF_GENOME} [--thread INT=10]
done
```
    
#### 2.1. Realign against masked genome
Paired-end FASTQs generated from BAM_PATH (stored in `fastq/`) are realigned to respective masked genomes stored in `masked_genome/`. Ploidy is estimated as follows. 

$ Depth_1 = Average depth of all reads in input BAM file $
$ Depth_2 = Average depth of high-quality reads in input BAM file $
$ Depth_3 = Average depth of all reads in the extracted BAM file $
$ Depth_4 = Depth_3 - Depth_2 $
$ Ploidy = 2 \times {Depth_4}/{Depth_1} $
    
Variants are called if and only if ploidy >= 2. Regions are omitted if otherwise.
    
#### 2.2. Variant calling
Variants are first called with the assumption of possible multiploidy cases, genotyped, left-aligned and trimmed. If polyploidy is detected, it is converted to diploid VCF(s). Variants are first called with the assumption of possible multiploidy events. The called variants are then genotyped and cleaned to generate a VCF file. Variants called are filtered by coverage, only those with low coverage are retained (<=15).
    
### 3. Compare with intrinsic VCFs
```{bash}
./6_postProcessing.sh [-h] --vcfpath VCF_DIR --regions REGIONS_DIR
```
Available option:
* `--vcfpath`: directory of VCFs in poor coverage regions (the one generated in step 2)
* `--regions`: directory of all BED files

For variants called from poor coverage regions, they are compared to an "intrinsic VCF". Intrinsic VCFs contain variants ...

By default, `masked_genome/` and VCFs in `vcf/` are deleted to save space.

## Simulations

## Contact and correspondance
