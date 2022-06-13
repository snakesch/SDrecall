# SD analysis

This pipeline extracts segmental duplication (SD) regions from a given genome, annotates the regions and extracts regions that intersect known causal genes of certain diseases. 

## Prerequisites
* [samtools](http://www.htslib.org/)
* [BEDTools](https://bedtools.readthedocs.io/en/latest/)
* [BISER](https://github.com/0xTCG/biser) v1.1
* GCC
* Python3
* A gene panel (See part 3 below)

## Default run
### Required input files
* Base Quality Score Recalibrated (BQSR) BAM file
* Gene annotation file (use NCBI RefSeq data if not specified)
* A panel of genes in BEDPE format (optional)

:crystal_ball: TODO: Quick run wrapper script :crystal_ball:

## Customized run
## Usage
### 0. Download reference files
#### 0.1 Reference genome (hg19)
Current algorithm only supports hg19 build. Please skip gene annotation and gene filtering (part 3) for other builds. 

Users can acquire the FASTA file of hg19 build by UCSC [here](https://github.com/creggian/ucsc-hg19-fasta).

#### 0.2 Gene annotation file
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

### 1. Extract SD regions by BISER
```{bash}
./1_biserFetch.sh [-h] --ref-genome REF_GENOME --out OUTPUT_PATH [--thread THREAD]
```
This script writes extracted regions to <OUTPUT_PATH>/SD_hg19.bed (for hg19 build). By default, 8 threads are used.

Note: Users are advised to use 6-10 threads for this step. (<6 threads will lead to unnecessarily long execution time)

### 2. Trim CIGAR strings of BISER output
```{bash}
./2_trimCIGAR.py [-h] -i INPUT_FN -o OUTPUT_FN [--fraglen | -f  FRAGLEN] [--gaplen | -g GAPLEN] [--verbose | -v VERBOSE]
```
Input BED file is trimmed by its CIGAR strings record by record. Resulting BED is written to OUTPUT_FN. By default, FRAGLEN = 300 and GAPLEN = 10. VERBOSE can be any logging level stipulated in [python's logging module](https://docs.python.org/3/library/logging.html#logging-levels).

#### 2.1 Trimming logic
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

### 3. Annotate and extract regions of interest 
```{bash}
./3_annotateExtract.py [-h] -i INPUT -r REF [-l LIST] [-c|--genecol GENECOL] [-v VERBOSE]
```
This step requires trimmed BED file from part 2 as INPUT. Users should also provide a file for gene annotation (REF), and a panel of genes of interest. VERBOSE follows the usage documented in part 2.

#### 3.1 Gene annotation
Gene annotation file (hg19) is acquired as discussed in section 0.2. This path should be provided via the `--ref|-r` argument.

#### 3.2 Extract region(s) of interest (optional)
A gene panel (in a list) should be provided via `--list|-l` argument. If this argument is not given, the script terminates after completing gene annotation. Users should manually inspect the annotated BED file for the 0-based column index of "Genetic defect" and specify it by `--genecol|-c` (default: 16).

The gene/region list should contain the following column:
```{latex}
Genetic defect
gene_1
gene_2
gene_3
...
```
Outputs are as follows:
* `*.homo.expanded.bed`: expanded two-way map of trimmed BED file
* `*.homo.expanded.geneanno.bed`: two-way map with annotations
* `*.homo.expanded.geneanno.PID.bed`: two-way map of regions intersecting the provided LIST
* `*.homo.expanded.geneanno.PID.condensed.bed`: intersected two-way map with gene names concatenated if all other columns are the same

## Simulations
<TODO>

<------------- For personal reference only ------------>

This repository implements segmental duplication (SD) pipeline from Xingtian. This is supposed to be implemented mainly in python3.

## 0. Input (TODO)
* BAM file ready for variant calling / paired-end FASTQ files
* SD region data fetched from BISER
* (optional) a gene panel

## 1. Major features
- [x] Fetch SD data from BISER
- [x] Trim CIGAR strings of SD regions
- [x] Gene annotation and panel filter
- [ ] Deploy
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
