# SD analysis

This pipeline extracts segmental duplication (SD) regions from a given genome, annotates the regions and extracts regions that intersect known causal genes of certain diseases. 

## Prerequisites
* [samtools](http://www.htslib.org/)
* [BEDTools](https://bedtools.readthedocs.io/en/latest/)
* [BISER](https://github.com/0xTCG/biser) v1.1
* [seqkit](https://github.com/shenwei356/seqkit) v2.2.0
* [seqtk](https://github.com/lh3/seqtk) v1.3
* [GATK](https://gatk.broadinstitute.org/hc/en-us) v4.2.5
* [bwa](https://github.com/lh3/bwa) v0.7.17
* [GCC](https://gcc.gnu.org/) v9.1.0
* [Python3](https://www.python.org/downloads/)
* A gene panel (See part 3 below)

## Required input files
* Base Quality Score Recalibrated (BQSR) BAM file
* Gene annotation file (use NCBI RefSeq data if not specified)
* A panel of genes in BEDPE format (optional)

## Usage
## Quick run
:crystal_ball: TODO: Quick run wrapper script :crystal_ball:

## Customized run
### 0. Prepare SD regions (BED format)
#### 0.1. Download reference files
##### 0.1.1. Reference genome (hg19)
Current algorithm only supports hg19 build. 

Users can acquire the FASTA file of hg19 build by UCSC [here](https://github.com/creggian/ucsc-hg19-fasta).

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
./1_biserFetch.sh [-h] --ref-genome REF_GENOME --out OUTPUT_PATH [--thread THREAD]
```
This script writes extracted regions to <OUTPUT_PATH>/SD_hg19.bed (for hg19 build). By default, 8 threads are used.

Note: Users are advised to use 6-10 threads for this step. (<6 threads will lead to unnecessarily long execution time)

#### 0.3. Trim CIGAR strings of BISER output
```{bash}
./2_trimCIGAR.py [-h] -i INPUT_FN -o OUTPUT_FN [--fraglen | -f  FRAGLEN] [--gaplen | -g GAPLEN] [--verbose | -v VERBOSE]
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
./3_annotateExtract.py [-h] -i INPUT -r REF -o OUTPATH [-l LIST] [-c|--genecol GENECOL] [-v VERBOSE]
```
This step requires trimmed BED file from part 2 as INPUT. Gene annotation is done with reference to REF specified by `-r`. Users can also provide a list of genes of interest. Current alogorithm will extract specified genes for analysis (`-l`); it takes the whole annotated BED file for analysis if otherwise. GENECOL is a 0-based column index of "Genetic defect" in the given table/list (default = 16). 

The gene/region list should contain the following column:
```html
<table>
    <tr>
        <th>Genetic defect</th>
    </tr>
    <tr>
        <th>gene_1</th>
        <th>gene_2</th>
        <th>gene_3</th>
        <th><span>&#8320;
            </span>
            </th>
    </tr>
</table>
```
The following files are written to the path of INPUT:
* `*.homo.expanded.bed`: expanded two-way map of trimmed BED file
* `*.homo.expanded.geneanno.bed`: two-way map with annotations
* `*.homo.expanded.geneanno.region.bed`: two-way map of regions intersecting the provided LIST
The following files are written to OUTPATH:
* `all_homo_regions.bed`: BED containing the coordinates of all genes of interest (equivalent to the annotated BED output if no list is given)
* `<region>_related_homo_region.bed`: extracted BED with only <region> data
Format:
```{latex}
\begin{tabular}{|c|c|c|c|}
\hline
<chr> & <start pos> & <end pos> & <gene> \\
\vdots & \vdots & \vdots & \vdots \\
\hline
\end{tabular}
```
### 1. Preparation
We first extract fragments that align to SD regions from the input BAM file, build a masked genome, and then build a file name map for easy reference.

For single region analysis,
```{bash}
./4_buildFileMap.sh --input-bam BAM --region-list REGION_BED --ref-genome REF_GENOME --ref-bed REF_BED --out OUTPATH [-mq INT] [--thread INT]
```
For multiple regions analysis,
```
find DIR -name "*.bed" | xargs -I{} -t bash ./4_buildFileMap.sh --input-bam BAM --region-list {} --ref-genome REF_GENOME --ref-bed REF_BED --out OUTPATH [-mq INT] [--thread INT]
```
BAM is the path of base quality score recalibrated (BQSR) BAM file, REF_GENOME is the indexed reference genome FASTA file, REF_BED contains 3 columns (in order): chromosome name, 0 (start pos) and length of the chromosome (in bp). `-mq` is the mapping quality filtering threshold. OUTPATH is the desired output directory.

For multiple regions analysis, DIR is the directory containing all BED files created in step 0.

In this step, we extract regions in REGION_BED from BAM and prepare masked genomes for each of the regions. A table of all files related to a certain region is also created. In the extraction step, we extract reads tagged with XA (multi-aligned) or with mapping quality MQ < _&alpha;_ where _&alpha;_ is the mapping quality threshold specified by `-mq`, which holds a default value of 30. Masked genomes are written to `OUTPATH/masked_genome/` and indexed.

### 2. Polyploid variant calling in the priority regions

### 3. Compare with intrinsic VCFs

## Simulations

## Contact and correspondance

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
