## Custom Run
### 0. Prepare SD regions (BED format)
#### 0.1. Download reference files
##### Reference genome (hg19)
Current algorithm only supports hg19 build. Users can acquire UCSC reference genome (hg19) [here](https://github.com/creggian/ucsc-hg19-fasta).

##### Gene annotation file
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
Available options:
`-om`: combined exon & intron annotation table
`-oe`: exon annotation table
`-oi`: intron annotation table

Please specify the **full paths** (directory + file name) for all 3 arguments.

#### 0.2. Extract SD regions by BISER
```{bash}
./1_biserFetch.sh [-h] --ref-genome REF_GENOME --out OUTPUT_PATH [--thread THREAD=8]
```
Available options:
* `--ref-genome|-r`: path of reference genome (hg19)
* `--out|-o`: output path, file is written to OUTPUT_PATH/SD_hg19.bed
* `--thread|-t`: number of threads (Recommended: 6-10)

#### 0.3. Trim CIGAR strings of BISER output
```{bash}
./2_trimCIGAR.py [-h] -i INPUT_FN -o OUTPUT_FN [--fraglen | -f  FRAGLEN=300] [--gaplen | -g GAPLEN=10] [--verbose | -v VERBOSE=INFO]
```
Available options:
* `--input|-i`: input path of BED files
* `--output|-o`: output directory
* `--fraglen|-f`: fragment length cutoff (Default: 300)
* `--gaplen|-g`: small gap cutoff (Default: 10)
* `--verbose|-v`: verbosity level (any logging level stipulated in [python's logging module](https://docs.python.org/3/library/logging.html#logging-levels))

##### 0.3.1. Trimming logic
Each CIGAR string is interpreted as a list of tuples `(A, B)`, where `A` represents fragment length (FRAGLEN) and `B` describes the nature of the fragment (e.g. insertion, deletion, clipped sequence, etc.), abbreviated as (D, I, N, S, M). We extract 1) contiguous blocks of match/mismatch (M) with fragment length >= FRAGLEN, and 2) discontinuous blocks with total gap length <= GAPLEN. If no M block with fragment length >= FRAGLEN is found, we only consider the latter discontinuous blocks. In case of ambiguity, the blocks with the longest fragment length (sum of INT when LABEL == M) is reported. Results are saved to an BED file.

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
./3_annotateExtract.py [-h] -i INPUT -r REF -o OUTPATH [-l LIST] [-v VERBOSE=INFO]
```
Available options:
* `--input|-i`: input path of trimmed BED file (from previous step)
* `--out|-o`: output path
* `--ref|-r`: reference annotation table (tab-separated)
* `--list|-l`: list of genes of interest (genes not listed will be excluded) (optional)
* `--verbose|-v`: verbosity level

Users should devise a gene panel from a subset of `genelist.txt` [here](genelist.txt). The gene list should be a column of a table that contains:

| Genetic defect |
| -------------- |
| gene_1 |
| gene_2 |
| gene_3 |
| ... |

The following files are written to `OUTPATH/homologous_regions`:
* `all_homo_regions.bed`: BED containing the coordinates of all genes of interest
* `<gene>_related_homo_region.bed`: extracted BED with only \<gene\> data

The following files are written to `OUTPATH/principal_components`:
* `<gene>.bed`: start and end coordinates of \<gene\>

Format:

| chr | start pos | end pos | gene |
| --- | --------- | ------- | ---- |
| ... |  ... | ... | ... |
    
**WARNING: Users should never attempt to skip `--list` option. The script will run for a total of 5548 genes if unspecified which may take unexpectedly long time.** 

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
Available options:
* `--input-bam|-i`: input path of base quality score recalibrated (BQSR) BAM file
* `--ref-genome|-r`: path of reference genome
* `--region-list|-l`: file list of BED files / path of BED file
* `--ref-bed|-rb`: BED of reference genome (3 columns (in order): chromosome name, start position and length of chromosome (in bp)
* `--mq-threshold|-mq`: mapping quality (MQ) threshold (Default: 30)
* `--out|-o`: output path
* `--thread|-t`: number of threads (Default: 8) (For multiple regions analysis, total number of threads used = XARGS_THREADS $\times$ THREAD

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
* `--thread|-t`: number of threads (Default: 10)

Note: The script only handles ONE BED file each time. To analyse multiple regions, users can iterate over each BED file:
```{bash}
for region in $(find BED_DIR -name "*.bed" -type f | sort)
do
    ${ROOT}/5_maskedAlignPolyVarCall.sh --input-bam BAM_PATH --data DATAPATH --region-bed $region --ref-genome ${REF_GENOME} [--thread INT=10]
done
```
    
#### 2.1. Realign against masked genome
Paired-end FASTQs generated from BAM_PATH (stored in `fastq/`) are realigned to respective masked genomes stored in `masked_genome/`. Ploidy is estimated as follows. 

Depth_1 = Average depth of all reads in input BAM file
    
Depth_2 = Average depth of high-quality reads in input BAM file
    
Depth_3 = Average depth of all reads in the extracted BAM file
    
Depth_4 = Depth_3 - Depth_2
    
Ploidy = 2 $\times$ Depth_4 / Depth_1
    
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
