## Custom Run

Note: All commands described in this documentation should be executed from the top directory of git repository. (eg. ~/SDrecall)

### 0. Reference files
#### Reference genome

Depending on different genomic assembly, users have to provide a reference genome file (fasta). For testing, users may download hg19 reference genome [here](https://github.com/creggian/ucsc-hg19-fasta).

#### Annotation table

Gene annotation data are downloaded from [NCBI RefSeq FTP server](ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt) with `0_geneAnnotation.py`.

```{bash}
usage: 0_geneAnnotation.py [-h] -om OUTPUT_MERGED -oe OUTPUT_EXON -oi OUTPUT_INTRON [-v VERBOSE]
                           [-tg TARGET_GENES]

optional arguments:
  -h, --help            show this help message and exit
  -om OUTPUT_MERGED, --output_merged OUTPUT_MERGED
                        Output path of merged annotation table (exon & intron)
  -oe OUTPUT_EXON, --output_exon OUTPUT_EXON
                        Output path of exon annotation table
  -oi OUTPUT_INTRON, --output_intron OUTPUT_INTRON
                        Output path of intron annotation table
  -v VERBOSE, --verbose VERBOSE
                        Verbosity level
  -tg TARGET_GENES, --target_genes TARGET_GENES
                        Target gene list (one gene per row)
```

### 1. Get relevant SD regions

Segmental duplication (SD) regions are first ascertained with BISER, then trimmed by their CIGAR strings (see Trimming logic below). The trimmed SD regions are annotated with the annotation table described in section 0. Homologous region (all SD pairs) and principal component (SD within gene) coordinates are extracted to the specified output directory in BED format (in `out/homologous_regions` and `out/principal_components`). Only genes in the given `genelist.txt` will be included.

```{bash}
usage: 1_getTrimmedSD.py [-h] -r REF_GENOME -o OUTPUT -b BUILD -l LIST -a TABLE [-t THREAD]
                         [-f FRAGLEN] [-g GAPLEN] [--keep_trimmed] [-v VERBOSE]

Extract trimmed SD regions from reference genome.

Options:
  -h, --help            show this help message and exit
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
  -f FRAGLEN, --fraglen FRAGLEN
                        expected fragment length (default: 300)
  -g GAPLEN, --gaplen GAPLEN
                        small gap cutoff value (default: 10)
  --keep_trimmed        keep trimmed SD BED file for debugging
  -v VERBOSE, --verbose VERBOSE
                        verbosity level (default: INFO)
```

##### 1.1. Trimming logic
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

### 2. Preparation

Multi-aligned reads with mapping quality lower than the specified threshold (default: 30; specified via `-mq`) mapped to gene-related homologous regions (BED files in the directory `ref/homologous_regions`) are extracted from the input BAM file for each gene and output in BAM and FASTQ format. Reference genome is masked against homologous regions within a gene (BED files in the directory `ref/principal_components`) and written to `out/masked_genome` by default. Intrinsic variants are also called from the input BAM file for downstream comparison.

```{bash}
usage: 2_preparation.py [-h] --input_bam INPUT_BAM --ref_genome REF_GENOME --outpath OUTPATH [--homo_dir HOMO_DIR] [--pc_dir PC_DIR] 
                        [-mq MAPPING_QUALITY] [--length LENGTH] [--thread THREAD] [-v VERBOSE]

Preparation for masked alignment.

Options:
  -h, --help            show this help message and exit
  --input_bam INPUT_BAM
                        input BAM file
  --homo_dir HOMO_DIR   directory of BED files of homologous regions
  --pc_dir PC_DIR       directory of BED files of principal components
  -mq MAPPING_QUALITY, --mapping_quality MAPPING_QUALITY
                        MQ threshold
  --ref_genome REF_GENOME
                        reference genome
  --outpath OUTPATH     output directory for FASTQ and masked genomes
  --length LENGTH       BED window size for extracting reads from homologous regions (default: 250)
  --thread THREAD       number of threads (default: 8)
  -v VERBOSE, --verbose VERBOSE
                        verbosity level (default: INFO)

```

#### 2.1. Intrinsic VCF generation

Homologous regions are first broken down into smaller regions with size <= LEN. For each block with LEN >= 100, its sequence is extracted from FASTA of reference genome and then converted to FASTQ format. A masked genome is prepared against all principal components (i.e. gene regions, including 5'-UTR and 3'-UTR), that is masking regions other than principal components. FASTQ files are realigned to the masked genome and intrinsic variants are called.

### 3. Masked (re-)alignment and multiploid variant calling

Input BAM is first re-aligned against respective masked genomes. Ploidy status is deduced by average read depths (see Ploidy estimation). Variants are called from the realigned BAM files with GATK, corrected for their ploidy status such that the resulting VCF is diploid (see Ploidy correction), and filtered by coverage such that low coverage variants (<=15) are retained. The resulting VCF is then compared against its intrinsic counterpart generated from the same input BAM file (see intrinsic variant labelling).  

```{bash}
usage: 3_maskedAlignPolyVar.py [-h] --input_bam INPUT_BAM [--bed_dir BED_DIR] --intrinsic_vcf INTRINSIC_VCF
                               [--lower LOWER] [--masked_genomes MASKED_GENOMES] [--fastq_dir FASTQ_DIR]
                               [--ref_genome REF_GENOME] [--output_vcf OUTPUT_VCF] [--thread THREAD] [--keep_vcf]
                               [-v VERBOSE]

Masked alignment and multi-ploidy variant calling.

Options:
  -h, --help            show this help message and exit
  --input_bam INPUT_BAM
                        input BAM file
  --bed_dir BED_DIR     directory of BED files
  --intrinsic_vcf INTRINSIC_VCF
                        path of intrinsic VCF
  --lower LOWER         lower bound for unlikely intrinsic variants
  --masked_genomes MASKED_GENOMES
                        directory of masked genomes
  --fastq_dir FASTQ_DIR
                        directory of FASTQ files
  --ref_genome REF_GENOME
                        path of reference genome
  --output_vcf OUTPUT_VCF
                        directory to write VCF
  --thread THREAD       number of threads
  --keep_vcf            keep intermediate VCFs
  -v VERBOSE, --verbose VERBOSE
                        verbosity level (default: INFO)
```

#### 3.1. Ploidy estimation

```
Depth_1 $=$ Average depth of all reads in input BAM file
    
Depth_2 $=$ Average depth of high-quality reads in input BAM file
    
Depth_3 $=$ Average depth of all reads in the extracted BAM file
    
Depth_4 $=$ Depth_3 - Depth_2
    
Ploidy $=$ 2 $\times$ Depth_4 / Depth_1
```

Note: Variants are called if and only if ploidy >= 2. Regions are omitted if otherwise.

#### 3.2. Ploidy correction 

For each variant with ploidy > 2, if the called genotype:
* has more than 2 ALT alleles: changed to `1/1`
* has only 1 ALT allele: changed to `0/1`
* is missing and has no ALT allele: changed to `.`
* not missing and has no ALT allele: changed to `0/0`.

For each haploid variant, if the called genotype is:
* `1`: changed to `1/1`
* `0`: changed to `0/0`
* `.`: no change

#### 3.3. Intrinsic variant labelling 

For each variant, we compare the query VCF (Query) to the intrinsic VCF (Intrinsic) and determine from VCFs the allelic depth (AD) of reference allele (REF) and alternate allele (ALT). From the allelic depths, we compute two ratios, `qra_ratio` and `ira_ratio`.

```
qra_ratio $=$ $\frac{Query REF AD}{Query ALT AD}$

ira_ratio $=$ $\frac{Intrinsic REF AD + 1}{Intrinsic ALT AD}$

```

If `ira_ratio` equals 0, the variant is considered unlikely intrinsic. If `$\frac{qra_ratio}{ira-ratio} <= $` lower limit specified via `--lower`, the variant is regarded as likely intrinsic. All other variants are labelled as unlikely intrinsic.

### 4. Merge with prioritized VCF (optional)

Users may compare the resulting VCF from step 3 (original VCF) with another VCF (prioritized VCF). Variants found in original VCF will be tagged with `ov_tag` in the FITLER field whereas variants found in prioritized VCF will be tagged with `pv_tag`. Variants found in both VCFs will be tagged twice.

```{bash}
usage: 4_mergePVCF.py [-h] --pvcf PVCF --ovcf OVCF --pv_tag PV_TAG --ov_tag OV_TAG --outpath OUTPATH
                      [--thread THREAD] [-v VERBOSE]

Merge prioritized VCF and original VCF.

Options:
  -h, --help            show this help message and exit
  --pvcf PVCF           prioritized VCF (gz)
  --ovcf OVCF           original VCF (gz)
  --pv_tag PV_TAG       tag used for variants from prioritized VCF
  --ov_tag OV_TAG       tag used for variants from original VCF
  --outpath OUTPATH     absolute output path of merged VCF (gz)
  --thread THREAD       number of threads (default: 8)
  -v VERBOSE, --verbose VERBOSE
                        verbosity level (default: INFO)

```


