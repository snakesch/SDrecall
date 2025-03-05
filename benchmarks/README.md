# SDrecall TEST BAM Files
## Download Information

The BAM files used for benchmarking SDrecall have been moved to Zenodo to reduce the repository size and improve distribution. These files contain aligned sequencing data needed for running the benchmarks.

**Zenodo Record:** [https://zenodo.org/records/14965091](https://zenodo.org/records/14965091)

## Available Files

- HG002.SD.deduped.hg19.bam - Human Genome sample HG002 aligned to hg19 reference
- HG002.SD.deduped.hg19.bam.bai - index to the above BAM file
- HG002.SD.deduped.hg38.bam - Human Genome sample HG002 aligned to hg38 reference
- HG002.SD.deduped.hg38.bam.bai - index to the above BAM file

## Download Instructions

You can download the BAM files directly using the following URLs:

```bash
# Download hg19 BAM file
wget https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg19.bam
# Download the hg19 BAM index
wget https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg19.bam.bai

# Download hg38 BAM file
wget https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg38.bam
# Download hg38 BAM index
wget https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg38.bam.bai
```

Or use the Zenodo API:

```bash
# Using curl
curl -O https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg19.bam
curl -O https://zenodo.org/records/14965091/files/HG002.SD.deduped.hg38.bam
```

## Generated from the BAM file of sample HG002 offered by Genome In A Bottle 

Original BAM file on GRCh37 was downloaded from `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam`

Original BAM file on GRCh38 was downloaded from `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam`

Reference Segmental Duplication pairs on hg19 was downloaded from `https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_doSchema=describe+table+schema`

Reference Segmental Duplication pairs on hg38 was downloaded from `https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_doSchema=describe+table+schema`

Showcase of shell commands used to generate the test BAM sample from the original BAM file with average depth of 300X.

```bash
# First Reference SD pairs are expanded to make all SD regions list in the first 3 columns
sambamba slice -L ${REF_SDs} ${RAW_BAM} | \
samtools view --subsample-seed 0 --subsample 0.1 | cut -f 1 | sort - | uniq - > ${QNAME_LIST} && \
samtools view -h -u -N ${QNAME_LIST} ${RAW_BAM} | \
samtools fastq -1 ${FORWARD_READs} -2 ${REVERSE_READs} -0 /dev/null -n && \
bwa mem -M -R "@RG\tID:HG002\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${HG002}" \
ucsc.hg19.fasta ${FORWARD_READs} ${REVERSE_READs} | \
samtools view -uSh - | samtools sort -O bam -o HG002.SD.bam && \
gatk \
--java-options "-Xmx28G" \
MarkDuplicates \
-I HG002.SD.bam \
-O HG002.SD.deduped.bam \
-M HG002.SD.deduped.duplicateindex \
-ASO coordinate \
--CREATE_INDEX \
--REMOVE_DUPLICATES \
--VALIDATION_STRINGENCY LENIENT \
-CO "GATK_MARKDUP_DONE" && \
samtools index HG002.SD.deduped.bam
```

## Citation

If you use these datasets in your research, please cite both the SDrecall repository and the Zenodo record:

```
@software{sdrecall_2023,
  author = {Yang, Xing Tian and [Other Authors]},
  title = {SDrecall: A Scalable Approach for Sensitive Variant Detection in Segmental Duplications},
  url = {https://github.com/snakesch/SDrecall},
  year = {2023}
}

@dataset{yang_2023_14965091,
  author = {Yang, Xing Tian and [Other Authors]},
  title = {SDrecall Benchmark BAM Files},
  year = {2023},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.14965091},
  url = {https://zenodo.org/records/14965091}
}
```


# SDrecall TEST target region
To allow users to quickly test the integrity and validity of the SDrecall workflow, we offer a small target region targeting the coding segments of around 500 genes known to be causal to Primary ImmunoDeficiencies

Files are offered in this github repository:
- benchmarks/hg19.coding.pid.pad20.bed
- benchmarks/hg38.coding.pid.pad20.bed

