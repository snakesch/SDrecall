#!/usr/bin/env bash
set -euo pipefail

base="${BASE:-/paedyl01/disk1/yangyxt/test_tmp/sex_haplotag_gtploidy_20260706}"
runner="${RUNNER:-/paedyl01/disk1/yangyxt/SDrecall-rust-migration/scripts/run_sex_gt_ploidy_haplotag.sh}"
queue="${QUEUE:-large}"
ncpus="${NCPUS:-25}"
mem="${MEM:-160gb}"
walltime="${WALLTIME:-04:00:00}"

mkdir -p "${base}/pbs_logs"

submit() {
    local assembly="$1"
    local ploidy_class="$2"
    local vcf="$3"
    local raw_bam="$4"
    local reference="$5"
    local out_dir="${base}/${assembly}/sex_${ploidy_class}"
    local parquet_out="${base}/parquet_by_ploidy/${assembly}.sex_${ploidy_class}.haplotags.parquet"
    local log="${base}/pbs_logs/${assembly}.sex_${ploidy_class}.log"
    local job_name="sexgt_${assembly}_${ploidy_class}"

    mkdir -p "${out_dir}" "$(dirname "${parquet_out}")" "$(dirname "${log}")"

    qsub -V \
        -N "${job_name}" \
        -q "${queue}" \
        -l "select=1:ncpus=${ncpus}:mem=${mem},walltime=${walltime}" \
        -- /bin/bash -lc \
        "source ~/.bashrc && conda activate ngs_pipeline && bash '${runner}' --assembly '${assembly}' --ploidy-class '${ploidy_class}' --vcf '${vcf}' --raw-bam '${raw_bam}' --reference '${reference}' --out-dir '${out_dir}' --parquet-out-dir '${parquet_out}' --sample HG002 --threads '${ncpus}' --skip-existing > '${log}' 2>&1"
}

submit_hg19() {
    submit \
        hg19 "$1" \
        /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/defrabbV0.020_phased/GRCh37_HG2-T2TQ100-V1.1_smvar.vcf.gz \
        /paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/HG002.hs37d5.300x.bam \
        /paedyl01/disk1/yangyxt/indexed_genome/GRCh37/human_g1k_v37_decoy.fasta
}

submit_hg38() {
    submit \
        hg38 "$1" \
        /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/defrabbV0.020_phased/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz \
        /paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/HG002.GRCh38.300x.bam \
        /paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta
}

submit_hg19 haploid
submit_hg19 diploid
submit_hg38 haploid
submit_hg38 diploid
