#!/bin/bash

# This script contains miscellaneous helper functions for proper presentation and error handling.
# Author: Louis She (2022-04)
# Contact: louisshe@hku.hk

# Helper functions
function timestamp () {
    echo -e "[$(date +%a) $(date +%b-'%-m') $(date +'%-I':%M:%S%P)]"
}

function index_genome () {
    # Possible arguments: bwtsw, fai, dict
    set -e
    local fasta="$1"
    shift
    while [[ $# -gt 0 ]]
    do
        local idx="$1"
        shift

        if [[ "${idx}" == "bwtsw" ]]
        then
            if  [[ ! -f ${fasta}.pac ]] || [[ ${fasta} -nt ${fasta}.pac ]] || \
                [[ ! -f ${fasta}.bwt ]] || [[ ${fasta} -nt ${fasta}.bwt ]] || \
                [[ ! -f ${fasta}.ann ]] || [[ ${fasta} -nt ${fasta}.ann ]] || \
                [[ ! -f ${fasta}.amb ]] || [[ ${fasta} -nt ${fasta}.amb ]] || \
                [[ ! -f ${fasta}.sa ]] || [[ ${fasta} -nt ${fasta}.sa ]]; then
                bwa index -a bwtsw -b 380000000 "${fasta}"
            fi
        fi

        if [[ "${idx}" == "fai" ]]
        then
            if [[ ! -f "${fasta}.fai" ]] || [[ "${fasta}" -nt "${fasta}.fai" ]]
            then
                samtools faidx "${fasta}"
            fi
        fi

        if [[ "${idx}" == "dict" ]]
        then
            gatk CreateSequenceDictionary -R "${fasta}"
        fi

        echo -e "$(timestamp) INFO: Successfully created index file for ${fasta}"
    done
    return 0
}

function validate_mgenome {
    local mgenome="$1"
    local target_region="$2"

    [[ -f "${mgenome}" ]] || return 3;
    [[ -f "${target_region}" ]] || return 3;

    # main checking
    N_RECORDS=$(bedtools getfasta -fi "${mgenome}" -bed "${target_region}" | seqkit locate -p N -V 0 -G -I | { grep -c N$ || true; })
    if [[ "${N_RECORDS}" -ne 0 ]]
    then
        echo -e >&2 "$(timestamp) ERROR: Masked genome ${mgenome} is invalid. Pls delete it and rerun."
        echo -e >&2 "$(timestamp) ERROR: ${N_RECORDS} N bases detected. (Expected = 0)"
        return -1
    else
        return 0
    fi
}

function getDepth () {
    local bamf=$1
    local target_region=$2
    local threads=$3
    local readType=$4 # HQ
    if [[ -f $bamf ]] && [[ -f $target_region ]] && [[ $(cat $target_region | wc -l) -gt 0 ]]
    then
        if [[ $readType == "HQ" ]]
        then
            samtools view -h -@ $threads $bamf | grep -v XA:Z | samtools depth -@ $threads -a -b $target_region -Q 30 -g DUP,UNMAP,QCFAIL,SECONDARY - | awk '{sum+=$3} END{ printf "%s",sum/NR;}'
        else
            samtools depth -@ $threads -a -b $target_region -g DUP,UNMAP,QCFAIL,SECONDARY $bamf | awk '{sum+=$3} END{ printf "%s",sum/NR;}'
        fi
    else
        return 127
    fi
}

function isValidVCF () {
    while [[ $# -gt 0 ]]
    do
        local vcf="$1"
        shift

        if [[ ! -f $vcf ]]; then return 1; fi

        if [[ $vcf == *.gz ]]; then
            zgrep -m 1 -v "^#" $vcf || { echo -e >&2 "$(timestamp) ERROR: No records found in $vcf" && return 1; }
        else
            grep -m 1 -v "^#" $vcf || { echo -e >&2 "$(timestamp) ERROR: No records found in $vcf" && return 1; }
            cat $vcf | bgzip -c > ${vcf}.gz
        fi

        if [ ! -f ${vcf}.tbi ]; then tabix -p vcf $vcf; fi
        gatk ValidateVariants -V $vcf --validation-type-to-exclude ALL 1>/dev/null 2>&1
    done
}

function pick_poorcov_region_bam {
    local OPTIND i t o d
    while getopts i:t::o::d:: args
    do
        case ${args} in
            i) local input_bam=$OPTARG ;;
            t) local target_region=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            d) local depth_threshold=$OPTARG;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
    
    samtools view -hb -@ ${NTHREADS} -P -L "${target_region}" "${input_bam}" > "${output_bed/.txt}.bam.tmp"
    if ! $(dirname $0)/checkBAM.sh "${output_bed/.txt}.bam.tmp"; then exit 3; fi
    mosdepth -Q 10 -b "${target_region}" "${input_bam/.bam/.highMQ}" "${input_bam}"
    [ -f "${input_bam/.bam/.highMQ}.per-base.bed.gz" ] || { echo -e >&2 "$(timestamp) ERROR: mosdepth out of memory." && exit 2; }
    
    # mosdepth should create .bed.gz here, exit otherwise
    zcat "${input_bam/.bam/.highMQ}.per-base.bed.gz" | awk -F '\t' '$4 <= '${depth_threshold}' {print;}' > ${output_bed/.gz/}
    rm -f "${output_bed/.txt}.bam.tmp"
}

# Alias
alias time="/usr/bin/time -f '$(timestamp) INFO: Elapsed time: %es\tMemory usage: %M KB\tCPU usage: %P'"
