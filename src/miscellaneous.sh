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

        if [ ! -f ${vcf}.csi ]; then bcftools index $vcf; fi
        gatk ValidateVariants -V $vcf --validation-type-to-exclude ALL 1>/dev/null 2>&1
    done
}

# - Empirical function - #
function randomID {
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}

function quick_get_bam_depth_MQ () {
    local OPTIND i t o q
    while getopts i:t::o::q:: args
    do
        case ${args} in
            i) local input_bam=$OPTARG ;;
            t) local target_region=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            q) local MQ_threshold=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
    if [[ -z ${MQ_threshold} ]]; then local MQ_threshold=10; fi

    if [[ ! -z ${target_region} ]]; then
        local target_region_ID=$(basename ${target_region} | awk -F '.' '{printf "%s", $1;}')
        local targeted_bam=${input_bam/.bam/}.${target_region_ID}.bam
        bash $(dirname $0)/select_bam_by_regions.sh \
        -b ${target_region} \
        -i ${input_bam} \
        -o ${targeted_bam} && \
        local input_bam=${targeted_bam}
    fi

    # Function used to output a bed file
    check_bam_validity ${input_bam}

    # Set MQ threshold at 10, because that's when reads are abandoned  by HaplotypeCaller in GATK
    if [[ ! -z ${target_region} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Command line is: ${mosdepth} -Q ${MQ_threshold} -b ${target_region} ${input_bam/.bam/.highMQ} ${input_bam}"
        mosdepth -Q ${MQ_threshold} -b ${target_region} ${input_bam/.bam/.highMQ} ${input_bam} && \
        ls -lh ${input_bam/.bam/.highMQ}.per-base.bed.gz || \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": ${mosdepth} still out of memory."
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Command line is: ${mosdepth} -Q ${MQ_threshold} ${input_bam/.bam/.highMQ} ${input_bam}"
        mosdepth -Q ${MQ_threshold} ${input_bam/.bam/.highMQ} ${input_bam} && \
        ls -lh ${input_bam/.bam/.highMQ}.per-base.bed.gz || \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": ${mosdepth} still out of memory."
    fi

    if [[ ! ${output_bed} =~ \.bed(\.gz)*$ ]]; then
        local ob_var_name=${output_bed}
        local output_bed=${input_bam/.bam/.highMQ}.per-base.bed.gz
    fi

    if [[ ! -z ${ob_var_name} ]]; then
        eval ${ob_var_name}="'${output_bed}'"
    fi
}

function pick_poorcov_region_bam {
    local OPTIND i t o q d
    while getopts i:t::o::q::d:: args
    do
        case ${args} in
            i) local input_bam=$OPTARG ;;
            t) local target_region=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            q) local MQ_threshold=$OPTARG ;;
            d) local depth_threshold=$OPTARG;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
    local tmp_bed=$(dirname $output_bed)/$(randomID).bed
    touch ${tmp_bed} && chmod +w ${tmp_bed}
    if [[ -z ${MQ_threshold} ]]; then local MQ_threshold=10; fi
    if [[ -z ${depth_threshold} ]]; then local depth_threshold=10; fi

    if [[ -z ${target_region} ]]; then
        quick_get_bam_depth_MQ -i ${input_bam} -o ${tmp_bed} -q ${MQ_threshold}
    else
        quick_get_bam_depth_MQ -i ${input_bam} -t ${target_region} -o ${tmp_bed} -q ${MQ_threshold}
    fi

    pick_poorcov_region_bed -i ${tmp_bed} -o ${output_bed} -d ${depth_threshold}
    rm -f "${tmp_bed}"
}

function check_bam_validity {
    local input=${1}
    local expected_lines=${2}
    if [[ -z ${expected_lines} ]]; then local expected_lines=10; fi

    if [[ ! -f ${input} ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, found ${input} not even exist."
        return 1
    fi

    local line_num=$( samtools view -c ${input} 2>&1 | awk '{print;}' - )

    if [[ ${line_num} -lt ${expected_lines} ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, found ${input} only has ${line_num} lines where we expect it to have ${expected_lines} lines."
        return 1
    fi
    if check_bam_qname_format ${input}; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${input} has proper query name format."
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${input} has inappropriate query name format."
        return 1
    fi

    check_bam_index ${input}
    local tmp_check_log=${input::-4}.$(randomID).check

    >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} format check log is ${tmp_check_log}"
    ${gatk} ValidateSamFile \
    -I ${input} \
    --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE \
    -MO 1000000 > ${tmp_check_log} 2>&1

    local -a err_types=($(awk '$1 ~ /^ERROR::/{print $1;}' ${tmp_check_log} | sort - | uniq - | awk -F ':' '{printf "%s ",$3;}'))
    if [[ ${#err_types[@]} -eq 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): No format errors detected in ${input_bam}"
    else
        for err in "${err_types[@]}"; do
            if [[ ! ${err} =~ "INVALID_INDEX_FILE_POINTER" ]] && [[ ! ${err} =~ "MATE_NOT_FOUND" ]] && [[ ! ${err} =~ "INVALID_PLATFORM_VALUE" ]] && [[ ! ${err} =~ "MISSING_READ_GROUP" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} is corrupted and cannot be simply fixed."
                return 1
            elif [[ ${err} =~ "INVALID_INDEX_FILE_POINTER" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} has corrupted index file, re-indexing it."
                #First delete the old index file
                rm -f ${input::-1}i || >&2 echo "In function ${FUNCNAME}: $(timestamp): Index file ${input::-1}i already deleted."
                rm -f ${input}.bai || >&2 echo "In function ${FUNCNAME}: $(timestamp): Index file ${input}.bai already deleted."
                samtools index ${input}
            elif [[ ${err} =~ "MATE_NOT_FOUND" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} has some singleton reads. Just leave them there."
            elif [[ ${err} =~ "INVALID_PLATFORM_VALUE" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} might be a simulated golden bam file, which does not have a valid PL value. Let it pass. The error message is:"
                >&2 awk '$1 ~ /^ERROR::INVALID_PLATFORM_VALUE/{print $0;}' ${tmp_check_log}
            elif [[ ${err} =~ "MISSING_READ_GROUP" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} does not have the read group line. Just ignore."
            fi
        done
        # IF the for loop above is executed without return 1, then no other error types are identified. return
    fi
    rm ${tmp_check_log} 2> /dev/null || >&2 echo "${tmp_check_log} already deleted"
}

function check_bam_qname_format {
    local input_bam=${1}
    local threads=${2}

    if [[ -z ${threads} ]]; then local threads=1; fi

    local issue_lines=$(samtools view ${input_bam} 2>&1 | head -20 | awk -F '\t' '$1 ~ /\/[1-2]$/{print;}' | wc -l)
    if [[ ${issue_lines} -eq 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): No malformed lines found."
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: The bam file ${input_bam} contains problematic query name format: "
        { samtools view ${input_bam} | head -10
        samtools view -h ${input_bam} | \
        awk -F '\t' '{gsub(/\/[1-2]$/, "", $1); print;}' | \
        samtools sort -O bam -o ${input_bam/.bam/.temp.bam} -@ ${threads} - ;} && \
        samtools index ${input_bam/.bam/.temp.bam} && \
        mv ${input_bam/.bam/.temp.bam} ${input_bam} && \
        mv ${input_bam/.bam/.temp.bam.bai} ${input_bam}.bai && \
        { rm -f ${input_bam/.bam/.bai} || >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: ${input_bam/.bam/.bai} not existed, does not have to delete it."; } && \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Manually remove the hanging /1 or /2 tag at the end of query names. New bam is $(ls -lh ${input_bam})" && \
        return || \
        { >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Failed to manually modify the query name format on ${input_bam}." && return 1; }
    fi
}

function check_bam_index {
    local bam=${1}
    if [[ ! -f ${bam}.bai ]] && [[ ! -f ${bam::-1}i ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Index file for ${bam} not existed."
        samtools index ${bam} || \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Indexing ${bam} failed. Need to sort it first."
    elif [[ -f ${bam}.bai ]] && [[ ${bam} -nt ${bam}.bai ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Index file for ${bam} is older than ${bam} itself. Reindexing"
        samtools index ${bam} && ls -lh ${bam}* || \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Indexing ${bam} failed. Need to sort it first."
    elif [[ -f ${bam::-4}.bai ]] && [[ ${bam} -nt ${bam::-4}.bai ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Index file for ${bam} is older than ${bam} itself. Reindexing"
        samtools index ${bam} || \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Indexing ${bam} failed. Need to sort it first."
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Index file for ${bam} is valid."
    fi
}

function pick_poorcov_region_bed {
    local OPTIND i t o d
    while getopts i:t::o::d:: args
    do
        case ${args} in
            i) local input_bed=$OPTARG ;;
            t) local target_region=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            d) local depth_threshold=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
    if [[ -z ${depth_threshold} ]]; then local depth_threshold=10; fi
    # Function used to output a bed file with the 4th column recording depth information

    if [[ -z ${target_region} ]]; then
        if [[ ! ${output_bed} =~ \.bed(\.gz)*$ ]]; then
            local ob_var_name=${output_bed}
            local output_bed=${input_bed/.bed*/.lowDepth.bed}
        fi
        if [[ ${input_bed} =~ \.gz$ ]]; then
            zcat ${input_bed} | \
            awk -F '\t' '$4 <= '${depth_threshold}' {print;}' > ${output_bed/.gz/}
        else
            awk -F '\t' '$4 <= '${depth_threshold}' {print;}' ${input_bed} > ${output_bed/.gz/}
        fi
    else
        local target_region_ID=$(basename ${target_region} | awk -F '.' '{printf "%s", $1;}')
        if [[ ! ${output_bed} =~ \.bed(\.gz)*$ ]]; then
            local ob_var_name=${output_bed}
            local output_bed=${input_bed/.bed*/.lowDepth}.${target_region_ID}.bed
        fi
        bedtools intersect -wa -a ${input_bed} -b ${target_region} | \
        awk -F '\t' '$4 <= '${depth_threshold}'{print;}' > ${output_bed/.gz/}
    fi

    if [[ ${output_bed} =~ \.gz$ ]]; then
        gzip -f -c ${output_bed/.gz/} > ${output_bed}
    fi

    if [[ ! -z ${ob_var_name} ]]; then
        eval ${ob_var_name}="'${output_bed}'"
    fi
}

# function check_vcf_validity {
#     local input_vcf=${1}
#     local expected_lines=${2}
#     local expected_samples=${3}

#     if [ -z ${expected_lines} ]; then expected_lines=1; fi
#     if [ ! -f ${input_vcf} ]; then
#         >&2 echo "$(timestamp): In function ${FUNCNAME}, ${input_vcf} does not even exist."
#         return 1
#     fi

#     if [[ ${input_vcf} =~ \.gz$ ]]; then
#         >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} should be bgzipped, check gz vcf validity"
#         if check_gz_vcf ${input_vcf} ${expected_lines}; then
#             if check_vcf_samples ${input_vcf} ${expected_samples}; then
#                 >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${input_vcf} has solid sample names as expected."
#             else
#                 return 1
#             fi
#         elif check_plain_vcf ${input_vcf::-3} ${expected_lines} && [[ ${input_vcf::-3} -nt ${input_vcf} ]]; then
#             if check_vcf_samples ${input_vcf} ${expected_samples}; then
#                 bgzip -f -c ${input_vcf::-3} > ${input_vcf} && \
#                 tabix -f -p vcf ${input_vcf} && return || \
#                 return 1
#             else
#                 return 1
#             fi
#         else
#             return 1
#         fi
#     elif [[ ${input_vcf} =~ \.vcf$ ]]; then
#         >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} should be plain text, check plain vcf validity"
#         if check_plain_vcf ${input_vcf} ${expected_lines}; then
#             if check_vcf_samples ${input_vcf} ${expected_samples}; then
#                 >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${input_vcf} has solid sample names as expected."
#             else
#                 return 1
#             fi
#         elif check_gz_vcf ${input_vcf}.gz ${expected_lines} && [[ ${input_vcf}.gz -nt ${input_vcf} ]]; then
#             if check_vcf_samples ${input_vcf} ${expected_samples}; then
#                 gunzip -f -c ${input_vcf}.gz > ${input_vcf} && return || \
#                 return 1
#             else
#                 return 1
#             fi
#         else
#             return 1
#         fi
#     else
#         return 2
#     fi
# }

# function check_gz_vcf {
#     local input_vcf=${1}
#     local expected_lines=${2}
#     if [ -z ${expected_lines} ]; then expected_lines=1; fi
#     if [ ! -f ${input_vcf} ]; then return 1; fi

#     if check_vcf_format ${input_vcf}; then  # Check input_vcf format
#         >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has solid vcf format. Check whether contains enough record number"
#         if [ $(zcat ${input_vcf} | awk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then  # Check input_vcf content
#             >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has enough records."
#         else
#             >&2 echo "${input_vcf} has $(zcat ${input_vcf} | awk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
#             return 3
#         fi
#     else
#         return 2
#     fi
# }

# function check_vcf_samples {
#     local input_vcf=${1}
#     local expected_samples=${2}
#     # expected samples should be delimited by comma

#     if [[ ${input_vcf} =~ \.vcf$ ]]; then
#         bgzip -f -c ${input_vcf} > ${input_vcf/.vcf/.tmp.vcf.gz} && \
#         tabix -f -p vcf ${input_vcf/.vcf/.tmp.vcf.gz} && \
#         local gz_vcf=${input_vcf/.vcf/.tmp.vcf.gz} && \
#         local asks=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}') && \
#         rm -f ${gz_vcf} || { >&2 echo "Failed to compress ${input_vcf}" && return 1; }
#         if [[ $? != "0" ]]; then return 1; fi
#     else
#         local gz_vcf=${input_vcf} && \
#         local asks=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}')
#     fi

#     if [[ -z ${expected_samples} ]]; then
#         >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): No expected samples specified. Quit checking samples"
#     else
#         local expected_samples=$(echo ${expected_samples} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' | sort - | awk '{printf "%s\n", $1;}')
#         if [[ ${expected_samples} == ${asks} ]]; then
#             >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): Expected samples are identical with actual samples in vcf file ${input_vcf}"
#         else
#             >&2 echo "${input_vcf} has these samples:"
#             >&2 echo "${asks}"
#             >&2 echo "While ${input_vcf} is expected to have these samples:"
#             >&2 echo "${expected_samples}"
#             return 1
#         fi
#     fi
# }

# function check_plain_vcf {
#     local input_vcf=${1}
#     local expected_lines=${2}
#     if [ -z ${expected_lines} ]; then expected_lines=1; fi
#     if [ ! -f ${input_vcf} ]; then return 1; fi

#     if check_vcf_format ${input_vcf}; then
#         >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has solid vcf format. Check whether contains enough record number"
#         if [ $(awk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then
#             >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has enough records."
#         else
#             >&2 echo "${input_vcf} has $(awk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
#             return 3
#         fi
#     else
#         return 2
#     fi
# }

# function check_vcf_format {
#     local input_vcf=${1}
#     local ref_genome="$(dirname $(dirname $0))/ref/ucsc.hg19.fasta"
#     if [ ! -f ${input_vcf} ]; then return 1; fi

#     if [[ ${input_vcf} =~ gz$ ]]; then
#         validate_vcf_index ${input_vcf}
#         if check_gz_file_validity ${input_vcf}; then
#             ${gatk} ValidateVariants \
#             -R ${ref_genome} \
#             -V ${input_vcf} \
#             --validation-type-to-exclude ALL \
#             --verbosity DEBUG || \
#             { >&2 echo "${input_vcf} format is corrupted." && return 1; }
#             return $?
#         else
#             return 2
#         fi
#     else
#         ${gatk} ValidateVariants \
#         -R ${ref_genome} \
#         -V ${input_vcf} \
#         --validation-type-to-exclude ALL \
#         --verbosity DEBUG || \
#         { >&2 echo "${input_vcf} format is corrupted." && return 1; }
#         return $?
#     fi
# }

# function validate_vcf_index {
#     local input_vcf=${1}
#     if [ ! -f ${input_vcf} ]; then >&2 echo "${input_vcf} not existed." && return 1; fi

#     if [[ ${input_vcf} =~ gz$ ]]; then
#         if [ -f ${input_vcf}.tbi ]; then
#             if [ ${input_vcf}.tbi -nt ${input_vcf} ]; then
#                 >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): TBI index file is new. Leave it there."
#             else
#                 >&2 echo "${input_vcf} index too old, trying to rebuild one."
#                 tabix -f -p vcf ${input_vcf}
#             fi
#         else
#             >&2 echo "${input_vcf} index file not existed, build one."
#             tabix -f -p vcf ${input_vcf}
#         fi
#     else
#         >&2 echo "${input_vcf} is plain vcf, does not have to have an index file."
#     fi
# }

# function bcftools_concatvcfs {
#     local OPTIND v o e s
#     while getopts v:o::e::s:: args
#     do
#         case ${args} in
#             v) local input_vcfs=$OPTARG ;;
#             o) local merged_vcf=$OPTARG ;;
#             e) local ignore_error=$OPTARG ;;
#             s) local samples=$OPTARG ;;  # Should be delimited by comma
#             *) echo "No argument passed, Pls at least pass sample path." ;;
#         esac
#     done

#     if [[ ${input_vcfs}  =~ \/ ]] && [[ ! ${input_vcfs} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${input_vcfs} =~ , ]]; then
#         local -a vcfs=($(awk '{printf "%s ", $1;}' < ${input_vcfs}))
#     else
#         local -a vcfs=($(echo ${input_vcfs} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
#     fi
#     >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: the vcfs are ${vcfs[*]}"
#     >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: The merged vcf is ${merged_vcf}"
#     # Check file existence.
#     local -a invalid_vcfs
#     for vcf in "${vcfs[@]}"; do
#         if isValidVCF ${vcf}; then
#             echo -e "Line "${LINENO}": In function "${FUNCNAME}": To be merged ${vcf} is valid."
#         else
#             >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": ${vcf} not existed or corrupted. Run bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh check_vcf_validity ${vcf} to see for yourself"
#             if [[ -z ${ignore_error} ]]; then
#                 return 1
#             else
#                 invalid_vcfs+=( ${vcf} )
#             fi
#         fi
#         if [[ ${vcf} =~ \.vcf$ ]]; then
#             >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Since ${vcf} is plain text format and bcftools -f need to use bgzipped vcfs, we compress the vcf and index it with tabix."
#             bcftools sort -Oz -o ${vcf}.gz ${vcf} && tabix -f -p vcf ${vcf}.gz && \
#             ls -lh ${vcf}.gz
#         fi
#         # if [[ ! -z ${samples} ]]; then
#         #     bcftools view -s "${samples}" -Oz -o ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
#         #     ls -lht ${vcf/.vcf*/.samp.vcf.gz} && \
#         #     check_vcf_validity ${vcf/.vcf*/.samp.vcf.gz} && \
#         #     mv ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
#         #     tabix -f -p vcf ${vcf/.vcf*/.vcf.gz}
#         # fi
#     done

#     if [[ ${#empty_vcfs[@]} -gt 0 ]]; then
#         >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): These files are empty, remove them before concating: ${empty_vcfs[*]}"
#         for ev in "${empty_vcfs[@]}"; do
#             vcfs=( "${vcfs[@]/$ev}" )
#         done
#     fi

#     local tmp_file_list=$(randomID).lst
#     echo "${vcfs[*]}" | awk -F '\t' 'BEGIN{RS=" ";} length($1) > 0{gsub(/\n$/, ""); if ($1 ~ /\.gz$/) printf "%s\n", $1; else printf "%s.gz\n", $1;}' > ${tmp_file_list}
#     >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Here is the temp list file storing the paths of to be concat vcfs:"
#     ls -lh ${tmp_file_list}
#     cat ${tmp_file_list}

#     if [[ ${merged_vcf} =~ \.vcf$ ]]; then
#         local plain_merged=${merged_vcf}
#         local merged_vcf=${merged_vcf}.gz
#     fi

#     # If using file list for input a list of vcfs, each one of them need to be bgzipped and tabix indexed
#     bcftools concat -o ${merged_vcf} -a --no-version -Oz -f ${tmp_file_list} && rm -f ${tmp_file_list} && \
#     bcftools sort -o ${merged_vcf/.vcf.gz/.sorted.vcf.gz} -Oz ${merged_vcf} && \
#     mv ${merged_vcf/.vcf.gz/.sorted.vcf.gz} ${merged_vcf}

#     if [[ ! -z ${plain_merged} ]]; then
#         gunzip -c ${merged_vcf} > ${plain_merged}
#     else
#         tabix -f -p vcf ${merged_vcf}
#     fi
# }

# Alias
alias time="/usr/bin/time -f '$(timestamp) INFO: Elapsed time: %es\tMemory usage: %M KB\tCPU usage: %P'"
