#!/usr/bin/env bash
### This is only a showcase script to display how we downsample the GIAB original BAM files to 30x coverage depth

central_scripts="<path>/ngs_scripts"
source ${central_scripts}/common_bash_utils.sh
self_script=$(realpath ${BASH_SOURCE[0]})
cwd=$(dirname ${self_script})
module load parallel

function prepare_bam_file_per_sample(){
    local raw_bam=${1}
    local target_region=${2}
    local output_bam=${3}
    local total_cpu=${4}
    local ref_genome=${5}

    if [[ -z ${ref_genome} ]]; then
        local ref_genome="<path>/indexed_genome/ucsc.hg19.fasta"
    fi

    local bamID=$(basename ${raw_bam} | cut -f 1 -d ".")
    if [[ ${bamID} == "HG006" ]] || [[ ${bamID} == "HG007" ]]; then
        local raw_bam=${raw_bam/300x/100x}
    fi

    local threads=$(determine_job_num -t ${total_cpu} -c 1)
    
    local freads=$(dirname ${output_bam})/${bamID}.SD.1.fastq.gz
    local rreads=$(dirname ${output_bam})/${bamID}.SD.2.fastq.gz
    local sorted_freads=$(dirname ${output_bam})/${bamID}.SD.1.sorted.fastq.gz
    local sorted_rreads=$(dirname ${output_bam})/${bamID}.SD.2.sorted.fastq.gz
    local ureads=$(dirname ${output_bam})/${bamID}.SD.unpaired.fastq.gz
    local sreads=$(dirname ${output_bam})/${bamID}.SD.singleton.fastq.gz
    local refgen_tag=$(basename ${ref_genome} | cut -f 2 -d ".")
    local tmp_tag=$(randomID)

    log "The input raw bam file is ${raw_bam}"
    log "The target region is ${target_region}"
    log "The output bam file is ${output_bam}"
    log "The total CPU is ${total_cpu}"
    log "The reference genome is ${ref_genome}"

    # Now if target region is a BED file with chr prefix, we need to generate a temp bed file without the chr prefix
    if [[ ${refgen_tag} == "hg19" ]]; then
        display_table ${target_region}
        mawk -F '\t' 'BEGIN{FS=OFS="\t";} \
                    {gsub(/^chr/, "", $1); \
                    print;}' ${target_region} > ${target_region/.bed/.stripchr.${tmp_tag}.bed} && \
        if check_file_content_identity ${target_region/.bed/.stripchr.${tmp_tag}.bed} ${target_region/.bed/.stripchr.bed}; then
            silent_remove_tmps ${target_region/.bed/.stripchr.${tmp_tag}.bed}
        else
            mv ${target_region/.bed/.stripchr.${tmp_tag}.bed} ${target_region/.bed/.stripchr.bed}
        fi
        display_table ${target_region/.bed/.stripchr.bed}
        local target_region=${target_region/.bed/.stripchr.bed}
    fi

    # RAW BAM file is extremely large, we need to downsamp it and select particular regions
    export TMPDIR=<path>/test_tmp

    if [[ ${raw_bam} -nt ${raw_bam}.bai ]]; then
        samtools index -@ ${threads} ${raw_bam}
    fi

    if [[ ${output_bam} -nt ${raw_bam} ]] && \
       [[ ${sorted_freads} -nt ${raw_bam} ]] && \
       [[ ${sorted_rreads} -nt ${raw_bam} ]] && \
       [[ $(file_size ${freads}) -gt $(file_size ${ureads}) ]] && \
       [[ $(file_size ${rreads}) -gt $(file_size ${ureads}) ]] && \
       [[ ${output_bam} -nt ${target_region} ]] && \
       [[ ${output_bam} -nt ${self_script} ]] && \
       check_fastq_pair_consistency ${sorted_freads} ${sorted_rreads} && \
       check_bam_validity ${output_bam}; then
        log "Output file already existed and updated so skip the rest of the function"
        return 0
    else
        log "Output bam file not updated or generated. Execute the rest of the function"
    fi

    local tmp_qname_lst=${TMPDIR}/$(randomID).lst
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The query name list is ${tmp_qname_lst}"

    if [[ ${freads} -nt ${raw_bam} ]] && \
       [[ ${rreads} -nt ${raw_bam} ]] && \
       [[ ${freads} -nt ${target_region} ]] && \
       [[ ${freads} -nt ${self_script} ]]; then
        >&2 echo $'\n\n'"Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): The generated Fastq files are updated comparing the modified time with ${raw_bam}"
    elif [[ ${raw_bam} =~ \.300x\. ]]; then
        sambamba slice -L ${target_region} ${raw_bam} | \
        samtools view -@ ${threads} --subsample-seed 0 --subsample 0.1 | cut -f 1 | sort - | uniq - > ${tmp_qname_lst} || return 1;
        ls -lh ${tmp_qname_lst}
        if [[ $(cat ${tmp_qname_lst} | wc -l) -gt 1 ]]; then
            samtools view -@ ${threads} -h -u -N ${tmp_qname_lst} ${raw_bam} | \
            samtools fastq -@ ${threads} -1 ${freads} -2 ${rreads} -0 ${ureads} -n || return 1;
        fi
        ls -lht ${freads} ${rreads}
    elif [[ ${raw_bam} =~ \.100x\. ]]; then
        sambamba slice -L ${target_region} ${raw_bam} | \
        samtools view -@ ${threads} --subsample-seed 0 --subsample 0.333 | cut -f 1 | sort - | uniq - > ${tmp_qname_lst} || return 1;
        ls -lh ${tmp_qname_lst}
        if [[ $(cat ${tmp_qname_lst} | wc -l) -gt 1 ]]; then
            samtools view -@ ${threads} -h -u -N ${tmp_qname_lst} ${raw_bam} | \
            samtools fastq -@ ${threads} -1 ${freads} -2 ${rreads} -0 ${ureads} -n || return 1;
        fi
        ls -lht ${freads} ${rreads}
    else
        # >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Not sure about the avg depth of the input raw BAM file, check the file: ${raw_bam}"
        sambamba slice -L ${target_region} ${raw_bam} | \
        samtools view -@ ${threads} | cut -f 1 | sort - | uniq - > ${tmp_qname_lst} || return 1;
        ls -lh ${tmp_qname_lst}
        if [[ $(cat ${tmp_qname_lst} | wc -l) -gt 1 ]]; then
            samtools view -@ ${threads} -h -u -N ${tmp_qname_lst} ${raw_bam} | \
            samtools fastq -@ ${threads} -1 ${freads} -2 ${rreads} -0 ${ureads} -n || return 1;
        fi
        ls -lht ${freads} ${rreads}
    fi

    if [[ ${freads} -ot ${self_script} ]] || \
       [[ ${freads} -ot ${raw_bam} ]] || \
       [[ ${freads} -ot ${target_region} ]]; then
        log "Running into error in generating valid fastq files for ${bamID} with ${raw_bam} and ${target_region}"
        return 1
    fi

    silent_remove_tmps ${sreads}

    # Now we just need to adjust sync the read 1 and read2 
    bash ${central_scripts}/common_bash_utils.sh \
    sync_fastq -1 ${freads} -2 ${rreads} -p $(dirname ${output_bam})/${bamID}.SD && \
    ls -lht ${sorted_freads} ${sorted_rreads}

    if check_fastq_pair_consistency ${sorted_freads} ${sorted_rreads}; then
        >&2 echo $'\n\n'"Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${sorted_freads} and ${sorted_rreads} seems to have consistent query names"
    else
        >&2 echo $'\n\n'"Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${sorted_freads} and ${sorted_rreads} seems to have inconsistent query names. Quit with error"
        return 1;
    fi

    # Now we need to align the reads to the reference genome
    if [[ $(dirname ${output_bam})/${bamID}.SD.bam -nt ${freads} ]] && \
       [[ $(dirname ${output_bam})/${bamID}.SD.bam -nt ${rreads} ]] && \
       [[ $(dirname ${output_bam})/${bamID}.SD.bam -nt ${target_region} ]] && \
       [[ $(dirname ${output_bam})/${bamID}.SD.bam -nt ${self_script} ]] && \
       check_bam_validity $(dirname ${output_bam})/${bamID}.SD.bam; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Output BAM already there and updated"
    else
        independent_bwa_alignment \
        -g ${ref_genome} \
        -s ${bamID} \
        -f ${sorted_freads} \
        -r ${sorted_rreads} \
        -o $(dirname ${output_bam})/${bamID}.SD.bam \
        -t ${threads} && \
        ls -lh $(dirname ${output_bam})/${bamID}.SD.bam
    fi

    # Now we need to markduplicates the BAM file
    if [[ $(dirname ${output_bam})/${bamID}.SD.bam -ot ${output_bam} ]] && \
       [[ $(dirname ${output_bam})/${bamID}.SD.bam -ot ${output_bam/.bam/.duplicatesindex} ]] && \
       [[ ${output_bam} -nt ${target_region} ]] && \
       [[ ${output_bam} -nt ${self_script} ]] && \
       check_bam_validity ${output_bam}; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Mark Duplicates BAM already there and updated"
    else
        local gatk="bash ${central_scripts}/common_bash_utils.sh gatk_wrapper"
        ${gatk} \
        --java-options "-Xmx28G -XX:+UseParallelOldGC -XX:+UseAdaptiveSizePolicy -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${threads}" \
        MarkDuplicates \
        -I $(dirname ${output_bam})/${bamID}.SD.bam \
        -O ${output_bam} \
        -M ${output_bam/.bam/.duplicatesindex} \
        -ASO coordinate \
        --CREATE_INDEX \
        --REMOVE_DUPLICATES \
        --VALIDATION_STRINGENCY LENIENT \
        -CO "GATK_MARKDUP_DONE" && \
        samtools index ${output_bam}
    fi
}


function parallel_prepare_bams (){
    local OPTIND b d r c o g
    while getopts b::d::r::c::o::g:: args
    do
        case ${args} in
            b) local bam_dir=$OPTARG ;;
            d) local sd_base=$OPTARG ;;
            r) local total_sd_region=$OPTARG ;;
            c) local total_cpu=$OPTARG ;;
            o) local output_dir=$OPTARG ;;
            g) local ref_genome=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    # To prepare bam file in hg38, you need to specify
    # The reference genome
    # The sd_base
    # The total_sd_region

	if [[ -z ${ref_genome} ]]; then
        local ref_genome="<path>/indexed_genome/ucsc.hg19.fasta"
    fi

    local refgen_tag=$(basename ${ref_genome} | cut -f 2 -d ".")

    if [[ -z ${bam_dir} ]]; then
        local bam_dir="<path>/wgs/GIAB_samples/raw_data/download_data"
    fi

    if [[ -z ${output_dir} ]]; then
        local output_dir="<path>/wgs/GIAB_samples/vcfs/${refgen_tag}"
    fi

    if [[ -z ${sd_base} ]]; then
        local sd_base="<path>/wgs/GIAB_samples/aligned_results/${refgen_tag}"
    fi

    if [[ -z ${total_sd_region} ]]; then
		if [[ ${refgen_tag} == "hg19" ]]; then
        	local total_sd_region=<path>/public_data/SD_from_SEDEF/hg19/test_BISER/WGAC.hg19.cigar.trimmed.homo.expanded.merged.bed
		elif [[ ${refgen_tag} == "hg38" ]]; then
			local total_sd_region="<path>/public_data/SD_from_SEDEF/hg19/test_BISER/WGAC.hg38.cigar.trimmed.homo.expanded.merged.bed"
		fi
    fi

    if [[ -z ${total_cpu} ]]; then
        local total_cpu=6
    fi

    if [[ ${refgen_tag} == "hg19" ]]; then
        local original_bam_temp="hs37d5.300x.bam"
    elif [[ ${refgen_tag} == "hg38" ]]; then
        local original_bam_temp="GRCh38.300x.bam"
    else
        log "Cannot extract the reference genome tag with ${ref_genome}"
        return 1
    fi

    local -a samples=( "HG002" "HG003" "HG004" "HG005" "HG006" "HG007" )

    parallel -j6 --dry-run \
    bash ${self_script} prepare_bam_file_per_sample \
    ${bam_dir}/{}.${original_bam_temp} \
    ${total_sd_region} \
    ${sd_base}/{}.SD.deduped.bam \
    ${total_cpu} \
    ${ref_genome} '>' $(dirname ${self_script})/{}.prepare_bam.refSD.${refgen_tag}.log '2>&1' ::: "${samples[@]}" && \
    parallel -j6 --tmpdir <path>/test_tmp \
    --joblog $(dirname ${self_script})/parallel_prepare_bam_GIAB_samples.refSD.${refgen_tag}.log \
    bash ${self_script} prepare_bam_file_per_sample \
    ${bam_dir}/{}.${original_bam_temp} \
    ${total_sd_region} \
    ${sd_base}/{}.SD.deduped.bam \
    ${total_cpu} \
    ${ref_genome} '>' $(dirname ${self_script})/{}.prepare_bam.refSD.${refgen_tag}.log '2>&1' ::: "${samples[@]}"

    check_parallel_joblog $(dirname ${self_script})/parallel_prepare_bam_GIAB_samples.refSD.${refgen_tag}.log

}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    declare -a func_names=($(typeset -f | awk '!/^main[ (]/ && /^[^ {}]+ *\(\)/ { gsub(/[()]/, "", $1); printf $1" ";}'))
    declare -a input_func_names=($(return_array_intersection "${func_names[*]}" "$*"))
    declare -a arg_indices=($(get_array_index "${input_func_names[*]}" "$*"))
    if [[ ${#input_func_names[@]} -gt 0 ]]; then
        >&2 echo "Seems like the command is trying to directly run a function ${input_func_names[*]} in this script ${BASH_SOURCE[0]}. The defined functions are ${func_names[*]}"
        >&2 echo "The input arguments are ${*}"
        first_func_ind=${arg_indices}
        >&2 echo "The identified first func name is at the ${first_func_ind}th input argument, while the total input arguments are: $*"
        following_arg_ind=$((first_func_ind + 1))
        >&2 echo "Executing: ${*:${following_arg_ind}}"
        "${@:${following_arg_ind}}"
    else
        set -T
        set -E
        set -o pipefail
        trap 'PreCommand ${LINENO} ${BASH_SOURCE} 2> /dev/null' DEBUG
        trap 'catch_exit_status $?' EXIT
        trap 'catch_error_report $? ${LINENO} ${BASH_COMMAND}' ERR
        return_to_conda_base

        # Run the main part
        parallel_prepare_bams "$@"
    fi
fi