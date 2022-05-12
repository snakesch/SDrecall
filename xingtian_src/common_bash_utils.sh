#!/bin/bash
central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts

function timestamp () {
    date
}

function echo_line_no () {
    grep -n "$1" $0 |  sed "s/echo_line_no//" 
    # grep the line(s) containing input $1 with line numbers
    # replace the function name with nothing 
}

function display_table () {
    if [[ -z $2 ]]; then local rows=4; else local rows=${2}; fi
    if [[ -z $3 ]]; then local delimiter="\t"; else local delimiter=${3}; fi

    python3 /paedyl01/disk1/yangyxt/ngs_scripts/display_table.py \
    -df $1 \
    -nr ${rows}
}

# Define a fucntion to judge if some element in an array
function containsElement () {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}

function randomID () {
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}


function catch_error_report () {
    local error_code=${1}
    local line_no=${2}
    local command_exec=${3}
    local email_addr='u3005579@connect.hku.hk'
    local err_message=$(perror ${error_code} | head -1)

    >&2 echo $(timestamp)": Hitting an error (code is "${error_code}"): Error message is "${err_message}

    local tmp_email_content=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).email.tmp
    local func_chains=$(awk -v all_funcs="${FUNCNAME[*]}" 'BEGIN{split(all_funcs, funcs, " "); \
                        for (i in funcs) printf "---> %s ", funcs[i];}')
    for func in "${FUNCNAME[@]}"; do 
        if [[ ${func} != "catch_error_report" ]]; then local exec_func=${func} && break; fi
    done
    for script in "${BASH_SOURCE[@]}"; do
        if [[ $(basename ${script}) != "common_bash_utils.sh" ]]; then local script_name=${script} && break; fi
    done
    local script_path="/paedyl01/disk1/yangyxt/ngs_scripts/"$(basename "${script_name}")
    if [ -f ${script_path} ]; then
        local exec_command=$(awk 'NR == '${line_no}'{print;}' < ${script_path})
    else
        local exec_command="Not in a script file written in-house. Quit trying to fetch the code line."
    fi
    echo "$(timestamp): The script causing this error is ${script_name}."
    echo "$(timestamp): All the invoked bash scripts are ${BASH_SOURCE[*]}."

    if [[ ${#sample_ID} -gt 0 ]]; then
        local sampID=${sample_ID}
    elif [[ ${#family_name} -gt 0 ]]; then
        local sampID=${family_name}
    fi
    
    awk -v ec="${error_code}" -v ln="${line_no}" -v sp="${script_path}" \
    -v ef="${exec_func}" -v fc="${func_chains}" -v ss="${BASH_SOURCE[*]}" \
    -v ecom="${exec_command}" -v envcom="${command_exec}" -v sid="${sampID}" \
    'BEGIN{printf "Dear Yangyxt,\n\
    An Error occurred (code %s) on line %s in the script file: %s.\n\
    The error happened in PBS job ID '${PBS_JOBID}'.\n\
    The sample ID (can be a FAM_ID or INDIVIDUAL_ID) being processed is %s.\n\
    The function being executed is %s. And the whole function chain being called are:\n\
    %s.\n\
    The scripts invoked are:\n\
    %s.\n\
    The specific command causing this error deduced by the LINENO variable is:\n\
    %s\n\
    The specific command causing error according to BASH_COMMAND variable is:\n\
    %s", ec, ln, sp, sid, ef, fc, ss, ecom, envcom;}' > ${tmp_email_content} && \
    mail -s "Error_when_dealing_with_${sample_batch}_in_script_$(basename ${script_path})" ${email_addr} < ${tmp_email_content} || { echo "Failed to send email" && awk '{print;}' < ${tmp_email_content} ;}

    awk '{print;}' < ${tmp_email_content}
    rm -f ${tmp_email_content}
    exit ${error_code} 
}


function catch_exit_status () {
    trap 'echo ${BASH_COMMAND} failed in line ${LINENO}' ERR
    local exit_code=${1}
    local email_cont_path=${2}

    for script in "${BASH_SOURCE[@]}"; do
        if [[ $(basename ${script}) != "common_bash_utils.sh" ]]; then local script_name=${script} && break; fi
    done
    
    local script_path=$(cd "$(dirname "${script_name}")" && pwd)/$(basename "${script_name}")

    if [[ ${exit_code} -eq 126 ]]; then
        >&2 echo $(timestamp)": Hitting an error causing script to exit. Exit status implied the command invoked cannot be executed. Could be a permission problem or command is not executable."
    elif [[ ${exit_code} -eq 127 ]]; then
        >&2 echo $(timestamp)": The exit status 127 stands for a typical error, command not found."
    elif [[ ${exit_code} -eq 128 ]]; then
        >&2 echo $(timestamp)": The exit status 128 stands for an illegal argument specified for a command."
    elif [[ ${exit_code} -eq 130 ]]; then
        >&2 echo $(timestamp)": The exit status 130 usually comes from the user's ctrl+c."
    elif [[ ${exit_code} -gt 128 ]]; then
        local fatal_error_code=$(( exit_code - 128 ))
        >&2 echo $(timestamp)": The script ${script_path} hit a fatal error with error code ${fatal_error_code}."
    elif [[ ${exit_code} -gt 0 ]]; then
        >&2 echo $(timestamp)": Hit an unknown kind error. Directly email report."
    elif [[ ${exit_code} -eq 0 ]]; then
        >&2 echo $(timestamp)": Whole script ${BASH_SOURCE} has been executed without an error."
        exit 0
    fi

    exit ${exit_code}
}

function send_email_to_clinicians () {
    local OPTIND r s c
    while getopts r::s::c:: args
    do 
        case ${args} in
            r) local recepients=$OPTARG ;; 
            s) local subject=$OPTARG ;;
            c) local log_content=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${recepients} ]]; then local recepients="dan.leung@connect.hku.hk,jingy09hk@gmail.com,jingy09@hku.hk,u3005579@connect.hku.hk"; fi
    local sample_batch=$(echo ${subject} | awk -F '_' '{printf $1;}')
    local subject=\"${subject}\"

    awk -v sb="${sample_batch}" 'BEGIN{printf "Dear Dr Jing and Daniel,\n"; \
    printf "The Bioinformatic processing of batch %s is finished. Here are the gdrive links for your reference:\n\n", sb;}' > ${log_content}.tmp

    cat ${log_content} >> ${log_content}.tmp && \
    mail -s ${subject} ${recepients} < ${log_content}.tmp || echo "Failed to send emails."

    trap "rm -f ${tmp_email_content}" RETURN
}


function check_return_code () {
    local return_code=$(echo $?)
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The last step's return code is "${return_code}
    if [[ ${return_code} -gt 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}:$(timestamp)":The last step does not finish in normal way. Exiting the whole script now." && \
        exit 1
    else
        >&2 echo $(timestamp)": The last step finished properly. Continue"; fi
}

function build_genomicsDB_from_scratch () {
    local genomicDB=${1}
    local chr=${2}
    local sample_map=${3}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    if [[ -d ${genomicDB} ]]; then rm -R -f ${genomicDB}; fi

    local batch_size=$(wc -l ${sample_map} | awk '{printf $1}')
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The batch size for probe ${probe} is ${batch_size}

    # After running the script above, all gvcfs corresponding to each sample should be already there. Now we generate a script, which is used to generate the script for consolidated genotyping
    time ${gatk} --java-options "-Xmx8g -Xms2g" GenomicsDBImport \
        --tmp-dir /paedyl01/disk1/yangyxt/test_tmp \
        --genomicsdb-workspace-path ${genomicDB} \
        -R ${ref_gen}/ucsc.hg19.fasta \
        --batch-size ${batch_size} \
        --sample-name-map ${sample_map} \
        --reader-threads 5 \
        --intervals chr${chr}
    check_return_code
}

function get_parent_dir () {
    local input=${1}
    if [[ -f ${input} ]]; then
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    elif [[ -d ${input} ]]; then
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    elif [[ ${input} =~ \.[a-z]+$ ]]; then
        cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
        echo ${cd}
    else
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    fi
}

function gnws () {
    # gnws stands for get name of the file without suffix. (still with the path)
    local input=${1}
    echo ${input} | awk -F '.' '{for(i=1;i<NF;i++) printf $i".";}' | awk '{gsub(/\.$/, ""); print}'
}


function get_name () {
    # Get name of the file without suffix and without path
    local input=${1}

    if [[ -d ${input} ]]; then
        local last=$(echo ${input} | awk -F '/' '{printf $NF;}')
        if [[ ${#last} -eq 0 ]]; then
            echo ${input} | awk -F '/' '{printf $(NF-1);}'
        else
            echo ${input} | awk -F '/' '{printf $NF;}'
        fi
    elif [[ -f ${input} ]]; then
        echo ${input} | awk -F '/' '{printf $NF;}' | awk -F '.' '{for(i=1;i<NF;i++) printf $i".";}' | awk '{gsub(/\.$/, ""); print}'
    fi
}

function join_by () { local IFS="$1"; shift; echo "$*"; }

# Define a fucntion to judge if some element in an array
function containsElement () {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}


function drop_redundant_rows () {
    export TMPDIR=/paedyl01/disk1/yangyxt/test_tmp/
    >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): Before dropping duplicated rows in ${1}, the table now is $(ls -lh ${1})"
    python3 -c 'import pandas as pd; pd.read_table("'${1}'", low_memory=False).dropna(how="all").drop_duplicates().to_csv("'${1}'", sep="\t", index=False);'
    >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): After dropping duplicated rows in ${1}, the table now is $(ls -lh ${1})"
}


# Define a function that check the size of a table once it below 2k the whole script exit 1.
function check_result_size() {
    local size=$(ls -l ${1} | awk '{gsub(/ +/, " "); print}' | awk '{print $5}')
    local lines=$(wc -l ${1} | awk '{print $1}')
    if [[ ${size} -le 2800 ]] || [[ ${lines} -le 1 ]]; then 
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The table "$1" after processing is just too oddly small, call it an end for now."
        ls -lh $1
        if [[ -z ${2} ]]; then
            exit 1
        else
            return 0
        fi
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The table "$1" after processing is not too small to only fit a header line. Continue"
        ls -lh $1
        if [[ -z ${2} ]]; then
            if [[ $1 =~ (\.bed$|\.txt$|\.ped$|\.tsv$|\.csv$|\.table$|\.tmp$|\.sam$) ]] && [[ ! $1 =~ \.vcf ]]; then display_table $1; fi
        else
            return 0
        fi
    fi
}


function check_gz_file_validity () {
    local input=${1}
    if [[ ${input} =~ \.[b]*gz$ ]]; then
        local validity=$(gzip -v -t ${input} 2>&1 | awk 'END{printf $NF;}')
        local size=$(ls -l ${input} | awk '{gsub(/[ ]+/, " "); print;}' | awk '{printf $5}')
        
        if [[ ${validity} == "OK" ]] && [[ ${size} -gt 30 ]]; then
            return 0
        elif [[ ${validity} == "OK" ]]; then
            >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} is validily compressed but the size is oddly small: ${size}"
            return 1
        elif [[ ${size} -gt 30 ]]; then
            >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} is not oddly small: ${size}. But it's not a valid gzipped file."
            return 1
        else
            >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} is neither validily compressed nor with a normal file size: ${size}"
            return 1
        fi
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp) Input file is not gzipped."
        return 1
    fi
}


function check_fq_file_validity() {
    local input_fq=${1}
    local expected_lines=${2}

    if [[ -z ${expected_lines} ]]; then local expected_lines=100; fi

    if [[ ${input_fq} =~ \.gz$ ]]; then
        local p_fq=${input_fq/.gz/}
        local gz_fq=${input_fq}
    else
        local p_fq=${input_fq}
        local gz_fq="${input_fq}.gz"
    fi

    local p_lines=0
    local gz_lines=0

    local p_fq_validity=1
    local gz_fq_validity=1

    if [[ -f ${p_fq} ]]; then
        local actual_lines=$(cat ${p_fq} | wc -l)
        local p_lines=${actual_lines}
        >&2 echo $(timestamp)": This fastq file ${p_fq} actually has ${actual_lines} lines"
        if [[ ${actual_lines} -gt ${expected_lines} ]] && [[ $((actual_lines % 4)) -eq 0 ]]; then
            local p_fq_validity=0
        elif [[ ${actual_lines} -gt ${expected_lines} ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${p_fq} exists, but it contains ${actual_lines} lines, remainder of 4 is not zero."
        elif [[ $((actual_lines % 4)) -eq 0 ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${p_fq} exists, but the length of it does not meet expectation. It only contains ${actual_lines} lines and is expected to have ${expected_lines} lines."
        fi
    fi

    if [[ -f ${gz_fq} ]]; then
        local actual_lines=$(zcat ${gz_fq} | wc -l)
        local gz_lines=${actual_lines}
        >&2 echo $(timestamp)": This fastq file ${gz_fq} actually has ${actual_lines} lines"
        if check_gz_file_validity ${gz_fq} && [[ ${actual_lines} -gt ${expected_lines} ]] && [[ $((actual_lines % 4)) -eq 0 ]]; then
            local gz_fq_validity=0
        elif [[ ${actual_lines} -gt ${expected_lines} ]] && [[ $((actual_lines % 4)) -eq 0 ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${gz_fq} exists, and its big enough, it contains ${actual_lines} lines, but it does not have valid gzip compression."
        elif check_gz_file_validity ${gz_fq} && [[ ${actual_lines} -gt ${expected_lines} ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${gz_fq} exists, but it contains ${actual_lines} lines, remainder of 4 is not zero."
        elif check_gz_file_validity ${gz_fq} && [[ $((actual_lines % 4)) -eq 0 ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): ${gz_fq} exists, but the length of it does not meet expectation. It only contains ${actual_lines} lines and is expected to have ${expected_lines} lines."
        fi
    fi

    if [[ ${input_fq} == ${p_fq} ]]; then
        if [[ ${p_fq_validity} -eq 0 ]] && [[ ${gz_fq_validity} != "0" ]]; then
            return 0
        elif [[ ${p_fq_validity} -eq 0 ]] && [[ ${gz_fq_validity} -eq 0 ]]; then
            if [[ ${p_lines} -ge ${gz_lines} ]]; then
                return 0
            else
                gunzip -f -c ${gz_fq} > ${input_fq} && \
                return 0
            fi
        elif [[ ${gz_fq_validity} -eq 0 ]]; then
            gunzip -f -c ${gz_fq} > ${input_fq} && \
            return 0
        else
            return 1
        fi
    elif [[ ${input_fq} == ${gz_fq} ]]; then
        if [[ ${p_fq_validity} != "0" ]] && [[ ${gz_fq_validity} -eq 0 ]]; then
            return 0
        elif [[ ${p_fq_validity} -eq 0 ]] && [[ ${gz_fq_validity} -eq 0 ]]; then
            if [[ ${p_lines} -le ${gz_lines} ]]; then
                return 0
            else
                gzip -f -c ${p_fq} > ${input_fq} && \
                return 0
            fi
        elif [[ ${p_fq_validity} -eq 0 ]]; then
            gzip -f -c ${p_fq} > ${input_fq} && \
            return 0
        else
            return 1
        fi
    fi
}


function quick_check_bam_validity () {
    local input=${1}
    samtools quickcheck -v ${input} && return 0 || return 1
}


function check_bam_qname_format() {
    local input_bam=${1}
    local threads=${2} 

    if [[ -z ${threads} ]]; then local threads=1; fi

    local issue_lines=$(samtools view ${input_bam} | head -20 | awk -F '\t' '$1 ~ /\/[1-2]$/{print;}' | wc -l)
    if [[ ${issue_lines} -eq 0 ]]; then
        return 0
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: The bam file ${input_bam} contains problematic query name format: "
        { samtools view ${input_bam} | head -10
        samtools view -h ${input_bam} | \
        mawk -F '\t' '{gsub(/\/[1-2]$/, "", $1); print;}' | \
        samtools sort -O bam -o ${input_bam/.bam/.temp.bam} -@ ${threads} - ;} && \
        samtools index ${input_bam/.bam/.temp.bam} && \
        mv ${input_bam/.bam/.temp.bam} ${input_bam} && \
        mv ${input_bam/.bam/.temp.bam.bai} ${input_bam}.bai && \
        { rm -f ${input_bam/.bam/.bai} || >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: ${input_bam/.bam/.bai} not existed, does not have to delete it."; } && \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Manually remove the hanging /1 or /2 tag at the end of query names. New bam is $(ls -lh ${input_bam})" && \
        return 0 || \
        { >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Failed to manually modify the query name format on ${input_bam}." && return 1; }
    fi
}

function check_bam_index() {
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


function check_bam_validity () {
    local input=${1}
    local expected_lines=${2}
    if [[ -z ${expected_lines} ]]; then local expected_lines=10; fi
    local gatk=$HOME/software/gatk-4.2.2.0/gatk

    if [[ ! -f ${input} ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, found ${input} not even exist."
        return 1
    fi

    local line_no=$( samtools view -c ${input} )

    if [[ ${line_no} -lt ${expected_lines} ]]; then 
        >&2 echo "$(timestamp): In function ${FUNCNAME}, found ${input} only has ${line_no} lines where we expect it to have ${expected_lines} lines."
        return 1
    fi
    if check_bam_qname_format ${input}; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${input} has proper query name format."
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${input} has inappropriate query name format."
        return 1
    fi

    check_bam_index ${input}

    ${gatk} ValidateSamFile \
    -I ${input} \
    --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE \
    -MO 1000000 > ${input::-4}.check

    local -a err_types=($(mawk '$1 ~ /^ERROR::/{print $1;}' ${input::-4}.check | sort - | uniq - | awk -F ':' '{printf "%s ",$3;}'))
    if [[ ${#err_types[@]} -eq 0 ]]; then
        return 0
    else
        for err in "${err_types[@]}"; do
            if [[ ! ${err} =~ "INVALID_INDEX_FILE_POINTER" ]] && [[ ! ${err} =~ "MATE_NOT_FOUND" ]] && [[ ! ${err} =~ "INVALID_PLATFORM_VALUE" ]]; then
                echo "In function ${FUNCNAME}: $(timestamp): ${input} is corrupted and cannot be simply fixed."
                return 1
            elif [[ ${err} =~ "INVALID_INDEX_FILE_POINTER" ]]; then
                echo "In function ${FUNCNAME}: $(timestamp): ${input} has corrupted index file, re-indexing it."
                #First delete the old index file
                rm -f ${input::-1}i || echo "In function ${FUNCNAME}: $(timestamp): Index file ${input::-1}i already deleted."
                rm -f ${input}.bai || echo "In function ${FUNCNAME}: $(timestamp): Index file ${input}.bai already deleted."
                samtools index ${input}
            elif [[ ${err} =~ "MATE_NOT_FOUND" ]]; then
                echo "In function ${FUNCNAME}: $(timestamp): ${input} has some singleton reads. Just leave them there."
            elif [[ ${err} =~ "INVALID_PLATFORM_VALUE" ]]; then
                echo "In function ${FUNCNAME}: $(timestamp): ${input} might be a simulated golden bam file, which does not have a valid PL value. Let it pass. The error message is:"
                mawk '$1 ~ /^ERROR::INVALID_PLATFORM_VALUE/{print $0;}' ${input::-4}.check
            fi
        done
        # IF the for loop above is executed without return 1, then no other error types are identified. return 0
        return 0
    fi
    trap "rm ${input::-4}.check || echo ${input::-4}.check alread removed" RETURN
}


function validate_vcf_index() {
    local input_vcf=${1}
    if [ ! -f ${input_vcf} ]; then >&2 echo "${input_vcf} not existed." && return 1; fi

    if [[ ${input_vcf} =~ gz$ ]]; then
        if [ -f ${input_vcf}.tbi ]; then
            if [ ${input_vcf}.tbi -nt ${input_vcf} ]; then
                return 0
            else
                >&2 echo "${input_vcf} index too old, trying to rebuild one."
                tabix -f -p vcf ${input_vcf}
            fi
        else
            >&2 echo "${input_vcf} index file not existed, build one."
            tabix -f -p vcf ${input_vcf}
        fi
    else
        >&2 echo "${input_vcf} is plain vcf, does not have to have an index file."
        return 0
    fi
}


function check_vcf_format() {
    local input_vcf=${1}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
    if [ ! -f ${input_vcf} ]; then return 1; fi
    if [ -z ${2} ]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; else local ref_genome=${2}; fi

    if [[ ${input_vcf} =~ gz$ ]]; then
        validate_vcf_index ${input_vcf}
        if check_gz_file_validity ${input_vcf}; then
            ${gatk} ValidateVariants \
            -R ${ref_genome} \
            -V ${input_vcf} \
            --validation-type-to-exclude ALL \
            --verbosity DEBUG || \
            { >&2 echo "${input_vcf} format is corrupted." && return 1; }
            return $?
        else
            return 2
        fi
    else
        ${gatk} ValidateVariants \
        -R ${ref_genome} \
        -V ${input_vcf} \
        --validation-type-to-exclude ALL \
        --verbosity DEBUG || \
        { >&2 echo "${input_vcf} format is corrupted." && return 1; }
        return $?
    fi
}


function check_plain_vcf() {
    local input_vcf=${1}
    local expected_lines=${2}
    if [ -z ${expected_lines} ]; then expected_lines=1; fi
    if [ ! -f ${input_vcf} ]; then return 1; fi

    if check_vcf_format ${input_vcf}; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has solid vcf format. Check whether contains enough record number"
        if [ $(mawk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then
            >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has enough records."
            return 0
        else
            >&2 echo "${input_vcf} has $(mawk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
            return 3
        fi
    else
        return 2
    fi
}


function check_gz_vcf(){
    local input_vcf=${1}
    local expected_lines=${2}
    if [ -z ${expected_lines} ]; then expected_lines=1; fi
    if [ ! -f ${input_vcf} ]; then return 1; fi

    if check_vcf_format ${input_vcf}; then  # Check input_vcf format
        >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has solid vcf format. Check whether contains enough record number"
        if [ $(zcat ${input_vcf} | mawk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then  # Check input_vcf content
            >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} has enough records."
            return 0
        else
            >&2 echo "${input_vcf} has $(zcat ${input_vcf} | mawk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
            return 3
        fi
    else
        return 2
    fi
}


function check_vcf_samples(){
    local input_vcf=${1}
    local expected_samples=${2}
    # expected samples should be delimited by comma

    if [[ ${input_vcf} =~ \.vcf$ ]]; then
        bgzip -f -c ${input_vcf} > ${input_vcf/.vcf/.tmp.vcf.gz} && \
        tabix -f -p vcf ${input_vcf/.vcf/.tmp.vcf.gz} && \
        local gz_vcf=${input_vcf/.vcf/.tmp.vcf.gz} && \
        local as=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}') && \
        rm -f ${gz_vcf} || { >&2 echo "Failed to compress ${input_vcf}" && return 1; }
        if [[ $? != "0" ]]; then return 1; fi 
    else
        local gz_vcf=${input_vcf} && \
        local as=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}')
    fi

    if [[ -z ${expected_samples} ]]; then
        return 0
    else
        local es=$(echo ${expected_samples} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' | sort - | awk '{printf "%s\n", $1;}')
        if [[ ${es} == ${as} ]]; then
            return 0
        else
            >&2 echo "${input_vcf} has these samples:"
            >&2 echo "${as}"
            >&2 echo "While ${input_vcf} is expected to have these samples:"
            >&2 echo "${es}"
            return 1
        fi
    fi
}


function check_vcf_refs () {
    local input_vcf=${1}
    local fix=${2}
    local ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"

    if [[ ! -f ${input_vcf} ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, ${input_vcf} not exist."
        return 1
    fi

    module load bcftools
    local result=$(bcftools filter -i 'TYPE="indel" || TYPE="snp"' ${input_vcf} | bcftools norm -c e -f ${ref_genome} - 2>&1 | tail -1)
    if [[ ${result} =~ Reference\ allele\ mismatch\ at ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, ${input_vcf} has False REF alleles: ${result}"
        if [[ ${fix} == "fix" ]]; then
            >&2 echo "$(timestamp): In function ${FUNCNAME}, directly fix false REF alleles in ${input_vcf} with command: bcftools norm -c ws -f ${ref_genome} -Oz -o ${input_vcf/.vcf*/.tmp.vcf.gz} ${input_vcf}"
            bcftools norm -c ws -f ${ref_genome} -Oz -o ${input_vcf/.vcf*/.tmp.vcf.gz} ${input_vcf}
            module unload bcftools
            if [[ ${input_vcf} =~ \.gz ]]; then
                mv ${input_vcf/.vcf*/.tmp.vcf.gz} ${input_vcf} && \
                tabix -f -p vcf ${input_vcf}
            else
                gunzip -f -c ${input_vcf/.vcf*/.tmp.vcf.gz} > ${input_vcf} && \
                rm -f ${input_vcf/.vcf*/.tmp.vcf.gz}
            fi
        else
            >&2 echo "$(timestamp): In function ${FUNCNAME}, Do not fix ${input_vcf} false REF alleles"
            module unload bcftools
            return 1
        fi
    else
        >&2 echo "$(timestamp): In function ${FUNCNAME}, ${input_vcf} has valid REF alleles: ${result}"
        module unload bcftools
        return 0
    fi
}


function check_vcf_validity () {
    local input_vcf=${1}
    local expected_lines=${2}
    local expected_samples=${3}

    if [ -z ${expected_lines} ]; then expected_lines=1; fi
    if [ ! -f ${input_vcf} ]; then 
        >&2 echo "$(timestamp): In function ${FUNCNAME}, ${input_vcf} does not even exist." 
        return 1
    fi

    if [[ ${input_vcf} =~ \.gz$ ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} should be bgzipped, check gz vcf validity"
        if check_gz_vcf ${input_vcf} ${expected_lines}; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                return 0
            else
                return 1
            fi
        elif check_plain_vcf ${input_vcf::-3} ${expected_lines} && [[ ${input_vcf::-3} -nt ${input_vcf} ]]; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                bgzip -f -c ${input_vcf::-3} > ${input_vcf} && \
                tabix -f -p vcf ${input_vcf} && \
                return 0 || return 1
            else
                return 1
            fi
        else
            return 1
        fi
    elif [[ ${input_vcf} =~ \.vcf$ ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}: ${input_vcf} should be plain text, check plain vcf validity"
        if check_plain_vcf ${input_vcf} ${expected_lines}; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                return 0
            else
                return 1
            fi
        elif check_gz_vcf ${input_vcf}.gz ${expected_lines} && [[ ${input_vcf}.gz -nt ${input_vcf} ]]; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                gunzip -f -c ${input_vcf}.gz > ${input_vcf} && \
                return 0 || return 1
            else
                return 1
            fi
        else
            return 1
        fi
    else
        return 2
    fi
}


function check_odd_tab() {
    local -a odd_tabs=($(bash /paedyl01/disk1/yangyxt/ngs_scripts/check_odd_tab_rows.sh $1))
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: odd_tabs = ${odd_tabs[*]}"
    if [[ ${#odd_tabs[@]} -ge 2 ]]
    then
        sed -i 's/\t*$//' $1
        local -a second_odd=($(bash /paedyl01/disk1/yangyxt/ngs_scripts/check_odd_tab_rows.sh $1))
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: second_odd = ${second_odd[*]}"
        if [[ ${#second_odd[@]} -ge 2 ]]
        then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The table "$1" has some odd rows with errorenous tabs."
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The rows having odd tabs are ${second_odd[*]}"
            exit 1
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The table"$1" is normal and continue."
            ls -lh $1
        fi
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The table"$1" is normal and continue."
        ls -lh $1
    fi 
}

function drop_duplicated_cols () {
    local table=$1
    local -a x_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /_x$/) printf i" ";} }' < ${table}))
    local -a y_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /_y$/) printf i" ";} }' < ${table}))
    local -a extra_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /\.1$/) printf i" ";} }' < ${table}))
    local -a tmp_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /^uniq_ID/) printf i" ";} }' < ${table}))

    if [[ ${#x_cols[@]} -eq ${#y_cols[@]} ]]; then
        if [[ ${#x_cols[@]} -gt 0 ]]; then
            local tobe_dropped_cols=$(echo "${y_cols[*]}" | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
            cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && \
            awk -F '\t' '{ if(NR==1) {for (i=1;i<NF;i++) {gsub(/_x$/, "", $i); printf $i"\t";} printf $NF"\n";} else print; }' < ${table}.tmp > ${table}
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "We dont have duplicated columns for this table: "${table}
            return 0
        fi
    elif [[ ${#x_cols[@]} -gt 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Seems that the _y suffix cols has been removed, need to remove _x suffix for those labels"
        awk -F '\t' '{ if(NR==1) {for (i=1;i<NF;i++) {gsub(/_x$/, "", $i); printf $i"\t";} printf $NF"\n";} else print; }' < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "We seems to detect duplicated rows but the number between _x suffix and _y suffix cols is not consistent. Check the table."
        display_table ${table}
    fi

    # if [[ ${#extra_cols[@]} -gt 0 ]]; then
    #     local tobe_dropped_cols=$(echo "${extra_cols[*]}" | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
    #     cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
    # fi

    # if [[ ${#tmp_cols[@]} -gt 0 ]]; then
    #     local tobe_dropped_cols=$(echo "${tmp_cols[*]}" | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
    #     cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
    # fi
}


function drop_certain_cols() {
    local table=${1}
    local -a cols=($(echo ${2} | mawk 'BEGIN{RS=",";} {printf "%s ", $1;}'))
    local output_table=${3}

    if [[ ${cols} =~ ^[0-9]+$ ]]; then
        local -a col_inds=("${cols[@]}")
    else
        local -a col_inds
        for col in "${cols[@]}"; do
            local col_ind=$(head -1 ${table} | mawk 'BEGIN{RS="\t";} $1 == "'${col}'"{gsub(/\n$/, ""); printf "%s", NR; exit 0;}')
            if [[ ${col_ind} =~ ^[0-9]+$ ]]; then
                col_inds+=( "${col_ind}" )
            fi
        done
    fi

    >&2 echo "$(timestamp): Line "${LINENO}": In function "${FUNCNAME}": ${cols[*]} respective column indices are ${col_inds[*]}"
    local col_args=$(echo "${col_inds[@]}" | mawk 'BEGIN{RS=" ";} { gsub(/\n$/, ""); printf "%s,", $1;}' | awk '{gsub(/,$/, ""); print;}')
    
    if [[ ${#col_inds[@]} -gt 0 ]]; then
        if [[ -z ${output_table} ]]; then
            cut -f ${col_args} --complement ${table} | uniq - > ${table}.tmp && \
            mv ${table}.tmp ${table} && ls -lh ${table} && \
            echo "$(timestamp): Line "${LINENO}": In function "${FUNCNAME}": ${table} has been removed these columns: ${col_inds[*]}."
        else
            cut -f ${col_args} --complement ${table} | uniq - > ${output_table} && \
            echo "$(timestamp): Line "${LINENO}": In function "${FUNCNAME}": ${output_table} is ${table} without these columns: ${col_inds[*]}."
        fi
    else
        >&2 echo "$(timestamp): Line "${LINENO}": In function "${FUNCNAME}": Cant find these columns : ${cols[*]} in table ${table}. Directly skip this function"
    fi
}


function awk_merge_two_tables () {
    # We assume that two tables are with headers.
    # Better merge two tables on a common sharing column with the same header label.
    local left_table=$1
    local right_table=$2
    local left_on_column=$3
    local right_on_column=$4
    local merged_table=$5
    local delimiter=$6

    if [[ -z ${delimiter} ]]; then local delimiter="\t"; fi

    if [[ ${left_on_column} =~ ^[0-9]+$ ]]
    then
        local left_col_ind=${left_on_column}
        local right_col_ind=${right_on_column}
    elif [[ ${left_on_column} == ${right_on_column} ]]
    then
        local left_col_ind=$(head -1 ${left_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
        local right_col_ind=$(head -1 ${right_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${right_on_column}'") printf i;}}')
    elif [[ ${left_on_column} != ${right_on_column} ]]
    then
        awk -F '\t' 'BEGIN{OFS=FS="\t";} {if (NR == 1) {gsub("'${right_on_column}'", "'${left_on_column}'"); print;} else print;}' < ${right_table} > ${right_table}.tmp && \
        mv ${right_table}.tmp ${right_table}
        check_return_code
        check_result_size ${right_table}
        local left_col_ind=$(head -1 ${left_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
        local right_col_ind=$(head -1 ${right_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: Left Table is $(ls -lh ${left_table})
    display_table ${left_table}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: Right Table is $(ls -lh ${right_table})
    display_table ${right_table}
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Left table target column is ${left_on_column} while the column index is ${left_col_ind}"
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Right table target column is ${right_on_column} while the column index is ${right_col_ind}"

    awk -v d="${delimiter}" 'BEGIN {FS=OFS=d;} \
    NR == FNR {aa[$'${right_col_ind}'] = $0;} \
    NR > FNR {print $0, aa[$'${left_col_ind}'];} \
    ' ${right_table} ${left_table} > ${merged_table}

    check_return_code
    # check_result_size ${merged_table}
}


function pandas_merge_table(){
    local OPTIND l r o a b t h c
    while getopts l:r:c::a::b::t::h::o: args
    do 
        case ${args} in
            l) local left_table=$OPTARG ;;
            r) local right_table=$OPTARG ;;
            c) local on_column=$OPTARG ;;
            a) local left_on_column=$OPTARG ;;
            b) local right_on_column=$OPTARG ;;
            t) local target_columns=$OPTARG ;;
            h) local how=$OPTARG ;;
            o) local output_table=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${how} ]]; then local how="left"; fi

    if [[ -z ${on_column} ]]; then
        python3 /paedyl01/disk1/yangyxt/ngs_scripts/pandas_merge_table_api.py \
        -lt ${left_table} \
        -rt ${right_table} \
        -lo ${left_on_column} \
        -ro ${right_on_column} \
        -ot ${output_table} \
        -hw ${how} \
        -tl ${target_columns}
    else
        python3 /paedyl01/disk1/yangyxt/ngs_scripts/pandas_merge_table_api.py \
        -lt ${left_table} \
        -rt ${right_table} \
        -ol ${on_column} \
        -ot ${output_table} \
        -hw ${how} \
        -tl ${target_columns}
    fi
        
    ls -lh ${output_table}
}


function awk_filter_out_columns () {
    #$2 should be target column labels joined by ","
    local input_table=${1}
    local tobe_dropped_cols=${2}

    # to-be dropped columns can be indices or headers. But they must be delimited by comma.
    local -a tobe_dropped_col_arr=($(echo ${tobe_dropped_cols} | awk -F ',' '{for (i=1;i<=NF;i++) printf "%s ", $i;}'))
    local index=1
    for col in "${tobe_dropped_col_arr[@]}"; do
        if [[ ${col} =~ \- ]]; then
            local -a range_arr=($(echo ${col} | awk -F '-' '{printf "%s %s", $1, $2;}'))
            for term in "${range_arr[@]}"; do
                if [[ ${term} =~ ^[0-9]$ ]]; then
                    local index=$(($index * 1))
                elif [[ ${term} =~ [A-Za-z] ]]; then
                    local index=$(($index * 0))
                fi
            done
        elif [[ ${col} =~ ^[0-9]$ ]]; then
            local index=$(($index * 1))
        elif [[ ${col} =~ [A-Za-z] ]]; then
            local index=$(($index * 0))
        fi
    done

    if [[ ${index} -eq 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Input cols are an array of column headers.
        local -a tobe_dropped_col_indices
        for col in "${tobe_dropped_col_arr[@]}"; do 
            local col_ind=$(head -1 ${input_table} | awk -F '\t' -v cl="${col}" '{for (i=1; i<=NF; i++) {if ((NR==1)&&($i == cl)) printf "%s ", i;}}' | awk '{printf "%s", $1;}')
            echo $'\t'"Line "${LINENO}": In function "${FUNCNAME}: In table ${input_table}, the column ${col} index is ${col_ind}.
            if [[ ${#col_ind} -gt 0 ]]; then tobe_dropped_col_indices+=( ${col_ind} ); fi
        done
        local tobe_dropped_cols=$(join_by , "${tobe_dropped_col_indices[@]}")
        time cut --complement -f ${tobe_dropped_cols} < ${input_table} > ${input_table}.tmp
        check_return_code
        check_result_size ${input_table}.tmp && \
        mv ${input_table}.tmp ${input_table}
    elif [[ ${index} -eq 1 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Input cols are an array of column indices.
        time cut --complement -f ${tobe_dropped_cols} < ${input_table} > ${input_table}.tmp
        check_return_code
        check_result_size ${input_table}.tmp && \
        mv ${input_table}.tmp ${input_table}
    fi

    check_result_size ${input_table}
}


function return_array_intersection () {
    # Both arrays should not contain item values with space
    local -a array1=($(echo ${1}))
    local -a array2=($(echo ${2}))
    local -a special_char=("." "+" "?" "^" "$" "(" ")" "[" "]" "{" "}" "|" ":")
    local spe=0

    for item in "${array1[@]}"; do
        local new_item=${item}
        for char in "${special_char[@]}"; do
            if [[ ${item} == *"${char}"* ]]; then local spe=1 && local new_item=$(echo ${new_item} | awk '{gsub("\\'${char}'","\\'${char}'",$0); print;}'); fi
        done
        # if [[ ${spe} -gt 0 ]]; then echo "Line "${LINENO}": In function "${FUNCNAME}: After adding escape symbol to special characters, iterating item ${item} now looks like ${new_item}; fi
        if [[ ${array2[*]} =~ ${new_item} ]]; then
            result+=(${item})
        fi
    done
    echo "${result[*]}"
}


function replace_contig_lines_in_vcfhead () {
    local vcf=${1}
    local ref_gen=${2}

    if [[ -z ${ref_gen} ]]; then
        local ref_gen=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.dict
    elif [[ ${ref_gen} =~ \.fasta$ ]]; then
        ${gatk} CreateSequenceDictionary -R ${ref_gen} -O ${ref_gen/.fasta/.dict} && \
        local ref_gen=${ref_gen/.fasta/.dict}
    fi
    
    # Create a contig header based on the order in ref genome
    awk -F '\t' 'NR>1 {printf "%s:%s\n", $2,$3;}' ${ref_gen} | awk -F ':' '{printf "##contig=<ID=%s,length=%s>\n", $2, $4;}' > ${ref_gen/.dict/.contig_header}

    if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
        gunzip -c -f ${vcf} > ${vcf/.gz/} && \
        local ori_vcf=${vcf}
        local vcf=${vcf/.gz/}
    fi

    awk -F '\t' '{if ($1 ~ /^##contig=/) next; else print}' < ${vcf} > ${vcf}.without_contig
    first_line=$(awk -F '\t' '{if ($1 ~ /^##contig=/) printf NR"\n"; else next;}' < ${vcf} | sort -n - | awk 'NR==1{printf $1;}')
    awk -F '\t' '{if (NR < '${first_line}') print; else exit 0;}' < ${vcf} > ${vcf}.before_contig
    awk -F '\t' '{if (NR >= '${first_line}') print;}' < ${vcf}.without_contig > ${vcf}.after_contig
    cat ${vcf}.before_contig /paedyl01/disk1/yangyxt/indexed_genome/backup_sorted_vcf_contig_header.txt ${vcf}.after_contig > ${vcf}.new && \
    rm ${vcf}.before_contig && rm ${vcf}.after_contig && rm ${vcf}.without_contig && \
    mv ${vcf}.new ${vcf} && \
    if [[ ! -z ${ori_vcf} ]]; then
        bgzip -f ${vcf} && \
        bcftools sort -Oz -o ${ori_vcf} ${ori_vcf} && \
        tabix -f -p vcf ${ori_vcf}
    fi
}


function insert_contig_header_by_genome () {
    local vcf=${1}
    local ref_gen=${2}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    if [[ -z ${ref_gen} ]]; then
        local ref_gen=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.dict
    elif [[ ${ref_gen} =~ \.fasta$ ]]; then
        if [[ -f ${ref_gen/.fasta/.dict} ]]; then
            local ref_gen=${ref_gen/.fasta/.dict}
        else
            ${gatk} CreateSequenceDictionary -R ${ref_gen} -O ${ref_gen/.fasta/.dict} && \
            local ref_gen=${ref_gen/.fasta/.dict}
        fi
    fi
    
    # Create a contig header based on the order in ref genome
    mawk -F '\t' 'NR>1 {printf "%s:%s\n", $2,$3;}' ${ref_gen} | awk -F ':' '{printf "##contig=<ID=%s,length=%s>\n", $2, $4;}' > ${ref_gen/.dict/.contig_header}

    # Insert the contig header part just above INFO header lines
    if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
        zcat ${vcf} | awk -F '\t' '$1 ~ /^##INFO/{exit 0 ;} $1 !~ /^##INFO/{print;}' > ${vcf/.vcf.gz/.tmp.vcf} 2> /dev/null && \
        cat ${ref_gen/.dict/.contig_header} >> ${vcf/.vcf.gz/.tmp.vcf} && \
        local info_row=$(zcat ${vcf} | awk -F '\t' '{if ($1 ~ /^##INFO/) {print NR; exit 0;}}' 2> /dev/null)
        echo ${info_row}
        zcat ${vcf} | tail -n +${info_row} >> ${vcf/.vcf.gz/.tmp.vcf} && \
        bgzip -f ${vcf/.vcf.gz/.tmp.vcf} && bcftools sort -Oz -o ${vcf} ${vcf/.vcf.gz/.tmp.vcf.gz} && tabix -f -p vcf ${vcf}
    elif [[ ${vcf} =~ \.vcf$ ]]; then
        awk -F '\t' '$1 ~ /^##INFO/{exit 0 ;} $1 !~ /^##INFO/{print;}' ${vcf} > ${vcf/.vcf/.tmp.vcf} && \
        cat ${ref_gen/.dict/.contig_header} >> ${vcf/.vcf/.tmp.vcf} && \
        local info_row=$(awk -F '\t' '{if ($1 ~ /^##INFO/) {print NR; exit 0;}}' ${vcf})
        cat ${vcf} | tail -n +${info_row} >> ${vcf/.vcf/.tmp.vcf} && \
        mv ${vcf/.vcf/.tmp.vcf} ${vcf}
    fi
}


function fix_vcf_contig_header_order_by_refgenome () {
    local vcf=${1}
    local ref_gen=${2}

    if [[ -z ${ref_gen} ]]; then
        local ref_gen=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.dict
    elif [[ ${ref_gen} =~ \.dict$ ]]; then
        awk -F '\t' 'NR>1 {print $2;}' ${ref_gen} | awk -F ':' '{printf "%s\n", $2;}' > ${ref_gen}.contig_order
    elif [[ ${ref_gen} =~ \.fasta$ ]]; then
        mawk '$0 ~ /^>/{print $1;}' ${ref_gen} | awk -F '>' '{printf "%s\n", $2;}' > ${ref_gen}.contig_order
    fi

    if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
        zcat ${vcf} | awk -F '\t' 'FNR == NR && $1 ~ /^##contig=<ID=/{con = gensub(/^##contig=<ID=([A-Za-z0-9]+_*[A-Za-z0-9]*_*[A-Za-z0-9]*)[,>].*/, "\\1", "g", $0); con_h_arr[con]=$0;} NR > FNR{print con_h_arr[$0]}' - ${ref_gen}.contig_order | mawk '$0 ~ /^##contig=<ID=/{print;}' > ${vcf::-3}.right_order.contig_lines && \
        cat ${vcf::-3}.right_order.contig_lines && \
        zcat ${vcf} | awk -F '\t' '$0 !~ /^##contig/{print;} $0 ~ /^##contig/{exit 0;}' > ${vcf::-3}.head.before_contig && \
        local lines=$(cat ${vcf::-3}.head.before_contig | wc -l ) && \
        zcat ${vcf} | awk -F '\t' -v row="${lines}" '$0 !~ /^##contig/ && FNR > row{print;}' > ${vcf::-3}.head.after_contig && \
        cat ${vcf::-3}.head.before_contig ${vcf::-3}.right_order.contig_lines ${vcf::-3}.head.after_contig > ${vcf::-3}
        bgzip -f ${vcf::-3} && bcftools sort -Oz -o ${vcf} ${vcf} && tabix -f -p vcf ${vcf}
    elif [[ ${vcf} =~ \.vcf$ ]]; then
        awk -F '\t' 'FNR == NR{if ($1 ~ /^##contig=<ID=/) {con = gensub(/^##contig=<ID=([A-Za-z0-9]+_*[A-Za-z0-9]*_*[A-Za-z0-9]*)[,>].*/, "\\1", "g", $0); con_h_arr[con]=$0;} else next;} NR > FNR{print con_h_arr[$0]}' ${vcf} ${ref_gen}.contig_order | mawk '$0 ~ /^##contig=<ID=/{print;}' > ${vcf::-4}.right_order.contig_lines && \
        head ${vcf::-4}.right_order.contig_lines && \
        awk -F '\t' '$0 !~ /^##contig/{print;} $0 ~ /^##contig/{exit 0;}' ${vcf} > ${vcf::-4}.head.before_contig && \
        local lines=$(cat ${vcf::-4}.head.before_contig | wc -l ) && \
        awk -F '\t' -v row="${lines}" '$0 !~ /^##contig/ && FNR > row{print;}' ${vcf} > ${vcf::-4}.head.after_contig && \
        cat ${vcf::-4}.head.before_contig ${vcf::-4}.right_order.contig_lines ${vcf::-4}.head.after_contig > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf}
    fi
}


function check_sex_from_bam(){
    local input_bam=${1}
    which samtools 2>/dev/null || module load samtools 2> /dev/null
    local xcov=$(echo "scale=4; $(samtools idxstats ${input_bam} | grep "X" | cut -f 3)/$(samtools idxstats ${input_bam} | grep "X" | cut -f 2)" | bc)
    local ycov=$(echo "scale=4; $(samtools idxstats ${input_bam} | grep "Y" | cut -f 3)/$(samtools idxstats ${input_bam} | grep "Y" | cut -f 2)" | bc)
    local rat=$(echo "scale=4; ${xcov}/${ycov}" | bc)
    >&2 echo "X:Y ratio: ${rat}"
    echo ${rat}
    module unload samtools 2> /dev/null || :
}


function return_array_substraction() {
    # Only return the items which array1 has and array2 does not have
    # Cannot return the items which array2 has and array1 does not have
    local -a array1=($(echo ${1}))
    local -a array2=($(echo ${2}))
    # echo array1 is ${array1[@]}
    # echo array2 is ${array2[@]}
    local -a special_char=("." "+" "?" "^" "$" "(" ")" "[" "]" "{" "}" "|" ":")
    local spe=0

    for item in "${array1[@]}"; do
        local new_item=${item}
        for char in "${special_char[@]}"; do
            if [[ ${item} == *"${char}"* ]]; then local spe=1 && local new_item=$(echo ${new_item} | awk '{gsub("\\'${char}'","\\'${char}'",$0); print;}'); fi
        done
        # if [[ ${spe} -gt 0 ]]; then echo "Line "${LINENO}": In function "${FUNCNAME}: After adding escape symbol to special characters, iterating item ${item} now looks like ${new_item}; fi
        if [[ ${array2[*]} =~ ${new_item} ]]; then
            result+=( ${item} )
        else
            subs+=( ${item} )
        fi
    done
    echo "${subs[*]}"
}


function merge_columns () {
    local delimiter=":"
    local header="no"

    local OPTIND d l i o c
    while getopts d::l::i:o::c: args
    do 
        case ${args} in
            d) local delimiter=$OPTARG ;;
            i) local input=$OPTARG ;;
            o) local output=$OPTARG ;;
            l) local merged_label=$OPTARG ;;
            c) local to_be_merged_columns=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ ${#output} -le 0 ]]; then local output=${input}; fi

    local -a col_arr=($(echo ${to_be_merged_columns} | awk -F ',' '{for(i=1;i<=NF;i++) printf $i" ";}'))
    if [[ ${col_arr} =~ [0-9]+ ]]; then
        :
    else
        for col in ${col_arr}; do
            local col_inds=${col_inds}$(head -1 ${input} | awk -F '\t' '{for(i=1;i<=NF;i++) {if($i == "'${col}'") printf i" ";}}')
        done
        local -a col_arr=(${col_inds})
    fi
    
    # Note here that the col_arr we have is 1-indexed, not 0-indexed.
    tmp_bed_right=$(gnws ${input}).$(randomID).bed
    tmp_bed_left=$(gnws ${input}).$(randomID).bed
    
    cut -f $(join_by , "${col_arr[@]}") ${input} | awk -F '\t' '{for(i=1;i<NF;i++) {printf $i"'${delimiter}'";} printf $NF"\n";}' > ${tmp_bed_right}
    cut --complement -f $(join_by , "${col_arr[@]}") ${input} > ${tmp_bed_left}
    
    awk 'BEGIN {FS=OFS="\t";} \
    NR == FNR {aa[NR] = $0; next;} \
    NR > FNR {for(i=1;i<NF;i++) {if (i != '${col_arr}') printf $i"\t"; else printf aa[FNR]"\t"$i"\t";} printf $NF"\n";}' ${tmp_bed_right} ${tmp_bed_left} > ${output} && \
    rm ${tmp_bed_left} && rm ${tmp_bed_right}
    
    echo ${output}
}

function insert_column () {
    local input=${1}
    local col_ind=${2}
    local insert_value=${3}

    if [[ -f ${insert_value} ]]; then
        awk 'BEGIN {FS=OFS="\t";} \
        NR == FNR {aa = $0; next;} \
        {for(i=1;i<NF;i++) {if (i != '${col_ind}') printf $i"\t"; else printf aa"\t"$i"\t";} printf $NF"\n";}' ${insert_value} ${input} > ${input}.tmp && \
        mv ${input}.tmp ${input}
    else
        awk -F '\t' '{for(i=1;i<NF;i++) {if(i != '${col_ind}') printf $i"\t"; else printf "'${insert_value}'\t"$i"\t";} printf $NF"\n";}' < ${input} > ${input}.tmp && \
        mv ${input}.tmp ${input}
    fi   
}

function merge_interval_bed () {
    local bed=${1}
    # The input bed file should follow the standard format.
    local tmp_bed=$(gnws ${bed}).$(randomID).bed

    sort -k1,1 -k2,2n ${bed} | bedtools merge -s -c 4 -o distinct -i - > ${tmp_bed} && \
    mv ${tmp_bed} ${bed}
}

function fetch_interval_bed_using_gene_symbols () {
    local genelist=${1}
    local __output_bed=${2}
    local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed

    bash ${central_scripts}/generate_bed_by_gene_names.sh -g ${genelist} -o ${tmp_bed}
    display_table ${tmp_bed}

    # Do this step is only to consider the fact that Gene symbol and exon are respectively stored in column4 and column 5
    merge_columns -i ${tmp_bed} -c 4,5

    # Do this step is for inserting integer value in column5
    insert_column ${tmp_bed} 5 0

    eval $__output_bed="'${tmp_bed}'"
}


function validify_gz_vcf() {
    local vcf=${1}
    if [[ ${input_vcf} =~ \.gz$ ]]; then 
        if check_gz_file_validity ${input_vcf} && [[ ${input_vcf}.tbi -nt ${input_vcf} ]]; then
            return
        elif check_gz_file_validity ${input_vcf}; then
            bcftools sort -Oz -o ${input_vcf} ${input_vcf} && \
            tabix -f -p vcf ${input_vcf}
            return
        else
            return 1
        fi
    else
        bgzip -f -c ${input_vcf} > ${input_vcf}.gz && \
        bcftools sort -Oz -o ${input_vcf}.gz ${input_vcf}.gz && \
        tabix -f -p vcf ${input_vcf}.gz
        return 0
    fi
}



function remove_sec_align_rec () {
    module load samtools

    local input=${1}
    local tmp_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).sam.head.txt

    samtools view -H ${input} > ${tmp_header}
    samtools view ${input} | awk -F '\t' '{if ($2 < 256) print;}' > ${input}.tmp.sam

    if [[ ${input} =~ \.sam$ ]]; then
        cat ${tmp_header} ${input}.tmp.sam > ${input} && rm ${tmp_header} && rm ${input}.tmp.sam
    elif [[ ${input} =~ \.bam$ ]]; then
        cat ${tmp_header} ${input}.tmp.sam > ${input::-4}.sam && rm ${tmp_header} && rm ${input}.tmp.sam && \
        samtools view -hSb ${input::-4}.sam > ${input} && rm ${input::-4}.sam
    fi
}

function remove_error_line_bed () {
    local bed=$1

    awk -F '\t' '{if ($2 >= $3) next; else print;}' < ${bed} > ${bed}.tmp && \
    mv ${bed}.tmp ${bed}
}

function convert_bam_table () {
    local input=${1}
    local input_dir=$(get_parent_dir ${input})

    if [[ -z ${2} ]]; then
        local output_table=${input_dir}/$(randomID).txt
    else
        local output_table=${2}
    fi

    samtools view ${input} > ${output_table}
    echo ${output_table}
}

function condense_on_col () {
    local input=${1}
    local col_name=${2}
    local tmp_input="${input}.$(randomID)"

    cp ${input} ${tmp_input} && \
    python3 /paedyl01/disk1/yangyxt/ngs_scripts/condense_on_certain_col.py \
    ${tmp_input} ${col_name} ${3} && \
    ls -lh ${tmp_input} && \
    check_result_size ${tmp_input} && \
    drop_redundant_rows ${tmp_input} && \
    mv ${tmp_input} ${input}

    trap "rm ${tmp_input}" ERR
}

function get_RG_SM () {
    local alignment=${1}
    samtools view -H ${alignment} | grep @RG | awk -F '\t' '{for (i=1;i<=NF;i++) if($i ~ /^SM\:/) printf $i;}' | awk -F ':' '{printf $NF;}'
}

function get_RG_ID () {
    local alignment=${1}
    samtools view -H ${alignment} | grep @RG | awk -F '\t' '{for (i=1;i<=NF;i++) if($i ~ /^ID\:/) printf $i;}' | awk -F ':' '{printf $NF;}'
}

function insert_align_records_to_another () {
    local OPTIND i t 
    while getopts i:t: args
    do 
        case ${args} in
            i) local insert_part_alignment=$OPTARG ;;
            t) local tobe_inserted_alignment=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local -a insert_queries=($(samtools view ${insert_part_alignment} | awk -F '\t' '{print $1}' | sort -k 1 | uniq - | awk '{printf $1" ";}')) && \
    samtools view -H ${tobe_inserted_alignment} > ${tobe_inserted_alignment}.tmp && \
    printf "@PG\tID:insert_align_records_to_another\tCL:records inserted from %s\n" "${insert_part_alignment}" >> ${tobe_inserted_alignment}.tmp && \
    samtools view ${tobe_inserted_alignment} | awk -F '\t' -v qs="${insert_queries[*]}" 'BEGIN{split(qs,qarr," "); for(i in qarr) q_arr[qarr[i]]=qarr[i];} {if($1 in q_arr) next; else print;}' >> ${tobe_inserted_alignment}.tmp && \
    samtools view ${insert_part_alignment} >> ${tobe_inserted_alignment}.tmp
    local insert_ID=$(get_RG_ID ${insert_part_alignment})
    local insert_SM=$(get_RG_SM ${insert_part_alignment})
    local tobe_insert_ID=$(get_RG_ID ${tobe_inserted_alignment})
    local tobe_insert_SM=$(get_RG_SM ${tobe_inserted_alignment})

    if [[ ${insert_ID} == ${tobe_insert_ID} ]] && [[ ${insert_SM} == ${tobe_insert_SM} ]]; then
        if [[ ${tobe_inserted_alignment} =~ \.bam$ ]]; then samtools sort -O bam -o ${tobe_inserted_alignment} < ${tobe_inserted_alignment}.tmp && samtools index ${tobe_inserted_alignment}; fi
        if [[ ${tobe_inserted_alignment} =~ \.sam$ ]]; then samtools sort -O sam -o ${tobe_inserted_alignment} < ${tobe_inserted_alignment}.tmp; fi
    fi
}

function identify_max_min_bash_array () {
    # When calling this function, input array must be enclosed in double quotes
    local OPTIND m n a
    while getopts m::n::a: args
    do 
        case ${args} in
            a) local array=$OPTARG ;; 
            m) local __max=$OPTARG ;;
            n) local __min=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local -a array=($(echo ${array}))

    local max=${array[0]}
    local min=${array[0]}

    for i in "${array[@]}"; do
        (( i > max )) && local max=$i
        (( i < min )) && local min=$i
    done

    if [[ ${#__max} -gt 0 ]]; then eval $__max="'${max}'"; fi
    if [[ ${#__min} -gt 0 ]]; then eval $__min="'${min}'"; fi
}

function add_PG_to_bam () {
    local OPTIND a p c 
    while getopts a:p:c:: args
    do 
        case ${args} in
            a) local alignment_file=$OPTARG ;;
            p) local PG_ID=$OPTARG ;;
            c) local PG_CL=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    samtools view -H ${alignment_file} > ${alignment_file}.head
    if [[ ${#PG_CL} -gt 0 ]]; then
        printf "@PG\tID:%s\tCL:%s\n" "${PG_ID}" "${PG_CL}" >> ${alignment_file}.head
    else
        printf "@PG\tID:%s\n" "${PG_ID}" >> ${alignment_file}.head
    fi

    samtools view ${alignment_file} >> ${alignment_file}.head && \
    if [[ ${alignment_file} =~ \.bam$ ]]; then
        samtools view ${alignment_file}.head -b -o ${alignment_file} && \
        samtools index ${alignment_file}
    elif [[ ${alignment_file} =~ \.sam$ ]]; then
        samtools view ${alignment_file}.head -o ${alignment_file}
    fi
}


function merge_logs () {
    local -a logs=($(echo $1 | awk -F ',' '{for(i=1;i<=NF;i++) printf $i" ";}'))
    local merged_log=$2

    for log in "${logs[@]}"; do
        awk '{print;} END{printf "\n\n";}' < ${log} > ${log}.txt && \
        mv ${log}.txt ${log}
    done

    cat "${logs[@]}" > ${merged_log} && \
    rm "${logs[@]}"
}


function gathervcfs () {
    # This function is allowing to gather vcf with different header info.
    # So the first thing to do is to merge the header info.
    if [[ ${1}  =~ \/ ]] && [[ ! ${1} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${1} =~ , ]]; then
        local -a vcfs=($(awk '{printf "%s ", $1;}' < ${1}))
    else
        local -a vcfs=($(echo ${1} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    fi
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: the vcfs are ${vcfs[*]}"
    # Check file existence.
    for vcf in "${vcfs[@]}"; do
        if [[ ! -f ${vcf} ]]; then echo "Line "${LINENO}": In function "${FUNCNAME}": ${vcf} not existed."; fi
    done
    
    local merged_vcf=${2}
    local output_shell_script=${3}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    local common_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf.header
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The temporary temp header file is ${common_header}

    if [[ ${1}  =~ \/ ]] && [[ ! ${1} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${1} =~ , ]]; then
        awk 'NR==1{print;}' < $(head -1 ${1}) > ${common_header}.firstline
    elif [[ $(echo ${1} | awk -F ',' '{printf $1}') =~ \.vcf\.gz$ ]]; then
        zcat $(echo ${1} | awk -F ',' '{printf $1}') | head -1 - > ${common_header}.firstline
    else
        head -1 $(echo ${1} | awk -F ',' '{printf $1}') > ${common_header}.firstline
    fi

    for vcf in "${vcfs[@]}"; do
        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk 'NR > 1{if($1 ~ /^##/) print; else next;}' >> ${common_header}
        else
            awk 'NR > 1{if($1 ~ /^##/) print; else next;}' < ${vcf} >> ${common_header}
        fi
    done

    cat ${common_header} | sort - | uniq - > ${common_header}.tmp
    cat ${common_header}.firstline ${common_header}.tmp > ${common_header} && \
    rm ${common_header}.tmp && \
    rm ${common_header}.firstline
    # Note here now the contig order is not consistent with the one in ucsc.hg19.fasta. We need to resort the order
    replace_contig_lines_in_vcfhead ${common_header}

    # Final check header whether there is a line does not start with ##
    local -a error_vcf_header_lines=($(awk '$0 !~ /^##/{printf "%s ", NR;}' < ${common_header}))
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": There are "${#error_vcf_header_lines[@]}" lines in the common header "${common_header}" that are malformed, and they are at line ${error_vcf_header_lines[*]}"
    awk '$0 ~ /^##/{print;}' < ${common_header} > ${common_header}.tmp && mv ${common_header}.tmp ${common_header}
    local -a error_vcf_header_lines=($(awk '$0 !~ /^##/{printf "%s ", NR;}' < ${common_header}))
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": After fixation, there are "${#error_vcf_header_lines[@]}" lines in the common header "${common_header}" that are malformed, and they are at line ${error_vcf_header_lines[*]}"

    # Check whether the first line in common header defines the VCF format version
    awk '{if($0 ~ /^##fileformat=VCFv4\./) print;}' < ${common_header} | sort - | uniq - > ${common_header}.vcf_format
    awk '{if($0 ~ /^##fileformat=VCFv4\./) next; else print;}' < ${common_header} > ${common_header}.main

    if [[ $(wc -l ${common_header}.vcf_format | awk '{printf $1;}') -eq 0 ]]; then
        awk 'BEGIN {printf "##fileformat=VCFv4.2\n";} {print;}' < ${common_header} > ${common_header}.tmp && \
        mv ${common_header}.tmp ${common_header}
    elif [[ $(wc -l ${common_header}.vcf_format | awk '{printf $1;}') -eq 1 ]]; then
        cat ${common_header}.vcf_format ${common_header}.main > ${common_header}
    else
        head -1 ${common_header}.vcf_format > ${common_header} && \
        cat ${common_header}.main >> ${common_header}
    fi
    
    rm ${common_header}.vcf_format ${common_header}.main

    # Substitude the common header to each vcf in the to_be merged list. Save their original head as a back up.
    local -a tmp_vcfs
    for vcf in "${vcfs[@]}"; do
        local tmp_vcf_body=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf
        echo $'\t'"Now dealing with file ${vcf}"

        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            local headerline=$(zcat ${vcf} | awk -F '\t' '$1 == "#CHROM" {printf NR;}')
            gunzip -f -c ${vcf} > ${vcf::-3}
            awk 'NR >= '${headerline}' {if(($1 !~ /^##/) && ($1 ~ /^#*[Cc]*[Hh]*[Rr]*[A-Z0-9_]+$/)) print;}' < ${vcf::-3} > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            # head -${headerline} ${vcf::-3} | awk 'NR < '${headerline}' {if($1 ~ /^##/) print;}' > ${vcf::-3}.backup.head
            # echo "Check the headerline of vcf file: ${vcf}"
            # head ${vcf::-3}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf::-3}.tmp && mv ${vcf::-3}.tmp ${vcf::-3} && \
            bgzip -f -c ${vcf::-3} > ${vcf::-7}.tmp.vcf.gz && tabix -f -p vcf ${vcf::-7}.tmp.vcf.gz && \
            ( rm ${tmp_vcf_body} ${vcf::-3} ${vcf::-3}.backup.head || echo "Three files: ${tmp_vcf_body} ${vcf::-3} ${vcf::-3}.backup.head already deleted." ) && \
            tmp_vcfs+=( ${vcf::-7}.tmp.vcf.gz )
        else
            local headerline=$(cat ${vcf} | awk -F '\t' '$1 == "#CHROM" {printf NR;}')
            awk 'NR >= '${headerline}' {if(($1 !~ /^##/) && ($1 ~ /^#*[Cc]*[Hh]*[Rr]*[A-Z0-9_]+$/)) print;}' < ${vcf} > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            # head -${headerline} ${vcf} | awk 'NR < '${headerline}' {if($1 ~ /^##/) print;}' > ${vcf}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf::-4}.tmp.vcf && \
            bgzip -f ${vcf::-4}.tmp.vcf && tabix -f -p vcf ${vcf::-4}.tmp.vcf.gz && \
            ( rm ${tmp_vcf_body} ${vcf::-3}.backup.head  || echo "Two files: ${tmp_vcf_body} ${vcf::-3}.backup.head already deleted." ) && \
            tmp_vcfs+=( ${vcf::-4}.tmp.vcf )
        fi
    done

    # Run GATK's SortVCF Tool
    # echo "${tmp_vcfs[*]}" | awk 'BEGIN{printf "\
    # time '${gatk}' --java-options \"-Xmx440G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" SortVcf \\\n";} {for (i=1;i<=NF;i++) {printf "\
    # -I "$i" \\\n";} printf "\
    # -O '${merged_vcf}'";}' | awk '{gsub(/^ {4}/, ""); print;}' > ${output_shell_script} && \
    # bash ${output_shell_script} && \
    # rm ${tmp_vcfs[@]}
    
    # Run GATK's GatherVCF tool.
    echo "${tmp_vcfs[*]}" | awk 'BEGIN{printf "\
    time '${gatk}' --java-options \"-Xmx440G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" GatherVcfs \\\n";} {for (i=1;i<=NF;i++) {printf "\
    -I "$i" \\\n";} printf "\
    -O '${merged_vcf}'";}' | awk '{gsub(/^ {4}/, ""); print;}' > ${output_shell_script} && \
    bash ${output_shell_script} && \
    rm "${tmp_vcfs[@]}"
}

function bcftools_concatvcfs () {
    local merged_vcf=${2}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    if [[ ${1}  =~ \/ ]] && [[ ! ${1} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${1} =~ , ]]; then
        local -a vcfs=($(awk '{printf "%s ", $1;}' < ${1}))
    else
        local -a vcfs=($(echo ${1} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    fi
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: the vcfs are ${vcfs[*]}"
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: The merged vcf is ${merged_vcf}"
    # Check file existence.
    local -a empty_vcfs
    for vcf in "${vcfs[@]}"; do
        if check_vcf_validity ${vcf} 1; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": To be merged ${vcf} is valid."
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": ${vcf} not existed or corrupted." && return 1
        fi
        if [[ ${vcf} =~ \.gz$ ]]; then
            if [[ $(zcat ${vcf} | awk '$1 !~ /^#/{print;}' | wc -l | awk '{print $1;}') -eq 0 ]]; then empty_vcfs+=( ${vcf} ); fi
        elif [[ ${vcf} =~ \.vcf$ ]]; then
            if [[ $(awk '$1 !~ /^#/{print;}' ${vcf} | wc -l | awk '{print $1;}') -eq 0 ]];
                then empty_vcfs+=( ${vcf} )
            else
                >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Since ${vcf} is plain text format and bcftools -f need to use bgzipped vcfs, we compress the vcf and index it with tabix."
                bcftools sort -Oz -o ${vcf}.gz ${vcf} && tabix -f -p vcf ${vcf}.gz && \
                ls -lh ${vcf}.gz
            fi
        fi
    done

    if [[ ${#empty_vcfs[@]} -gt 0 ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): These files are empty, remove them before concating: ${empty_vcfs[*]}"
        for ev in "${empty_vcfs[@]}"; do
            vcfs=( "${vcfs[@]/$ev}" )
        done
    fi
    
    local tmp_file_list=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).lst
    echo "${vcfs[*]}" | awk -F '\t' 'BEGIN{RS=" ";} length($1) > 0{gsub(/\n$/, ""); if ($1 ~ /\.gz$/) printf "%s\n", $1; else printf "%s.gz\n", $1;}' > ${tmp_file_list}
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Here is the temp list file storing the paths of to be concat vcfs:"
    ls -lh ${tmp_file_list}
    cat ${tmp_file_list}

    if [[ ${merged_vcf} =~ \.vcf$ ]]; then
        local plain_merged=${merged_vcf}
        local merged_vcf=${merged_vcf}.gz
    fi

    module load bcftools/1.10.2
    # If using file list for input a list of vcfs, each one of them need to be bgzipped and tabix indexed
    set -x
    bcftools concat -o ${merged_vcf} -a --no-version -Oz -f ${tmp_file_list} && rm -f ${tmp_file_list} && \
    bcftools sort -o ${merged_vcf/.vcf.gz/.sorted.vcf.gz} -Oz ${merged_vcf} && \
    mv ${merged_vcf/.vcf.gz/.sorted.vcf.gz} ${merged_vcf}
    set +x
    module unload bcftools/1.10.2
    
    if [[ ! -z ${plain_merged} ]]; then
        gunzip -c ${merged_vcf} > ${plain_merged}
    else
        tabix -f -p vcf ${merged_vcf}
    fi
}


function mergevcfs () {
    # This function is intended to merge variants info from multiple samples(each vcf contains a separate sample ID)
    # This function is allowing to gather vcf with different header info.
    # So the first thing to do is to merge the header info.
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
    local -a vcfs=($(echo ${1} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: the vcfs are ${vcfs[*]}"
    local merged_vcf=${2}
    local output_shell_script=${3}

    local common_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf.header
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The temporary temp header file is ${common_header}

    if [[ $(echo ${1} | awk -F ',' '{printf $1}') =~ \.vcf\.gz$ ]]; then
        zcat $(echo ${1} | awk -F ',' '{printf $1}') | head -1 > ${common_header}.firstline
    else
        head -1 $(echo ${1} | awk -F ',' '{printf $1}') > ${common_header}.firstline
    fi

    for vcf in "${vcfs[@]}"; do
        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk 'NR > 1{if($1 ~ /^##/) print; else exit 0;}' >> ${common_header}
        else
            awk 'NR > 1{if($1 ~ /^##/) print; else exit 0;}' < ${vcf} >> ${common_header}
        fi
    done

    cat ${common_header} | sort - | uniq - > ${common_header}.tmp
    cat ${common_header}.firstline ${common_header}.tmp > ${common_header} && rm ${common_header}.tmp && rm ${common_header}.firstline
    # Note here now the contig order is not consistent with the one in ucsc.hg19.fasta. We need to resort the order
    replace_contig_lines_in_vcfhead ${common_header}

    # Substitude the common header to each vcf in the to_be merged list. Save their original head as a back up.
    for vcf in "${vcfs[@]}"; do
        local tmp_vcf_body=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf

        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            gunzip -f -c ${vcf} > ${vcf::-3} && \
            awk '{if($1 !~ /^##/) print;}' < ${vcf::-3} > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            awk '{if($1 ~ /^##/) print;}' < ${vcf::-3} > ${vcf::-3}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf::-3}.tmp && \
            mv ${vcf::-3}.tmp ${vcf::-3} && \
            bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf} && rm ${tmp_vcf_body}
        else
            cat ${vcf} | awk '{if($1 !~ /^##/) print;}' > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            cat ${vcf} | awk '{if($1 ~ /^##/) print;}' > ${vcf}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf} && \
            check_vcf_validity ${vcf} && \
            rm ${tmp_vcf_body}
        fi
    done

    # check tabix index
    for vcf in "${vcfs[@]}"; do
        if [[ ${vcf} =~ \.vcf$ ]]; then
            bgzip -c ${vcf} > ${vcf}.gz && tabix -f -p vcf ${vcf}.gz
        elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
            tabix -f -p vcf ${vcf}
        fi
    done

    local -a gz_vcfs
    for vcf in "${vcfs[@]}"; do
        if [[ ${vcf} =~ \.vcf$ ]]; then
            gz_vcfs+=("${vcf}.gz")
        elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
            gz_vcfs+=("${vcf}")
        fi
    done

    # Run vcf-merge from  vcftools.
    if [[ ${merged_vcf} =~ \.gz$ ]]; then
        echo "${gz_vcfs[*]}" | awk 'BEGIN{printf "\
        module load vcftools \n\
        time vcf-merge ";} {for (i=1;i<=NF;i++) {printf \
        $i" ";} printf "\
        | bgzip -c > '${merged_vcf}'\n";}' | awk '{gsub(/^ {4}/, ""); print;} END{printf "\nmodule unload vcftools"}' > ${output_shell_script} && \
        bash ${output_shell_script}
    elif [[ ${merged_vcf} =~ \.vcf$ ]]; then
        echo "${gz_vcfs[*]}" | awk 'BEGIN{printf "\
        module load vcftools \n\
        time vcf-merge ";} {for (i=1;i<=NF;i++) {printf \
        $i" ";} printf "\
        > '${merged_vcf}'\n";}' | awk '{gsub(/^ {4}/, ""); print;} END{printf "\nmodule unload vcftools"}' > ${output_shell_script} && \
        bash ${output_shell_script}
    fi
}

function check_chaos_sorting () {
    local gene_col=$(awk -F '\t' '{if(NR==1) {for(i=1;i<=NF;i++) if($i == "Gene.refGene") printf i}}' < $1)
    local aachange_col=$(awk -F '\t' '{if(NR==1) {for(i=1;i<=NF;i++) if($i == "AAChange.refGene") printf i}}' < $1)
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The column index of Gene.refGene is ${gene_col}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The column index of AAChange.refGene is ${aachange_col}
    local cant_match_rows=$(awk -F '\t' '{if (($'${aachange_col}' == ".")||($'${aachange_col}' == "UNKNOWN")) next; else printf $'${gene_col}'"\t"$'${aachange_col}'"\n";}' < $1 | awk -F '\t' '{regex = "^"$1":"; if (($2 !~ regex)&&($1 !~ /;/)&&($2 !~ /^[0-9]+-[0-9]+:/)) print NR}' | wc -l | awk '{print $1}')
    local total_rows=$(wc -l $1 | awk '{print $1}')
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The number of rows that do not have a correct correspondence between AAChange.refGene and Gene.refGene is ${cant_match_rows}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The number of total rows in the inspecting table is ${total_rows}
    if [[ ${cant_match_rows} -gt 1 ]]; then awk -F '\t' '{if (($'${aachange_col}' == ".")||($'${aachange_col}' == "UNKNOWN")) next; else printf $'${gene_col}'"\t"$'${aachange_col}'"\n";}' < $1 | awk -F '\t' '{regex = "^"$1":"; if (($2 !~ regex)&&($1 !~ /;/)&&($2 !~ /^[0-9]+-[0-9]+:/)) printf NR"\t"$1"\t"$2"\n"}' > ${1::-4}.chaos_order_rows.txt; fi
    local ratio=$(echo "${cant_match_rows}/${total_rows}" | bc -l)
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The ratio is ${ratio}
    local decision=$(awk 'BEGIN{ ratio='${cant_match_rows}'/'${total_rows}'; if ( ratio > 0.1 ) {printf 1;} else printf 0}')
    if [[ ${decision} -eq 1 ]]; then echo "The table $1 order is in a mess. Exiting now" && exit 1; else echo "The table $1 row order is normal. Continue"; fi
}


function filter_on_GT_from_sample () {
    local OPTIND v g s
    while getopts v:g:s: args
    do 
        case ${args} in
            v) local vcf=$OPTARG ;;
            g) local genotype=$OPTARG ;;
            s) local sample=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local genotype=$(echo ${genotype} | awk '{gsub("/", "\\/"); gsub(/\|/, "\\|"); print;}')
    echo The genotype to be filtered out is ${genotype}
    if [[ ${vcf} =~ \.vcf$ ]]; then
        local sample_ind=$(awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}' < ${vcf}) && \
        awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' ~ /^'${genotype}':/) next; else print;}' < ${vcf} > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf}
    elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
        local sample_ind=$(zcat ${vcf} | head -1000 | awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}') && \
        zcat ${vcf} | awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' ~ /^'${genotype}':/) next; else print;}' > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf::-3} && \
        bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf}
    fi
}


function extract_on_GT_from_sample () {
    local OPTIND v g s
    while getopts v:g:s: args
    do 
        case ${args} in
            v) local vcf=$OPTARG ;;
            g) local genotype=$OPTARG ;;
            s) local sample=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local genotype=$(echo ${genotype} | awk '{gsub("/", "\\/"); gsub(/\|/, "\\|"); print;}')
    echo The genotype to be extracted is ${genotype}
    if [[ ${vcf} =~ \.vcf$ ]]; then
        local sample_ind=$(awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}' < ${vcf}) && \
        awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' !~ /^'${genotype}':/) next; else print;}' < ${vcf} > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf}
    elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
        local sample_ind=$(zcat ${vcf} | head -1000 | awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}') && \
        zcat ${vcf} | awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' !~ /^'${genotype}':/) next; else print;}' > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf::-3} && \
        bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf}
    fi
}


function get_intersect_length () {
    local bed_a=${1}
    local bed_b=${2}

    # Noted here bed_b can be a string of multiple bed files path delimited by comma
    module load BEDTools/2.27.1 2> /dev/null
    bedtools intersect -a ${bed_a} -b ${bed_b} | awk -F '\t' '{len=$3-$2; sum += len;} END{printf sum;}'
}

function create_asso_array () {
    local OPTIND t k v a
    while getopts t:k:v:a: args
    do 
        case ${args} in
            t) local table_path=$OPTARG ;;
            k) local key_col_ind=$OPTARG ;;
            v) local val_col_ind=$OPTARG ;;
            a) local __output_array=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bam sample path." ;;
        esac
    done

    # Make sure the key value is uniq.
    local uniq_len=$(awk -F '\t' '{printf $'${key_col_ind}'"\n";}' < ${table_path} | sort - | uniq - | wc -l | awk '{printf $1}')
    local table_len=$(wc -l ${table_path} | awk '{printf $1}')

    if [[ ${uniq_len} -lt ${table_len} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The key value column is not uniquified. Make sure the key values are unique before using this func."
        exit 1
    fi

    awk -F '\t' 'BEGIN{printf "( "} {printf "[\""$'${key_col_ind}'"\"]=\""$'${val_col_ind}'"\" ";} END{printf ")";}' < ${table_path}
}

function common_annovar () {
    local vcf=$1
    local output_dir=$2
    local andir=/paedyl01/disk1/yangyxt/annovar
    local annvar=$andir/annotate_variation.pl
    local tablan=$andir/table_annovar.pl

    module load Perl 

    time perl ${tablan} ${vcf} $andir/humandb/ \
    -buildver hg19 \
    -out ${output_dir} -remove \
    -protocol refGene,genomicSuperDups,gnomad_genome,1000g2015aug_all,exac03,esp6500siv2_all,cg69,dbscsnv11,dbnsfp33a,clinvar_20210501 \
    -operation g,r,f,f,f,f,f,f,f,f -otherinfo -nastring . -vcfinput
    
    check_return_code
    module unload Perl
}


function mergebams () {
    local -a tobe_merged_aligns=($(echo ${1} | awk -F ',' '{printf $1" ";}'))
    local output_align=${2}
    local delete_originals=${3}

    local first_line=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).head
    samtools view -H ${tobe_merged_aligns} | head -1 > ${first_line}
    local SQ_lines=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).sq.head
    samtools view -H ${tobe_merged_aligns} | awk -F '\t' '$1 ~ /^@SQ/{print;}' > ${SQ_lines}
    local other_headers=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).other.head
    if [[ -f ${other_headers} ]]; then rm -f ${other_headers} && touch ${other_headers}; else touch ${other_headers}; fi
    for align in "${tobe_merged_aligns[@]}"; do
        samtools view -H ${align} | awk -F '\t' '$1 !~ /^@[HS]{1}Q/{print;}' >> ${other_headers}
    done
    local common_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).common.head
    cat ${first_line} ${SQ_lines} > ${common_header}
    sort ${other_headers} | uniq - >> ${common_header}

    local total_body=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).sam
    touch ${total_body}
    for align in "${tobe_merged_aligns[@]}"; do
        samtools view ${align} >> ${total_body}
    done
    sort ${total_body} | uniq - > ${total_body}.tmp && mv ${total_body}.tmp ${total_body}

    if [[ ${output_align} =~ \.bam$ ]]; then
        cat ${common_header} ${total_body} > ${output_align}.sam && \
        rm ${first_line} ${SQ_lines} ${other_headers} ${common_header} ${total_body} && \
        samtools sort ${output_align}.sam -T ${output_align::-4} -O bam -o ${output_align} && \
        samtools index ${output_align}
    elif [[ ${output_align} =~ \.sam$ ]]; then
        cat ${common_header} ${total_body} > ${output_align} && \
        rm ${first_line} ${SQ_lines} ${other_headers} ${common_header} ${total_body}
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Output file format not correctly specified: ${output_align}"
    fi

    if [[ ${delete_originals} == "yes" ]] && [[ -f ${output_align} ]]; then
        # Compare size of original align file and merged align file
        if [[ $(ls -l ${output_align} | awk '{gsub(/ +/, " "); print}' | awk '{print $5}') -gt $(ls -l ${tobe_merged_aligns} | awk '{gsub(/ +/, " "); print}' | awk '{print $5}') ]]; then
            rm "${tobe_merged_aligns[@]}"
        fi
    fi  
}

function check_bam_pairing () {
    local bam=${1}

    sambamba index -q ${bam}

    local -a read1s=($(sambamba view -q -F "first_of_pair" -f sam ${bam} | awk -F '\t' '{print $1;}' | sort - | uniq - | awk '{printf $1" ";}'))
    local -a read2s=($(sambamba view -q -F "second_of_pair" -f sam ${bam} | awk -F '\t' '{print $1;}' | sort - | uniq - | awk '{printf $1" ";}'))

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): ${bam} has ${#read1s[@]} first reads
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): ${bam} has ${#read2s[@]} second reads

    local -a only_read1s=($(return_array_substraction "${read1s[*]}" "${read2s[*]}"))
    local -a only_read2s=($(return_array_substraction "${read2s[*]}" "${read1s[*]}"))
    local -a paired_reads=($(return_array_intersection "${read1s[*]}" "${read2s[*]}"))

    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${bam} has ${#only_read1s[*]} qnames with only read1 and they are: ${only_read1s[*]}"
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${bam} has ${#only_read2s[*]} qnames with only read2 and they are: ${only_read2s[*]}"
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): ${bam} has ${#paired_reads[*]} qnames with both read1 alignment and read2 alignment."
}

function mod_chr_prefix_in_bam () {
    # Never use printf $0, always printf "%s", $0 instead so your code doesn't fail if/when your input contains printf formatting characters. Ditto when using any other input values too.
    local input_align=${1}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The input alignment file is ${input_align}
    if [[ -z ${2} ]]; then local operation="add_chr_prefix"; else local operation=${2}; fi

    if [[ ${operation} == 'add_chr_prefix' ]]; then
        samtools view -H ${input_align} | awk 'BEGIN{FS=OFS="\t";} $1 ~ /^@SQ/{replaced=gensub(/SN:([0-9XYMT]{1,2}|[A-Z0-9\.]+)$/, "SN:chr\\1", "g", $2); $2=replaced; print;} $1 !~ /^@SQ/{print;}' > ${input_align}.newhead
        echo Now the new head is ${input_align}.newhead
        ls -lh ${input_align}.newhead
        samtools view ${input_align} | \
        awk 'BEGIN{FS=OFS="\t";} {\
        printf "%s\t%s\tchr%s", $1, $2, $3; \
        for (i=4;i<=NF;i++) {\
            if ($i ~ /^[SX]A:Z:/) {\
                rep=gensub(/^([SX]A):Z:([0-9XYMT]{1,2}|[A-Z0-9\.]+),/, "\\1:Z:chr\\2,", "g", $i); $i=rep;\
                rep=gensub(/;([0-9XYMT]{1,2}|[A-Z0-9\.]+),/, ";chr\\1,", "g", $i); $i=rep;\
                printf "\t%s", $i;} \
            else {\
                printf "\t%s", $i;}} \
        printf "\n";}' > ${input_align}.newbody
    elif [[ ${operation} == 'rm_chr_prefix' ]]; then
        samtools view -H ${input_align} | awk 'BEGIN{FS=OFS="\t";} $1 ~ /^@SQ/{rep=gensub(/^SN:chr([0-9XYMT]{1,2}|[A-Z0-9\.]+)$/, "SN:\\1", "g", $2); $2=rep; print;} $1 !~ /^@SQ/{print;}' > ${input_align}.newhead
        samtools view ${input_align} | \
        awk 'BEGIN{FS=OFS="\t";} {\
        gensub(/^chr([0-9XYMT]+)$/, "\\1", "g", $3); \
        printf "%s\t%s\t%s", $1 ,$2 ,$3 ;\
        for (i=4;i<=NF;i++) {\
            if ($i ~ /^[SX]A:Z:/) {\
                rep=gensub(/^([SX]A):Z:chr([0-9XYMT]{1,2}|[A-Z0-9\.]+),/, "\\1:Z:\\2,", "g", $i); $i=rep;\
                rep=gensub(/;chr([0-9XYMT]{1,2}|[A-Z0-9\.]+),/, ";\\1,", "g", $i); $i=rep;\
                printf "\t%s", $i;} \
            else \
                printf "\t%s", $i;} \
        printf "\n";}' > ${input_align}.newbody
    fi

    if [[ $(basename ${input_align}) =~ \.bam$ ]]; then
        cat ${input_align}.newhead ${input_align}.newbody > ${input_align}.sam && \
        samtools sort -O bam -o ${input_align} ${input_align}.sam && \
        samtools index ${input_align} && \
        rm ${input_align}.sam
    elif [[ $(basename ${input_align}) =~ \.sam$ ]]; then
        mv ${input_align}.newhead ${input_align}
    fi
}

function change_chr_name_in_bam () {
    local input_align=${1}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The input alignment file is ${input_align}
    if [[ -z ${2} ]]; then local tobe_changed="MT"; else local tobe_changed=${2}; fi
    if [[ -z ${3} ]]; then local changed_to="M"; else local changed_to=${3}; fi

    local -a tobe_changed_contigs=($(echo ${tobe_changed} | awk -F ',' '{for(i=1;i<=NF;i++) printf "%s ", $i;}'))
    local -a changed_to_contigs=($(echo ${changed_to} | awk -F ',' '{for(i=1;i<=NF;i++) printf "%s ", $i;}'))

    samtools view -H ${input_align} > ${input_align}.newhead
    for ind in "${!tobe_changed_contigs[@]}"; do
        awk 'BEGIN{FS=OFS="\t";} $1 ~ /^@SQ/{rep=gensub(/:('${tobe_changed_contigs[$ind]}')$/, ":'${changed_to_contigs[$ind]}'", "g", $2); $2=rep; print;} $1 !~ /^@SQ/{print;}' < ${input_align}.newhead > ${input_align}.newhead.tmp && \
        mv ${input_align}.newhead.tmp ${input_align}.newhead
    done

    samtools view ${input_align} > ${input_align}.newbody
    for ind in "${!tobe_changed_contigs[@]}"; do
        awk 'BEGIN{FS=OFS="\t";} { \
        if ($3 ~ /'${tobe_changed_contigs[$ind]}'$/) \
            printf "%s\t%s\t'${changed_to_contigs[$ind]}'", $1, $2; \
        else \
            printf "%s\t%s\t%s", $1, $2, $3; \
        for (i=4;i<=NF;i++) {\
            if ($i ~ /^[SX]A:Z:/) {\
                rep=gensub(/^([SX]A):Z:([CcHhRr]*)'${tobe_changed_contigs[$ind]}',/, "\\1:Z:\\2'${changed_to_contigs[$ind]}',", "g", $i); $i=rep;\
                rep=gensub(/;([CcHhRr]*)'${tobe_changed_contigs[$ind]}',/, ";\\1'${changed_to_contigs[$ind]}',", "g", $i); $i=rep;\
                printf "\t%s", $i;} \
            else \
                printf "\t%s", $i;} \
        printf "\n";}' < ${input_align}.newbody > ${input_align}.newbody.tmp && \
        mv ${input_align}.newbody.tmp ${input_align}.newbody
    done

    cat ${input_align}.newhead ${input_align}.newbody > ${input_align}.new

    if [[ $(basename ${input_align}) =~ \.bam$ ]]; then
        mv ${input_align}.new ${input_align}.sam && \
        samtools sort -O bam -o ${input_align} ${input_align}.sam && \
        samtools index ${input_align} && \
        rm ${input_align}.sam ${input_align}.newbody ${input_align}.newhead
    elif [[ $(basename ${input_align}) =~ \.sam$ ]]; then
        mv ${input_align}.new ${input_align} && \
        rm ${input_align}.newbody ${input_align}.newhead
    fi
}

function upload_gdrive_per_file_or_dir(){
    local OPTIND p i n d e
    while getopts p:i::n::d::e:: args
    do 
        case ${args} in
            p) local path=$OPTARG ;;
            i) local parent_dir_ID=$OPTARG ;;
            n) local __newIDvar=$OPTARG ;;
            d) local directory=$OPTARG ;;
            e) local email_cont_path=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass input_anno_table and var_gene_interaction." ;;
        esac
    done
	
	local name=$(basename ${path})
	local success=1

    if [[ -z ${parent_dir_ID} ]]; then local parent_dir_ID="1TiavZWvojroz7Oh2JibJTjuXJXGLy14B"; fi

    local parent_dir_url="https://drive.google.com/drive/u/0/folders/${parent_dir_ID}"
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Here is the URL for parent folder, ${parent_dir_url}"

	if [[ -d ${path} ]]; then
		local -a existed_IDs=($(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3=="dir")&&($2=="'${name}'")) printf $1" "}'))
		if [[ ${#existed_IDs[@]} -gt 1 ]]; then 
			for existed_ID in "${existed_IDs[@]::${#existed_IDs[@]}-1}"; do 
                gdrive sync upload --keep-local --delete-extraneous --dry-run ${path} ${existed_ID} && \
                gdrive sync upload --keep-local --delete-extraneous ${path} ${existed_ID}
            done
		elif [[ ${#existed_IDs[@]} -eq 0 ]]; then
			while [[ ${success} -gt 0 ]]; do
				gdrive upload -p ${parent_dir_ID} -r ${path} && local success=0 || local success=$((success+1))
				if [[ ${success} -ge 10 ]]; then
                    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Tried 10 times to upload the directory but failed. Check network."
                    return 1
                else
                    continue
                fi
			done
		fi
		
		# Get the accessing URL 
		local new_ID=$(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3=="dir")&&($2=="'${name}'")) printf $1" "}')
		local viewurl=$(gdrive info ${new_ID} | awk -F ': ' '{if ($1=="ViewUrl") print $2}')
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Here is the viewURL for ${name}: ${viewurl}"
        local dir_name=$(basename ${path})
        if [[ -z ${email_cont_path} ]]; then
            awk -v vu="${viewurl}" -v dir="${dir_name}" -v pd="${parent_dir_url}" \
            'BEGIN{printf "Under folder %s:\nFor folder %s:\n viewing URL is %s.\n\n", pd, dir, vu;}'
        else
            awk -v vu="${viewurl}" -v dir="${dir_name}" -v pd="${parent_dir_url}" \
            'BEGIN{printf "Under folder %s:\nFor folder %s:\n viewing URL is %s.\n\n", pd, dir, vu;}' >> ${email_cont_path}
        fi
		if [[ ! -z ${__newIDvar} ]] && [[ ! -z ${new_ID} ]]; then eval $__newIDvar="'${new_ID}'"; fi
	elif [[ -f ${path} ]]; then
		local -a existed_IDs=($(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3!="dir")&&($2=="'${name}'")) printf $1" "}'))
		if [[ ${#existed_IDs[@]} -gt 0 ]]; then 
			for existed_ID in "${existed_IDs[@]}"; do
                >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): File $(basename ${path}) existed, the existing ID is ${existed_ID}. Update the existing one with command:"
                echo $'\t'"python3 ${central_scripts}/gdrive_python_api.py -id ${existed_ID} -pt ${path}"
                python3 ${central_scripts}/gdrive_python_api.py -id ${existed_ID} -pt ${path}
                # gdrive update -p ${parent_dir_ID} --name $(basename ${path}) --no-progress ${existed_ID} ${path}
                >&2 echo "gdrive update return code is $?"
            done
		else
            while [[ ${success} -gt 0 ]]; do
                gdrive upload -p ${parent_dir_ID} ${path} && local success=0 || local success=$((success+1))
                if [[ ${success} -ge 10 ]]; then
                    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Tried 10 times to upload the file but failed. Check network."
                    return 1
                else 
                    continue
                fi
            done
        fi
		local new_ID=$(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3!="dir")&&($2=="'${name}'")) printf $1" "}')
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The new_ID of uploaded file ${path} is ${new_ID}"
		if [[ ! -z ${new_ID} ]]; then
			local viewurl=$(gdrive info ${new_ID} | awk -F ': ' '{if ($1=="ViewUrl") print $2}')
			local downloadurl=$(gdrive info ${new_ID} | awk -F ': ' '{if ($1=="DownloadUrl") print $2}')
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Here is the viewURL for ${name}: ${viewurl}"
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Here is the Downloading URL for ${name}: ${downloadurl}"
		fi
        local file_name=$(basename ${path})
        if [[ -z ${email_cont_path} ]]; then
            awk -v du="${downloadurl}" -v vu="${viewurl}" -v file="${file_name}" -v pd="${parent_dir_url}" \
            'BEGIN{printf "Under folder %s:\nFor file %s:\n Downloading URL is %s.\n Viewing URL is %s.\n\n", pd, file, du, vu;}'
        else
            awk -v du="${downloadurl}" -v vu="${viewurl}" -v file="${file_name}" -v pd="${parent_dir_url}" \
            'BEGIN{printf "Under folder %s:\nFor file %s:\n Downloading URL is %s.\n Viewing URL is %s.\n\n", pd, file, du, vu;}' >> ${email_cont_path}
        fi
		if [[ ! -z ${__newIDvar} ]] && [[ ! -z ${new_ID} ]]; then eval $__newIDvar="'${new_ID}'"; fi
	elif [[ ${directory} == "dir" ]]; then
		local -a existed_IDs=($(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3=="dir")&&($2=="'${name}'")) printf $1" ";}'))
		if [[ ${#existed_IDs[@]} -eq 0 ]]; then
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Trying to establish a new folder named ${path}, NO same name folder established before."
			while [[ ${success} -gt 0 ]]; do
				gdrive mkdir -p ${parent_dir_ID} ${path} && local success=0 || local success=$((success+1))
				if [[ ${success} -ge 10 ]]; then return; fi
			done
		elif [[ ${#existed_IDs[@]} -ge 1 ]]; then
			# for existed_ID in "${existed_IDs[@]::${#existed_IDs[@]}-1}"; do
			# 	gdrive delete -r ${existed_ID}
			# done
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Already has a folder named ${name} in the gdrive."
		fi
		
		local new_ID=$(gdrive list -m 1000 --name-width 400 --order modifiedTime --query "'${parent_dir_ID}' in parents" | awk '{gsub(/ {2,}/, "\t"); print}' | awk -F '\t' '{if(($3=="dir")&&($2=="'${name}'")) {printf "%s",$1; exit 0;}}')
		if [[ ! -z ${__newIDvar} ]] && [[ ! -z ${new_ID} ]]; then eval $__newIDvar="'$new_ID'"; fi
    else
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): To be uploaded file or dir ${path} does not exist."
	fi
}


function pad_bed () {
    local OPTIND b p g o
    while getopts b:p:g::o:: args
    do 
        case ${args} in
            b) local bed_file=$OPTARG ;;
            p) local pad_length=$OPTARG ;;
            g) local genome_size=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bed path and padding length." ;;
        esac
    done
    
    if [[ -z ${genome_size} ]]; then local genome_size=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.contigsize.genome; fi
    if [[ -z ${output_bed} ]]; then local output_bed=${bed_file}; fi
    
    module load BEDTools/2.27.1

    if [[ $(awk -F '\t' 'END{printf NF;}' < ${bed_file}) -ge 4 ]]; then
        mergeBed -i ${bed_file} -c 4 -o distinct -s > ${bed_file}.tmp.bed
    else
        mergeBed -i ${bed_file} -s > ${bed_file}.tmp.bed
    fi

    if [[ ${pad_length} =~ ^0\.[0-9]+$ ]]; then
        slopBed -i ${bed_file}.tmp.bed -g ${genome_size} -b ${pad_length} -pct > ${bed_file}.tmp
    elif [[ ${pad_length} =~ ^[0-9]+$ ]]; then
        slopBed -i ${bed_file}.tmp.bed -g ${genome_size} -b ${pad_length} > ${bed_file}.tmp
    fi

    module unload BEDTools/2.27.1

    mv ${bed_file}.tmp ${output_bed} && rm ${bed_file}.tmp.bed
}

function vcf_errorline_correction () {
    local input_vcf=${1}

    ls -lh ${input_vcf}

    if [[ ${input_vcf} =~ \.vcf\.gz$ ]]; then
        if [[ -f ${input_vcf::-3} ]] && [[ ! -f ${input_vcf} ]]; then
            local plain_vcf=${input_vcf::-3}
        else
            gunzip -c ${input_vcf} > ${input_vcf::-3} || zcat ${input_vcf} > ${input_vcf::-3}
            if [[ $? -eq 0 ]]; then
                local plain_vcf=${input_vcf::-3}
            else
                exit 1
            fi
        fi
    else
        local plain_vcf=${input_vcf}
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Plain text vcf file is ${plain_vcf}."
    ls -lh ${plain_vcf}    

    while true; do
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Plain text vcf file is ${plain_vcf}."
        ls -lh ${plain_vcf}
        if [[ ${#plain_vcf} -eq 0 ]]; then exit 1; fi 
        local -a error_lines=($(awk -F '\t' '$1 !~ /^chr/ && $1 !~ /^#/{printf NR" "; exit 0;}' < ${plain_vcf}))
        if [[ ${#error_lines[@]} -eq 0 ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": No error formated lines detected in this vcf file: ${input_vcf}."
            break
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Line num ${error_lines} is with error format in file ${input_vcf}.
            local err_line_cont=$(awk 'NR=='${error_lines}'{printf $0;}' < ${plain_vcf})
            awk -v mis_part="${err_line_cont}" -v err_line_num="${error_lines}" 'NR == err_line_num - 1 {printf $0""mis_part"\n";} NR != err_line_num && NR != err_line_num - 1{print;}' < ${plain_vcf} > ${plain_vcf}.corrected && \
            mv ${plain_vcf}.corrected ${plain_vcf}
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Line num ${error_lines} and its previous line in ${input_vcf} now looks like: $'\n' $(awk -v err_line_num="${error_lines}" 'NR == err_line_num || NR == err_line_num - 1{print;}' < ${plain_vcf})
        fi
    done

    if [[ ${input_vcf} =~ \.vcf\.gz$ ]] && [[ -f ${plain_vcf} ]]; then
        if [[ $(wc -l ${plain_vcf} | awk '{printf $1;}') -gt 1000 ]]; then
            bgzip -f ${plain_vcf} && tabix -f -p vcf ${input_vcf}
        fi
    fi
}

function get_array_index () {
    # Only applied to non-associative array
    local -a values=($(echo ${1}))
    local -a array=($(echo ${2}))
    local -a indices

    for value in "${values[@]}"; do
        for i in "${!array[@]}"; do
            if [[ "${array[$i]}" == "${value}" ]]; then
                indices+=("${i}")
            fi
        done
    done

    echo "${indices[*]}"
}

function backup_data_to_staging () {
	local usage=$(du -h --max-depth=0 /staging | awk -F '\t' '{prinf $1;}')
    du -h --max-depth=0 /staging 2> /dev/null || local ssd_avail="not_available"
	if [[ ${ssd_avail} == "not_available" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This job is not running on paedyl02 node, abandon migrate the process to staging."
		return 0 
	elif [[ ${usage} =~ ^(5\.[5-9]+T|6\.[0-9]+T) ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": /staging folder is nearly full. Already use ${usage} space."
	else
		# backup probe_dir
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wes/healthy_bams_for_CNV /staging/wes
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wesplus/healthy_bams_for_CNV /staging/wesplus
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wgs/healthy_bams_for_CNV /staging/wgs

        # backup another folder
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wes/backup_gvcfs /staging/wes
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wesplus/backup_gvcfs /staging/wesplus
        rsync --rsync-path=/bin/rsync -avu /paedyl01/disk1/yangyxt/wgs/backup_gvcfs /staging/wgs
	fi
}


function access_full_path_by_sampID () {
    local input_IDs=${1}
    local file_type=${2}
    local __output_var_name=${3}
    local skip=${4}

    if [[ -f ${input_IDs} ]]; then
        local -a IDs=($(awk '{printf $1" ";}' < ${input_IDs}))
    else
        local -a IDs=($(echo ${input_IDs}))
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input sample IDs are ${IDs[*]}."

    local -a seq_types=(wes wgs wesplus)
    local tmp_pair_info=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tmp
    touch ${tmp_pair_info}

    for seq_type in "${seq_types[@]}"; do
        awk -F '\t' -v ids="${IDs[*]}" 'BEGIN{split(ids, idarr, " ");} NR > 1{for(i in idarr) {if($2 == idarr[i]) printf $2"\t"$7"\t'${seq_type}'\n";}}'< /paedyl01/disk1/yangyxt/${seq_type}/${seq_type}_total_ped.ped >> ${tmp_pair_info}
    done

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Temp file storing sample ID and batch info is ${tmp_pair_info}. It has $(wc -l ${tmp_pair_info} | awk '{printf $1;}') rows. Its content are: "
    cat ${tmp_pair_info}

    if [[ ${file_type} == "vcf" ]]; then
        module load vcftools
        module load parallel/20210122
        
        parallel --halt now,fail=1 --dry-run -j3 -k --link cat /paedyl01/disk1/yangyxt/{3}/{2}/vcfs/{2}.1KGPH3.HC.postCGP.Gfiltered.vcf '|' vcf-subset -c {1} '>' /paedyl01/disk1/yangyxt/{3}/{2}/vcfs/{1}.vcf ::: $(awk -F '\t' '{printf "%s ", $1;}' < ${tmp_pair_info}) ::: $(awk -F '\t' '{printf "%s ", $2;}' < ${tmp_pair_info}) ::: $(awk -F '\t' '{printf "%s ", $3;}' < ${tmp_pair_info}) && \
        parallel --halt now,fail=1 --joblog /paedyl01/disk1/yangyxt/tmp_log/access_${file_type}_by_sampID.log -j2 -k --link cat /paedyl01/disk1/yangyxt/{3}/{2}/vcfs/{2}.1KGPH3.HC.postCGP.Gfiltered.vcf '|' vcf-subset -c {1} '>' /paedyl01/disk1/yangyxt/{3}/{2}/vcfs/{1}.vcf ::: $(awk -F '\t' '{printf "%s ", $1;}' < ${tmp_pair_info}) ::: $(awk -F '\t' '{printf "%s ", $2;}' < ${tmp_pair_info}) ::: $(awk -F '\t' '{printf "%s ", $3;}' < ${tmp_pair_info}); \
        cat /paedyl01/disk1/yangyxt/tmp_log/access_${file_type}_by_sampID.log || echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Job log file /paedyl01/disk1/yangyxt/tmp_log/access_${file_type}_by_sampID.log not existed."
        module unload parallel/20210122
        module unload vcftools

        local -a output_paths=($(awk -F '\t' '{printf "/paedyl01/disk1/yangyxt/"$3"/"$2"/vcfs/"$1".vcf ";}' < ${tmp_pair_info}))
    elif [[ ${file_type} == "bam" ]]; then
        local -a output_paths=($(awk -F '\t' '{printf "/paedyl01/disk1/yangyxt/"$3"/"$2"/aligned_results/"$1".bqsr.bam ";}' < ${tmp_pair_info}))
    elif [[ ${file_type} == "fastq" ]]; then
        local -a output_paths=($(awk -F '\t' '{printf "/paedyl01/disk1/yangyxt/"$3"/"$2"/raw_data/"$1"_1.fastq.gz "; \
                                               printf "/paedyl01/disk1/yangyxt/"$3"/"$2"/raw_data/"$1"_2.fastq.gz ";}' < ${tmp_pair_info}))
    fi

    rm -f ${tmp_pair_info}

    # Check output file paths
    for path in "${output_paths[@]}"; do
        if [[ ! -f ${path} ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": ${path} does not exist."
            if [[ ${skip} == "yes" ]]; then
                local -a output_paths=( "${output_paths[@]/$path}" )
            else
                return 1
            fi
        fi
    done

    eval $__output_var_name="'${output_paths[*]}'"
}


function find_involved_vars () {
    local __declare_content__=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tmp
    declare | awk '$0 !~ /^__declare_content__=/{print;}' > ${__declare_content__}
    local -a raw_var_names=($(cat ${__declare_content__} | awk '$0 ~ /^([A-Za-z0-9_]+)=/{vn=gensub(/^([A-Za-z0-9_]+)(=| \(\)).*/, "\\1", "g", $1); printf "%s ",vn;}'))
    local -a var_names
    for var_name in "${raw_var_names[@]}"; do
        if [[ ${var_name} != "_" ]]; then var_names+=( "${var_name}" ); fi
    done
    local -a func_names=($(cat ${__declare_content__} | awk '$0 ~ /^([A-Za-z0-9_]+) \(\)/{vn=gensub(/^([A-Za-z0-9_]+)(=| \(\)).*/, "\\1", "g", $1); printf "%s ",vn;}'))
    # local -a func_names=( "${func_names[@]/${FUNCNAME[0]}}" )

    # Deal with some special cases env_parallel
    local parallel_func="env_parallel"
    local -a filtered_func_names
    for func in "${func_names[@]}"; do
        if [[ ${func} != "env_parallel" ]]; then filtered_func_names+=( "${func}" ); fi
    done
    local -a func_names=(${filtered_func_names[@]})
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Current defined var_names are ${var_names[*]}"
    >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Current defined func_names are ${func_names[*]}"

    local -a target_func_names=(${1})
    local used_function_tmp_storage=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tmp
    local used_variable_tmp_storage=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tmp
    local all_used_funcs=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tmp
    echo "${target_func_names[*]}" | awk 'BEGIN{RS=" ";} { gsub(/\n$/, ""); printf "%s ", $1;}' > ${all_used_funcs}
    touch ${used_function_tmp_storage}
    touch ${used_variable_tmp_storage}
    
    # Recursively find what functions and variables are used within the target functions
    while true; do
        : > ${used_function_tmp_storage}
        : > ${used_variable_tmp_storage}
        for func_name in "${target_func_names[@]}"; do
            # Locate the function body content
            local function_start=$(cat ${__declare_content__} | awk '$0 ~ /^'${func_name}' \(\) /{printf NR}')
            local function_start=$((function_start + 1))
            local function_end=$(cat ${__declare_content__} | awk 'NR > '${function_start}'{if ($0 == "}") printf NR" ";}' | awk '{printf $1;}')

            cat ${__declare_content__} | awk -v fns="${func_names[*]}" 'BEGIN{split(fns, fnarr, " ");} NR > '${function_start}' && NR < '${function_end}' {for (i in fnarr) { exc_reg="[a-zA-Z_]+\.sh "fnarr[i]; inc_reg="[\t\( ]+"fnarr[i]"[ \);]+"; if (($0 ~ inc_reg) && ($0 !~ exc_reg )) printf "%s\n", fnarr[i];}}' > ${used_function_tmp_storage}.tmp 2> /dev/null
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": These functions ( $(awk '{printf "%s ", $1;}' < ${used_function_tmp_storage}.tmp)) are called within this function: ${func_name}."
            cat ${used_function_tmp_storage}.tmp >> ${used_function_tmp_storage}
            # remove function named env_parallel
            awk '$1 != "env_parallel" {print;}' < ${used_function_tmp_storage} | awk -v tfs="${target_func_names[*]}" 'BEGIN{split(tfs, tfsa, " ");} {fetched=1; for (i in tfsa) {if (tfsa[i] != $1) fetched=fetched * 1; else fetched=fetched * 0;} if (fetched == 1) printf "%s\n", $1; else next;}' | sort - | uniq - > ${used_function_tmp_storage}.tmp && mv ${used_function_tmp_storage}.tmp ${used_function_tmp_storage}
            cat ${__declare_content__} | awk -v vns="${var_names[*]}" 'BEGIN{split(vns, vnarr, " ");} NR > '${function_start}' && NR < '${function_end}' {for (i in vnarr) {regex="\\$[\\{]{0,1}"vnarr[i]; if ($0 ~ regex) printf "%s\n", vnarr[i];}}' >> ${used_variable_tmp_storage}
            # remove env variables since we do not need to import them
            awk '$1 !~ /^[A-Z_]+$/{print;}' < ${used_variable_tmp_storage} > ${used_variable_tmp_storage}.tmp && mv ${used_variable_tmp_storage}.tmp ${used_variable_tmp_storage}
        done
        cat ${used_variable_tmp_storage} | sort - | uniq - > ${used_variable_tmp_storage}.tmp && mv ${used_variable_tmp_storage}.tmp ${used_variable_tmp_storage} && \
        cat ${used_function_tmp_storage} | sort - | uniq - > ${used_function_tmp_storage}.tmp && mv ${used_function_tmp_storage}.tmp ${used_function_tmp_storage} && \
        cat ${used_variable_tmp_storage} >> ${all_used_funcs} && \
        cat ${used_function_tmp_storage} >> ${all_used_funcs}
    
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": These variables(funcs not included) ( $(cat ${used_variable_tmp_storage} | sort - | uniq - | awk '{printf "%s ", $1;}')) are called within these functions: ${target_func_names[*]}."
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": These functions ( $(awk '{printf "%s ", $1;}' < ${used_function_tmp_storage})) are called within these functions: ${target_func_names[*]}."
        
        local -a target_func_names=($(cat ${used_function_tmp_storage} | awk -v tfs="${target_func_names[*]}" 'BEGIN{split(tfs, tfsa, " ");} {fetched=1; for (i in tfsa) {if (tfsa[i] != $1) fetched=fetched * 1; else fetched=fetched * 0;} if (fetched == 1) printf "%s\n", $1; else next;}' | sort - | uniq - | awk '{printf "%s ", $1;}'))
        if [[ ${#target_func_names[@]} -eq 0 ]]; then break; fi
    done

    local -a all_var_names=($(cat ${all_used_funcs} | sort - | uniq -))
    # cat ${__declare_content__} | awk -v avn="${all_var_names[*]}" 'BEGIN{split(avn, avs, " ");} {for (i in avs) if $0}'
    rm -f ${all_used_funcs} ${used_function_tmp_storage} ${used_variable_tmp_storage} ${__declare_content__}
    local __output_var=${2}
    if [[ ${#__output_var} -gt 0 ]]; then eval $__output_var="'${all_var_names[*]}'"; else echo "${all_var_names[*]}"; fi
}


function get_fasta_from_bed () {
    local fasta_file=${1}
    local bed_file=${2}
    local output_fasta_file=${3}

    module load BEDTools/2.27.1
    bedtools getfasta -fi ${fasta_file} -bed ${bed_file} -fo ${output_fasta_file} && \
    module unload BEDTools/2.27.1
}


function remove_indel_vcf () {
	module load vcftools
	local vcf=${1}
	local SNV_vcf=${vcf/.vcf/.SNV.vcf}

	if [[ ${vcf} =~ \.vcf\.[b]*gz$ ]]; then
		vcftools --gzvcf ${vcf} --remove-indels --recode --recode-INFO-all --out ${SNV_vcf}
	else
		vcftools --vcf ${vcf} --remove-indels --recode --recode-INFO-all --out ${SNV_vcf}
	fi
	module unload vcftools
}

function preserve_common_var_vcf () {
	#Only preserve common variants
	local andir=/paedyl01/disk1/yangyxt/annovar
	local annvar=$andir/annotate_variation.pl
	local tablan=$andir/table_annovar.pl
	local vcf=${1}
	local __output_var=${2}
	local vcf_basename=$(echo ${vcf} | awk -F '.' '{for(i=1;i<NF;i++) {if ($i != "vcf") printf "%s.", $i;}}' | awk '{gsub(/\.$/, ""); print;}')
	local common_vcf=${vcf/.vcf/.common.vcf}

    >&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp)": Input vcf file is ${vcf}"

	module load Perl
	perl $tablan ${vcf} $andir/humandb/ -buildver hg19 -out ${vcf_basename} -remove -protocol gnomad_genome -operation f -otherinfo -nastring . -vcfinput && \
	ls -lh ${vcf_basename}.hg19_multianno.txt ${vcf_basename}.hg19_multianno.vcf && \
    rm ${vcf_basename}.avinput
	module unload Perl

	local Annotation_table=${vcf_basename}.hg19_multianno.txt
	local common_v_table=${Annotation_table/.txt/.common.txt}
	local common_v_bed=${common_v_table/.txt/.tsv}
	local af_col=$(awk -F '\t' 'FNR==1{for(i=1;i<=NF;i++) if($i ~ /^gnomAD_genome_ALL/) {printf i; exit 0;}}' < ${Annotation_table})

	>&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp): annontation table is ${Annotation_table}
	awk -F '\t' 'NR == 1 {print;} NR > 1 && $'${af_col}' > 0.01 {print;}' < ${Annotation_table} > ${common_v_table}
	awk -F '\t' 'NR > 1 {printf "%s\t%s\t%s\n", $1, $2, $3;}' < ${common_v_table} > ${common_v_bed}
	
	ls -lh ${common_v_bed}
	head ${common_v_bed}

	module load bcftools
	# Be aware the region file cannot have a header
	if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
		bcftools filter ${vcf} -R ${common_v_bed} -Oz -o ${common_vcf}
	elif [[ ${vcf} =~ \.vcf\.bgz$ ]]; then
		bcftools filter ${vcf} -R ${common_v_bed} -Ob -o ${common_vcf}
	else
		bgzip -c -f ${vcf} > ${vcf}.gz && tabix -f -p vcf ${vcf}.gz
		bcftools filter ${vcf}.gz -R ${common_v_bed} -Ov -o ${common_vcf} && \
		rm ${vcf}.gz ${vcf}.gz.tbi
	fi
	module unload bcftools

	rm -f ${common_v_bed}
	rm -f ${common_v_table}
	rm -f ${Annotation_table}
	rm -f ${Annotation_table/.txt/.vcf}
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Output Variable Name is $__output_var

	if [[ ${#__output_var} -gt 0 ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Start eval commands."$'\n'
		# The variable name cannot be the same as local variable name
		eval $__output_var="'$common_vcf'"
	fi
}


function basic_annovar() {
    local andir=/paedyl01/disk1/yangyxt/annovar
	local annvar=$andir/annotate_variation.pl
	local tablan=$andir/table_annovar.pl
    local input_vcf=${1}
    local __output=${2}
    local vcf_basename=${input_vcf/.vcf*/}

    module load Perl
	perl $tablan ${input_vcf} $andir/humandb/ -buildver hg19 -out ${vcf_basename} -remove -protocol refGene,gnomad_exome,gnomad_genome -operation g,f,f -otherinfo -nastring . -vcfinput && \
	ls -lh ${vcf_basename}.hg19_multianno.txt ${vcf_basename}.hg19_multianno.vcf && rm ${vcf_basename}.avinput
	module unload Perl

    if [[ ! -z ${__output} ]]; then
        eval $__output="'${vcf_basename}.hg19_multianno.vcf'"
    fi
}


function basic_vep() {
    local conservation_file=/paedyl01/disk1/yangyxt/public_data/conservation_score/phylocsf_gerp.sql
	local vep_gerp_score=/paedyl01/disk1/yangyxt/public_data/conservation_score/All_hg19_RS.bw
	local human_ancestor_file=/paedyl01/disk1/yangyxt/public_data/human_ancestral_fasta/human_ancestor.fa.gz
    local input_vcf=${1}

    source "$(which conda | awk -F '/' '{for (i=1;i<NF-1;i++) printf "%s/", $i;}')etc/profile.d/conda.sh"
	conda activate vep

    local forks=$(determine_job_num -m 10 -c 1)
    if [[ ${forks} -gt 4 ]]; then
        local forks=4
    fi

    vep -i ${1::-4}.vcf.gz \
	--format vcf \
    --assembly GRCh37 \
    --database \
    --allow_non_variant \
    --dont_skip \
	--fork ${forks} \
    --buffer_size 1000 \
    --fasta /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta  \
    --tab \
	--verbose \
    --force_overwrite \
    --max_af \
    --check_svs \
    -o ${input_vcf/.vcf*/.vep.tsv}

    conda deactivate

    mawk -F '\t' '{if ($1 ~ /^##/) next; else {print;}}' ${input_vcf/.vcf*/.vep.tsv} | \
    mawk -F '\t' 'BEGIN{OFS=FS="\t";} {if (NR == 1) {gsub("#Uploaded_variation", "uniq_ID"); print;} else print;}' > ${input_vcf/.vcf*/.vep.tmp} && \
    mv ${input_vcf/.vcf*/.vep.tmp} ${input_vcf/.vcf*/.vep.tsv}
}


function preserve_common_exon_var_vcf () {
	#Only preserve common variants
	local andir=/paedyl01/disk1/yangyxt/annovar
	local annvar=$andir/annotate_variation.pl
	local tablan=$andir/table_annovar.pl
	local vcf=${1}
	local __output_var=${2}
	local vcf_basename=$(echo ${vcf} | awk -F '.' '{for(i=1;i<NF;i++) {if ($i != "vcf") printf "%s.", $i;}}' | awk '{gsub(/\.$/, ""); print;}')
	local common_vcf=${vcf/.vcf/.common_exon.vcf}

	>&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp)": Input vcf file is ${vcf}"
    
    module load Perl
	perl $tablan ${vcf} $andir/humandb/ -buildver hg19 -out ${vcf_basename} -remove -protocol gnomad_exome -operation f -otherinfo -nastring . -vcfinput && \
	ls -lh ${vcf_basename}.hg19_multianno.txt ${vcf_basename}.hg19_multianno.vcf && rm ${vcf_basename}.avinput
	module unload Perl

	local Annotation_table=${vcf_basename}.hg19_multianno.txt
	local common_v_table=${Annotation_table/.txt/.common_exon.txt}
	local common_v_bed=${common_v_table/.txt/.tsv}
	local af_col=$(awk -F '\t' 'FNR==1{for(i=1;i<=NF;i++) if($i ~ /^gnomAD_exome_ALL/) {printf i; exit 0;}}' < ${Annotation_table})

	>&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp): annontation table is ${Annotation_table}
	awk -F '\t' 'NR == 1 {print;} NR > 1 && $'${af_col}' > 0.01 {print;}' < ${Annotation_table} > ${common_v_table}
	awk -F '\t' 'NR > 1 {printf "%s\t%s\t%s\n", $1, $2, $3;}' < ${common_v_table} > ${common_v_bed}
	
	ls -lh ${common_v_bed}
	head ${common_v_bed}

	module load bcftools
	# Be aware the region file cannot have a header
	if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
		bcftools filter ${vcf} -R ${common_v_bed} -Oz -o ${common_vcf}
	elif [[ ${vcf} =~ \.vcf\.bgz$ ]]; then
		bcftools filter ${vcf} -R ${common_v_bed} -Ob -o ${common_vcf}
	else
		bgzip -c -f ${vcf} > ${vcf}.gz && tabix -f -p vcf ${vcf}.gz
		bcftools filter ${vcf}.gz -R ${common_v_bed} -Ov -o ${common_vcf} && \
		rm ${vcf}.gz ${vcf}.gz.tbi
	fi
	module unload bcftools

	rm -f ${common_v_bed}
	rm -f ${common_v_table}
	rm -f ${Annotation_table}
	rm -f ${Annotation_table/.txt/.vcf}
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Output Variable Name is $__output_var

	if [[ ${#__output_var} -gt 0 ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Start eval commands."$'\n'
		# The variable name cannot be the same as local variable name
		eval $__output_var="'$common_vcf'"
	fi
}


function preserve_exon_var_vcf() {
	#Only preserve exonic variants
	local andir=/paedyl01/disk1/yangyxt/annovar
	local annvar=$andir/annotate_variation.pl
	local tablan=$andir/table_annovar.pl
	local vcf=${1}
	local __output_var=${2}
	local vcf_basename=$(echo ${vcf} | awk -F '.' '{for(i=1;i<NF;i++) {if ($i != "vcf") printf "%s.", $i;}}' | awk '{gsub(/\.$/, ""); print;}')
	local exon_vcf=${vcf/.vcf/.exon.vcf}

	>&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp)": Input vcf file is ${vcf}"
    
    module load Perl
	perl $tablan ${vcf} $andir/humandb/ -buildver hg19 -out ${vcf_basename} -remove -protocol refGene -operation g -otherinfo -nastring . -vcfinput && \
	rm ${vcf_basename}.avinput
	module unload Perl

	local Annotation_table=${vcf_basename}.hg19_multianno.txt
	local exon_v_table=${Annotation_table/.txt/.exon.txt}
	local exon_v_bed=${exon_v_table/.txt/.tsv}
	local func_col=$(awk -F '\t' 'FNR==1{for(i=1;i<=NF;i++) if($i ~ /^Func.refGene/) {printf i; exit 0;}}' < ${Annotation_table})

	>&2 echo "Line "${LINENO}", In function "${FUNCNAME}, $(timestamp): annontation table is ${Annotation_table}
	awk -F '\t' 'NR==1 {print;} NR > 1 && $'${func_col}' ~ /exonic/ {print;}' < ${Annotation_table} > ${exon_v_table}
	awk -F '\t' 'NR > 1 {printf "%s\t%s\t%s\n", $1, $2, $3;}' < ${exon_v_table} > ${exon_v_bed}
	
	ls -lh ${exon_v_bed}
	head ${exon_v_bed}

	module load bcftools
	# Be aware the region file cannot have a header
	if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
		bcftools filter ${vcf} -R ${exon_v_bed} -Oz -o ${exon_vcf}
	elif [[ ${vcf} =~ \.vcf\.bgz$ ]]; then
		bcftools filter ${vcf} -R ${exon_v_bed} -Ob -o ${exon_vcf}
	else
		bgzip -c ${vcf} > ${vcf}.gz && tabix -f -p vcf ${vcf}.gz
		bcftools filter ${vcf}.gz -R ${exon_v_bed} -Ov -o ${exon_vcf} && \
		rm ${vcf}.gz ${vcf}.gz.tbi
	fi
	module unload bcftools

	rm -f ${exon_v_bed}
	rm -f ${exon_v_table}
	rm -f ${Annotation_table}
	rm -f ${Annotation_table/.txt/.vcf}
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Output Variable Name is $__output_var
	if [[ ${#__output_var} -gt 0 ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Start eval commands."$'\n'
		# The variable name cannot be the same as local variable name
		eval $__output_var="'$exon_vcf'"
	fi
}


function split_vcf_by_sample () {
	local input_vcf=${1}
	local vcf_basename=$(basename ${input_vcf} | awk -F '.' '{for(i=1;i<NF;i++) {if ($i != "vcf") printf "%s.", $i;}}' | awk '{gsub(/\.$/, ""); print;}')
	local sample_file=$(echo ${input_vcf} | awk -F '.' '{for(i=1;i<NF;i++) {if ($i != "vcf") printf "%s.", $i;} printf "txt";}')
	local __output_list=${2}

    # trap "rm -f ${sample_file} || echo ${sample_file} already removed" RETURN

	module load bcftools
	
	if [[ ${input_vcf} =~ \.vcf\.[b]*gz$ ]]; then
		local -a samples=($(zcat ${input_vcf} | head -2000 | mawk -F '\t' '$1 == "#CHROM" {for(i=10;i<=NF;i++) printf "%s ", $i;}'))
		zcat ${input_vcf} | head -2000 | mawk -F '\t' '$1 == "#CHROM" {for(i=10;i<=NF;i++) printf "%s\n", $i;}' > ${sample_file}
	elif [[ ${input_vcf} =~ \.vcf$ ]]; then
		local -a samples=($(cat ${input_vcf} | head -2000 | mawk -F '\t' '$1 == "#CHROM" {for(i=10;i<=NF;i++) printf "%s ", $i;}'))
		head -2000 ${input_vcf} | mawk -F '\t' '$1 == "#CHROM" {for(i=10;i<=NF;i++) printf "%s\n", $i;}' > ${sample_file}
	fi

    if [[ ${#samples[@]} -eq 1 ]]; then
        cp -f ${input_vcf} $(dirname ${input_vcf})/${samples}.vcf.gz
        local -a per_sample_vcfs=( $(dirname ${input_vcf})/${samples}.vcf.gz )
        local output_vcf_list=$(join_by , "${per_sample_vcfs[@]}")
    else
        local line_num=$(wc -l ${sample_file} | awk '{printf $1;}')
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": How many samples are there in this vcf ${input_vcf}: ${line_num}."
        if [[ ${line_num} -le 1024 ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Running: bcftools +split -Oz -o $(get_parent_dir ${input_vcf}) -S ${sample_file} ${input_vcf}"
            bcftools +split -Oz -o $(get_parent_dir ${input_vcf}) -S ${sample_file} ${input_vcf} || { \
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Running: bcftools split failed." && return 1; }
        else
            local -a line_brks=($(seq 1 +1000 ${line_num}))
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": "$(timestamp)": line breaks are ${line_brks[*]}"
            if [ ${line_brks[-1]} -lt ${line_num} ]; then line_brks+=( ${line_num} ); fi
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": "$(timestamp)": line breaks are ${line_brks[*]}"
            for i in "${!line_brks[@]}"; do
                local next=$(awk -v ind="${i}" 'BEGIN{printf "%s", ind+1;}')
                local chunk_file=${sample_file::-4}.chunk${next}.txt
                if [[ ${line_brks[-1]} -eq ${line_brks[${i}]} ]]; then break; fi 
                head -n ${line_brks[${next}]} ${sample_file} | tail -n +${line_brks[${i}]} > ${chunk_file}
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Running bcftools +split -Oz -o $(get_parent_dir ${input_vcf}) -S ${chunk_file} ${input_vcf}"
                bcftools +split -Oz -o $(get_parent_dir ${input_vcf}) -S ${chunk_file} ${input_vcf} && rm -f ${chunk_file} || { \
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Running: bcftools split failed at this batch of samples: ${chunk_file}." && return 1; }
            done
        fi

        module unload bcftools

        local -a per_sample_vcfs=($(awk -F '\t' '{printf "%s.vcf.gz ", $1}' < ${sample_file}))
        local output_vcf_list=$(join_by , "${per_sample_vcfs[@]}")
    fi

    if [[ ${#samples[@]} -eq 1 ]]; then
		# rm -f ${sample_file}
		# rm -f "${per_sample_vcfs[@]}"
		if [[ ${#__output_list} -gt 0 ]]; then eval $__output_list="'${input_vcf}'"; fi
	else
		if [[ ${#__output_list} -gt 0 ]]; then eval $__output_list="'${output_vcf_list}'"; fi
		# rm -f ${sample_file}
	fi
}

function identify_all_field_ID () {
    module load bcftools/1.10.2 2>/dev/null

    local input_vcf=${1}
    local field_name=${2}

    if [[ ${input_vcf} =~ \.vcf\.[b]*gz$ ]]; then
        local -a IDs=($(bcftools view -h ${input_vcf} | awk -F '=<' -v fn="##${field_name}" '$1 == fn {split($2, na, ","); for (i in na) {if ($i ~ /^ID=[A-Za-z_]+$/) {print $i;}}' | awk -F '=' '{printf $2" ";}'))
    else
        local -a IDs=($(awk '$0 ~ /^##/{print;}' ${input_vcf} | awk -F '=<' -v fn="##${field_name}" '$1 == fn {split($2, na, ","); for (i in na) {if ($i ~ /^ID=[A-Za-z_]+$/) {print $i;}}' | awk -F '=' '{printf $2" ";}'))
    fi

    echo "${IDs[*]}"

    module unload bcftools/1.10.2 2>/dev/null
}


function remove_samples_from_vcf () {
    module load bcftools/1.10.2
    local input_vcf=${1}
    # Samples to be removed should be seperated by comma Or a path leading to a file containing one sample ID per row.
    if [[ -f ${2} ]]; then
        local rm_samples=$(awk 'BEGIN{printf "^";} {printf "%s,", $1;}' < ${2})
    else
        local rm_samples=$(echo "^"${2})
    fi

    if [[ -z ${3} ]]; then
        local output_vcf=${input_vcf}
    elif [[ ${3} =~ \/.+\. ]]; then
        local output_vcf=${3}
    else
        local base=$(echo ${input_vcf} | awk -F '.' '{for (i=1; i<NF; i++) {if ($i != "vcf") printf "%s.", $i; else exit 0;}}')
        local extension=$(echo ${input_vcf} | awk -F '.' '{for (i=1; i<NF; i++) {if ($i == "vcf") printf "%s.%s", $i, $(i+1);}}')
        local output_vcf=${base::-1}.part.${extension}
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": "$(timestamp)": To be removed samples are "${rm_samples::-1}
    export TMPDIR=/paedyl01/disk1/yangyxt/test_tmp

    set -x
    if [[ ${input_vcf} =~ \.vcf\.gz$ ]]; then
        bcftools view --force-samples -Oz -s ${rm_samples::-1} ${input_vcf} | bcftools sort -o ${input_vcf%.vcf.gz}.tmp.vcf.gz -T ${TMPDIR} && \
        mv ${input_vcf%.vcf.gz}.tmp.vcf.gz ${output_vcf} && \
        tabix -f -p vcf ${output_vcf}
    elif [[ ${input_vcf} =~ \.vcf\.bgz$ ]]; then
        bcftools view --force-samples -Ob -s ${rm_samples::-1} ${input_vcf} | bcftools sort -o ${input_vcf%.vcf.gz}.tmp.vcf.gz -T ${TMPDIR} && \
        mv ${input_vcf%.vcf.gz}.tmp.vcf.gz ${output_vcf} && \
        tabix -f -p vcf ${output_vcf}
    else
        bcftools view --force-samples -Ov -s ${rm_samples::-1} ${input_vcf} | bcftools sort -o ${input_vcf%.vcf}.tmp.vcf -T ${TMPDIR} && \
        mv ${input_vcf%.vcf}.tmp.vcf ${output_vcf}
    fi
    set +x

    module unload bcftools/1.10.2
}

function split_vcf_by_batch () {
    local wkd=${2}
    local sample_batch=$(basename ${wkd})
    local seq_type=$(echo ${wkd} | awk -F '/' '{printf $(NF-1);}')
    local ped=${wkd}/pedigree
    local bgvcf=/paedyl01/disk1/yangyxt/${seq_type}/backup_gvcfs
    local vcf_dir=${wkd}/vcfs
    local ped_file_name=${3}

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): Starting splitting vcf files for this batch ${sample_batch}"
    local sample_batch_IDs=$(awk -F '\t' '{if(NR > 1) printf "%s,",$2}' < $ped/$ped_file_name.ped)
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): The batch IDs are ${sample_batch_IDs}"

    # Cut off the last character
    local sample_ID_list=${sample_batch_IDs::-1}
    # Extract sample_batch_vcf to the vcf directory in sample_batch
    module load vcftools

    local vcf=${1}
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): The tobe splitted vcf is ${vcf}"

    if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
        time zcat ${vcf} | vcf-subset -c ${sample_ID_list} > ${vcf_dir}/${sample_batch}.HC.VQSR.vcf && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.HC.*.gz* && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.HC.*.recal* && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.g.vcf.gz && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.g.vcf.gz.tbi
    elif [[ ${vcf} =~ \.vcf$ ]]; then
        time cat ${vcf} | vcf-subset -c ${sample_ID_list} > ${vcf_dir}/${sample_batch}.HC.VQSR.vcf && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.HC.* && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.HC.*.recal* && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.g.vcf.gz && \
        rm -f $bgvcf/all_${seq_type}_samples_plus_${sample_batch}.g.vcf.gz.tbi
    fi

    module unload vcftools

    # Compress the vcf file and index it.
    bgzip -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf && \
    tabix -f -p vcf ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Finishing splitting vcf files for this batch"$'\n\n'
}


function extract_vcf_by_samples () {
    local vcf_file=${1}
    local samples=${2}  # Input samples should be delimited by comma or a file path containing one sample ID per line
    local subset_vcf=${3}

    if [[ ${vcf_file} =~ \.vcf$ ]]; then
        bgzip -f -c ${vcf_file} > ${vcf_file}.gz && \
        tabix -f -p vcf ${vcf_file}.gz && \
        local vcf_file=${vcf_file}.gz
    fi

    if [[ ${subset_vcf} =~ \.vcf$ ]]; then
        local plain_subset=${subset_vcf}
        local subset_vcf=${subset_vcf}.gz
    fi
    
    if [[ ! ${samples} =~ , ]] && [[ -f ${samples} ]]; then
        local tmp_samples_file=${samples}
    else
        local tmp_samples_file=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).lst
        touch ${tmp_samples_file}
        echo ${samples} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' > ${tmp_samples_file}
    fi

    module load bcftools
    bcftools view -O z --no-version -S ${tmp_samples_file} -o ${subset_vcf::-7}.tmp.vcf.gz --no-update ${vcf_file}
    if [[ ${tmp_samples_file} =~ ^\/paedyl01\/disk1\/yangyxt\/test_tmp\/.+\.lst$ ]]; then rm -f ${tmp_samples_file}; fi
    module unload bcftools

    if [[ ! -z ${plain_subset} ]]; then
        zcat ${subset_vcf::-7}.tmp.vcf.gz > ${plain_subset} && \
        rm -f ${subset_vcf::-7}.tmp.vcf.gz
    else
        mv ${subset_vcf::-7}.tmp.vcf.gz ${subset_vcf} && tabix -f -p vcf ${subset_vcf}
    fi
}


function batch_vcf_qc () {
    local batch_vcf=${1}
    local wkd=${2}
    local ped_file_name=${3}
    local sample_batch=$(basename ${wkd})
    local ped=${wkd}/pedigree
    local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts

    ls -lh ${batch_vcf}
    local new_folderID
    local sub_folderID

    local expected_samples=$(tail -n +2 ${ped}/${ped_file_name}.ped | awk -F '\t' '{printf "%s,", $2}' | awk '{gsub(/,$/, ""); print;}')

    if check_vcf_validity ${batch_vcf} 1000 ${expected_samples}; then 
        # Create a folder to store the batch_info
        upload_gdrive_per_file_or_dir \
        -p ${sample_batch} \
        -i 1TiavZWvojroz7Oh2JibJTjuXJXGLy14B \
        -n new_folderID \
        -d dir

        # Create a subfolder to store the vcf qc info
        upload_gdrive_per_file_or_dir \
        -p ${sample_batch}_vcf_QC \
        -i ${new_folderID} \
        -n sub_folderID \
        -d dir

        # Do the QC the pass the sub_folderID
        bash ${central_scripts}/vcf_QC.sh ${batch_vcf} ${ped}/${ped_file_name}.ped ${sub_folderID}

        # Check heterozygosity
        bash ${central_scripts}/heterozygosity_check.sh -v ${batch_vcf}
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Specified batch vcf file ${batch_vcf} not existed. Quit now."
        return 1
    fi
}


function cal_posterior_and_filter () {
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
    local ref_gen=/paedyl01/disk1/yangyxt/indexed_genome
    local wkd=${1}
    local ped_file_name=${2}
    local vcf_dir=${wkd}/vcfs
    local ped=${wkd}/pedigree
    local sample_batch=$(basename ${wkd})

    if [[ -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz ]] && \
    [[ -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz.tbi ]] && \
    [[ ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz.tbi -nt ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz ]] && \
    check_gz_file_validity ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz; then
        :
    elif [[ -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz ]] && check_gz_file_validity ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input file is compressed but not indexed or index file is older than gz file. Reindexing"
        tabix -f -p vcf ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz
    elif [[ -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input file not in gz format, compress it and index it."
        bgzip -f ${vcf_dir}/${sample_batch}.HC.VQSR.vcf && tabix -f -p vcf ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Starting CalculatingGenotypePosteriors"
    time $gatk CalculateGenotypePosteriors \
    -R $ref_gen/ucsc.hg19.fasta \
    -supporting $ref_gen/with_chr_1000G_phase3_v4_20130502.sites.vcf.gz \
    -ped $ped/$ped_file_name.ped \
    -V ${vcf_dir}/${sample_batch}.HC.VQSR.vcf.gz \
    -O ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.vcf && \
    echo "** Calculating posterior probabilities done. **" && \
    bgzip -f ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.vcf && \
    tabix -f -p vcf ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.vcf.gz

    check_return_code

    # Filter low quality genotypes
    time $gatk VariantFiltration \
    -R $ref_gen/ucsc.hg19.fasta \
    -V ${vcf_dir}/$sample_batch.1KGPH3.HC.postCGP.vcf.gz \
    -G-filter "GQ < 20.0" \
    -G-filter-name lowGQ \
    -O ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf && \
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "** Filtering low GQ variants DONE **" && \
    bgzip -f ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf && \
    tabix -f -p vcf ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz

    check_return_code

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Finishing CalculatingGenotypePosteriors$'\n\n'
}


function split_vcf_by_fam () {
    local wkd=${1}
    local sample_batch=$(basename ${wkd})
    local vcf_dir=${wkd}/vcfs
    local ped=${wkd}/pedigree
    local gvcf=${wkd}/gvcfs
    local scrf=${wkd}/respective_scripts
    local batchf=${wkd}/job_qsub_scripts

    if [[ -f ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf ]]; then
        bgzip -f -c ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf > ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz && \
        tabix -f -p vcf ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz
    elif [[ -f ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz ]] && check_gz_file_validity ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz; then
        zcat ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf.gz > ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf
    fi
    
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Starting Splitting vcf files by family
    cat $ped/python_modified.ped | awk -F '\t' \
    -v sa=$sample_batch -v wkd=${wkd} \
    -v v=${vcf_dir} \
    -v p=$ped -v g="$gvcf" \
    -v s="$scrf" \
    'BEGIN{printf "    wkd=%s\n\
    vcf_dir=%s\n\
    ped=%s\n\
    vcf_version=1KGPH3.HC.postCGP.Gfiltered\n\
    module load vcftools\n",wkd,v,p} \
    {printf "    cat ${vcf_dir}/%s.\"$vcf_version\".vcf | vcf-subset -c %s | bgzip -f -c > \"$vcf_dir\"/%s.vcf.gz && tabix -f -p vcf \"$vcf_dir\"/%s.vcf.gz\n", sa, $2, $1, $1;} \
    END{printf "\nmodule unload vcftools\n";}' | awk '{gsub(/^\t/, ""); gsub(/^ {3,4}/, ""); print}' > ${batchf}/autogen_extracting_vcf_script.sh && \
    cat ${batchf}/autogen_extracting_vcf_script.sh && \
    time bash ${batchf}/autogen_extracting_vcf_script.sh && \
    rm -f ${vcf_dir}/${sample_batch}.1KGPH3.HC.postCGP.Gfiltered.vcf
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $'\t\t'******Splitting vcf done********
    check_return_code
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): Finish Splitting vcf files by family$'\n\n'
}


function check_mend_error_and_left_align () {
    local vcf_dir=${1}
    local ped=${2}
    local ped_file_name=${3}
    local PICARD=/home/yangyxt/software/picard/picard.jar
    local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts

    cd ${vcf_dir}
    local -a finished_fam_array=($(ls *.vcf.gz | awk -F '.' '{printf "%s\n", $1;}' | sort - | uniq - | awk '{printf "%s ", $1;}'))
    local -a ideal_fam_array=($(cat ${ped}/${ped_file_name}.ped | awk -F '\t' 'NR > 1{printf "%s\n", $1;}' | sort - | uniq - | awk '{printf "%s ", $1;}'))
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The fam array for this batch is ${ideal_fam_array[*]}"
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The fam array for this batch that actually finished calling vcfs is ${finished_fam_array[*]}"
    local -a fam_array=($(return_array_intersection "${finished_fam_array[*]}" "${ideal_fam_array[*]}"))

    for fam in "${fam_array[@]}"
    do
        fam_mem_no=$(awk -F '\t' -v f=${fam} '{if($1 == f) print $2}' < $ped/python_modified.ped | awk -F ',' '{printf NF}')
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": "$'\t\t'"There are ${fam_mem_no} fam members in ${fam}"
        if [[ ${fam_mem_no} -ge 3 ]]
        then
            awk -F '\t' '{if ($1 == "'${fam}'") print}' < ${ped}/${ped_file_name}.ped > ${ped}/${ped_file_name}.${fam}.ped
            check_return_code

            check_vcf_validity ${vcf_dir}/${fam}.vcf.gz

            time java -jar ${PICARD} FindMendelianViolations \
            INPUT=${vcf_dir}/${fam}.vcf.gz \
            TRIOS=${ped}/${ped_file_name}.${fam}.ped \
            OUTPUT=${vcf_dir}/${fam}.mendelian_violation_metrics \
            MIN_DP=4 || (tabix -f -p vcf ${vcf_dir}/${fam}.vcf.gz && time java -jar ${PICARD} FindMendelianViolations \
            INPUT=${vcf_dir}/${fam}.vcf.gz \
            TRIOS=${ped}/${ped_file_name}.${fam}.ped \
            OUTPUT=${vcf_dir}/${fam}.mendelian_violation_metrics \
            MIN_DP=4)
            check_return_code
        fi
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": $(timestamp) Left_Align Indels and Split Multi allelic sites for ${fam}.vcf.gz by running deal_with_multi_allelic_records_at_vcf_level.sh"
        bash ${central_scripts}/deal_with_multi_allelic_records_at_vcf_level.sh ${vcf_dir}/${fam}
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp) Finish Left aligning Indels and Splitting Multi allelic sites.
    done
}

function independent_bwa_alignment () {
	local OPTIND f r s g o p t c e
    while getopts f::r::p::s::g::o::t::c::e:: args
    do
        case ${args} in
            f) local forward_reads=$OPTARG ;;
            r) local reverse_reads=$OPTARG ;;
            c) local cover=$OPTARG ;;
            p) local merged_pair_reads=$OPTARG ;;
            g) local ref_genome=$OPTARG ;;
            s) local samp_ID=$OPTARG ;;
            e) local expected_lines=$OPTARG ;;
            o) local output_align=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${samp_ID} ]]; then local samp_ID=$(basename ${forward_reads} | awk '{gsub(/_[a-z]*1\.f[ast]*q[\.gz]*$/, "", $0); printf "%s", $0;}'); fi
    if [[ -z ${output_align} ]]; then local output_align=${samp_ID}.bam; fi
    if [[ -z ${cover} ]]; then local cover="Yes"; fi
    if [[ -z ${expected_lines} ]]; then local expected_lines=10; fi

    rm -f ${output_align}.tmp.*.bam 2> /dev/null && \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are some temporary files left in last alignment, remove them before alignment " || \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are no temporary tmp.00x.bam file left in last alignment. "

    local cpu=$(get_pbs_cpu)
    local mem=$(get_pbs_mem)
    if [[ ! -z ${mem} ]]; then echo "In function ${FUNCNAME}: $(timestamp): Under current PBS job ${PBS_JOBID}, the total cpu is ${cpu}, and the total mem is ${mem}Gb."; fi
    if [[ -z ${cpu} ]]; then local cpu=8; fi
    if [[ -z ${threads} ]]; then local threads=${cpu}; fi
    if [[ ${threads} -gt 10 ]]; then local threads=10; fi
    
    if check_bam_validity ${output_align} ${expected_lines} && [[ ${cover} == "No" ]]; then
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: Output alignment file ${output_align} already generated. And it's big enough."
    else
        module load bwa
        module load samtools
        
        set -x
        /usr/bin/time bwa mem -t ${threads} -M -R "@RG\tID:${samp_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}" ${ref_genome} \
        ${forward_reads} ${reverse_reads} | samtools view -uSh -@ ${threads} - | samtools sort -O bam -@ ${threads} -o ${output_align} - && \
        samtools index ${output_align} && echo "BWA-MEM alignment done"
        set +x

        module unload samtools
        module unload bwa
    fi
}

function independent_bwamem2_alignment () {
	local OPTIND f r s g o p t
    while getopts f::r::p::s::g::o::t:: args
    do
        case ${args} in
            f) local forward_reads=$OPTARG ;;
            r) local reverse_reads=$OPTARG ;;
            p) local merged_pair_reads=$OPTARG ;;
            g) local ref_genome=$OPTARG ;;
            s) local samp_ID=$OPTARG ;;
            o) local output_align=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${samp_ID} ]]; then local samp_ID=$(basename ${forward_reads} | awk '{gsub(/_[a-z]*1\.f[ast]*q[\.gz]*$/, "", $0); printf "%s", $0;}'); fi
    if [[ -z ${output_align} ]]; then local output_align=${samp_ID}.bam; fi
    if [[ -z ${threads} ]]; then local threads=10; fi
    if [[ ${threads} -gt 10 ]]; then local threads=10; fi

    local bwa2=/home/yangyxt/software/bwa-mem2-2.1_x64-linux/bwa-mem2

    rm -f ${output_align}.tmp.*.bam 2> /dev/null && \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are some temporary files left in last alignment, remove them before alignment " || \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are no temporary tmp.00x.bam file left in last alignment. "
    
    module load samtools

    set -x
    /usr/bin/time ${bwa2} mem -t ${threads} -M -R "@RG\tID:${samp_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}" ${ref_genome} \
    ${forward_reads} ${reverse_reads} | samtools view -uSh -@ ${threads} - | samtools sort -O bam -@ ${threads} -o ${output_align} - && \
    samtools index ${output_align} && echo "BWA-MEM alignment done"
    set +x

    module unload samtools
}

function migrate_to_staging_temp () {
	local to_be_copied=${1}
	local usage=$(du -h --max-depth=0 /staging | awk -F '\t' '{printf $1;}')
	du -h --max-depth=0 /staging 2> /dev/null || local ssd_avail="not_available"

	if [[ ${ssd_avail} == "not_available" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This job is not running on paedyl02 node, abandon migrate the process to staging."
		ls -lh ${to_be_copied}
		return 0
	else
        >&2 echo "SSD available. Migrating Data Now."
		local parent_dir=$(get_parent_dir ${to_be_copied})
		local file_name=$(basename ${to_be_copied})
        set -x
		rsync --rsync-path=/bin/rsync -avu --delay-updates --include="${file_name}" --exclude="*" ${parent_dir}/ /staging/test_tmp/ && ls -lh /staging/test_tmp/${file_name} || (rm -f /staging/test_tmp/${file_name} 2> /dev/null; ls -lh ${to_be_copied})
        set +x
	fi
}

function migrate_from_staging_temp () {
	local mov_file=${1}
	local mov_filename=$(basename ${mov_file})
	local dest_dir=${2}
	
	local usage=$(du -h --max-depth=0 /staging | awk -F '\t' '{printf $1;}')
    du -h --max-depth=0 /staging 2> /dev/null || local ssd_avail="not_available"

    if [[ ${ssd_avail} == "not_available" ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This job is not running on paedyl02 node, abandon migrate the process to staging."
        return 0
    else
        set -x
        rsync --rsync-path=/bin/rsync --dry-run --remove-source-files -avu --delay-updates --include="${mov_filename}" --exclude="*" /staging/test_tmp/ ${dest_dir}
        rsync --rsync-path=/bin/rsync --remove-source-files -avu --delay-updates --include="${mov_filename}" --exclude="*" /staging/test_tmp/ ${dest_dir} || echo "We dont have ${mov_filename} under /staging/test_tmp"
        set +x
    fi
}

function select_vcf_by_interval () {
    local gz_vcf=${1}
    local region_tsv=${2}
    local output_gzvcf=${3}

    # Prepare index if there aren't any
    if [[ ! -f ${gz_vcf}.tbi ]]; then tabix -f -p vcf ${gz_vcf}; fi 

    # In bcftools only TSV file can be seen as 1-based index and inclusive on interval coordinates
    module load bcftools
    bcftools view -R ${region_tsv} -Oz -o ${output_gzvcf} ${gz_vcf}
    module unload bcftools

    tabix -f -p vcf ${output_gzvcf}
    ls -lh ${output_gzvcf}
}

function bcftools_merge_vcfs () {
    local OPTIND i m f c o t s e r
    while getopts i:m::f::c::o:t::s::e::r:: args
    do
        case ${args} in
            i) local input_files=$OPTARG ;;
            m) local missing=$OPTARG ;;
            f) local force_sample=$OPTARG ;;
            c) local collapse=$OPTARG ;;
            e) local merge=$OPTARG;;
            r) local info_rules=$OPTARG;;
            o) local output_vcf=$OPTARG ;;
            t) local threads=$OPTARG ;;
            s) local skip=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    if [[ -f ${input_files} ]]; then
        # Input argument maybe vcf path map (strictly according to the format shown here: http://samtools.github.io/bcftools/bcftools.html#merge)
        # One file path per line
        :
    else
        # Input argument maybe multiple vcf file paths separated by comma
        local tmp_file_lst=$TMPDIR/$(randomID).lst
        echo ${input_files} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' > ${tmp_file_lst}
        local input_files=${tmp_file_lst}
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Now input list file is ${tmp_file_lst}"
    fi

    if [[ -z ${collapse} ]]; then
        local collapse_arg=""
    else
        local collapse_arg="-c ${collapse} "
    fi

    if [[ -z ${info_rules} ]]; then
        local info_rules_arg=""
    else
        local info_rules_arg="-i ${info_rules} "
    fi

    if [[ -z ${merge} ]]; then
        local merge_arg="-m none"
    else
        local merge_arg="-m ${merge} "
    fi
    
    if [[ -z ${force_sample} ]]; then
        local fs_arg=""
    else
        local fs_arg="--force-samples "
    fi

    if [[ -z ${missing} ]]; then
        local m_arg=""
    else
        local m_arg="--missing-to-ref "
    fi

    if [[ -z ${threads} ]]; then
        local t_arg=""
    else
        local t_arg="--threads ${threads} "
    fi

    local -a input_file_arr=($(awk '{printf "%s ", $1;}' < ${input_files}))

    trap "rm ${tmp_file_lst}" RETURN

    if [[ ${skip} == "skip_check_validity" ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Function ${FUNCNAME} was called with argument specified to skip checking input vcfs validity."
    elif [[ ${skip} ]]; then
        for file in "${input_file_arr[@]}"; do
            if [[ -f ${file} ]] && check_vcf_validity ${file} 2> /dev/null; then
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": input vcf ${file} is valid"
            else
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": input vcf ${file} is not valid compressed. drop this one."
                awk -v f="${file}" '$0 != f{print;}' < ${input_files} > ${input_files}.tmp && \
                mv ${input_files}.tmp ${input_files}
            fi
        done
    else
        for file in "${input_file_arr[@]}"; do
            if [[ -f ${file} ]] && check_vcf_validity ${file} 2> /dev/null; then
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": input vcf ${file} is valid"
            else
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": input vcf ${file} is not valid compressed. Since skip is not specified. Quit merging."
                return 1
            fi
        done
    fi


    module load bcftools

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Running this command: bcftools merge ${t_arg}--no-version -l ${input_files} ${fs_arg}-Oz ${m_arg}-o ${output_vcf}"
    bcftools merge ${info_rules_arg} ${t_arg}--no-version -l ${input_files} ${fs_arg}-Oz ${m_arg}-o ${output_vcf} ${merge_arg} && \
    bcftools sort -Oz -o ${output_vcf} ${output_vcf} && \
    bcftools index ${output_vcf} && \
    ls -lh ${output_vcf}

    if [[ ! -z ${tmp_file_lst} ]]; then rm -f ${tmp_file_lst}; fi
    module unload bcftools
}

function convert_to_superpop () {
    local report=${1}
    local superpop_map=${2}
    local pop_super_map=/paedyl01/disk1/yangyxt/public_data/pop_superpop_map.tsv

    awk_merge_two_tables ${report} ${pop_super_map} 2 1 ${superpop_map} 
     
    tail -n +2 ${superpop_map} | cut --complement -f 5,8-10 - > ${superpop_map}.tmp && \
    mv ${superpop_map}.tmp ${superpop_map}
}


function select_qname_to_paired_fastq_from_bam(){
	local bam=${1}
	local qnames_lst=${2}
	local fastq_1=${3}
	local fastq_2=${4}
	local tmp_bam=${bam::-4}.$(randomID).bam
    local samtools="/home/yangyxt/software/samtools-1.14/samtools"
	
	module load BEDTools

    samtools sort -n -O bam -o ${tmp_bam::-4}.tmp.bam ${bam} && \
	${samtools} view -N ${qnames_lst} -o ${tmp_bam} ${tmp_bam::-4}.tmp.bam && \
	>&2 echo "$(timestamp): In function ${FUNCNAME}: Use qname list ${qnames_lst} to extract from bam file ${bam} to temp bam ${tmp_bam}, with ${samtools}"
	ls -lh ${tmp_bam}

	if [[ ${fastq_1} =~ \.gz$ ]]; then
		bedtools bamtofastq -i ${tmp_bam} -fq ${fastq_1::-3} -fq2 ${fastq_2::-3} && \
		gzip -f ${fastq_1::-3} && \
		gzip -f ${fastq_2::-3}
	else
		bedtools bamtofastq -i ${tmp_bam} -fq ${fastq_1} -fq2 ${fastq_2}
	fi

	>&2 echo "$(timestamp): In function ${FUNCNAME}: Use qname list ${qnames_lst} to extract from bam file ${bam} to paired fastq file ${fastq_1} and ${fastq_2}"
	ls -lh ${fastq_1}
	ls -lh ${fastq_2}

	module unload BEDTools
    trap "rm -f ${tmp_bam} ${tmp_bam::-4}.tmp.bam" RETURN
}


function extract_read_from_fastq_by_qname () {
    local input_fastq=${1}
    local qname_list=${2}
    local output_fq=${3}
    # Output will be directly to stdout
    local qnames=$(mawk '{printf "%s,", $1;}' < ${qname_list} | mawk '{gsub(/,$/, ""); print;}')
    if [[ ! -z ${output_fq} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "These query names are going to be extracted from ${input_fastq} to ${output_fq}: "${qnames}
    fi
    # If qnames are too big we can't use -v arg of awk to import the data, otherwise it will cause an error of argument list too long.

    if [[ ${input_fastq} =~ \.gz$ ]]; then
        # Check input fastq gz format, correct it if not valid
        if gzip -t ${input_fastq}; then
            :
        elif [[ -f ${input_fastq::-3} ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: " Input fastq file not correctly compressed. ${input_fastq}"
            mv ${input_fastq} ${input_fastq}.backup
            gzip -f ${input_fastq::-3}
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: " Input fastq file not correctly compressed. ${input_fastq}"
        fi

        if [[ ! -z ${output_fq} ]]; then
            # zcat ${input_fastq} | mawk 'FNR == NR{qa[FNR] = $1;} NR > FNR && $1 ~ /^@/{for(i in qa) {q="@"qa[i]; if ($1 == q) {print;getline;print;getline;print;getline;print;}}}' ${qname_list} - > ${output_fq}
            zcat ${input_fastq} | mawk 'FNR == NR && $1 ~ /^@/{qn = $1; fl = $0; getline; sl = $0; getline; tl = $0; getline; ll = $0; total=fl"\n"sl"\n"tl"\n"ll; qa[qn] = total;} NR > FNR{printf "%s\n", qa["@"$1];}' - ${qname_list} | mawk 'NF > 0 {print;}' > ${output_fq} && \
            ls -lh ${output_fq}
        else
            # zcat ${input_fastq} | mawk 'FNR == NR{qa[FNR] = $1;} NR > FNR && $1 ~ /^@/{for(i in qa) {q="@"qa[i]; if ($1 == q) {print;getline;print;getline;print;getline;print;}}}' ${qname_list} - 
            zcat ${input_fastq} | mawk 'FNR == NR && $1 ~ /^@/{qn = $1; fl = $0; getline; sl = $0; getline; tl = $0; getline; ll = $0; total=fl"\n"sl"\n"tl"\n"ll; qa[qn] = total;} NR > FNR{printf "%s\n", qa["@"$1];}' - ${qname_list} | mawk 'NF > 0 {print;}' 
        fi
    else
        if [[ ! -z ${output_fq} ]]; then
            # mawk 'FNR == NR{qa[FNR] = $1;} NR > FNR && $1 ~ /^@/{for(i in qa) {q="@"qa[i]; if ($1 == q) {print;getline;print;getline;print;getline;print;}}}' ${qname_list} ${input_fastq} > ${output_fq}
            mawk 'FNR == NR && $1 ~ /^@/{qn = $1; fl = $0; getline; sl = $0; getline; tl = $0; getline; ll = $0; total=fl"\n"sl"\n"tl"\n"ll; qa[qn] = total;} NR > FNR{printf "%s\n", qa["@"$1];}' ${input_fastq} ${qname_list} | mawk 'NF > 0 {print;}' > ${output_fq} && \
            ls -lh ${output_fq}
        else
            # mawk 'FNR == NR{qa[FNR] = $1;} NR > FNR && $1 ~ /^@/{for(i in qa) {q="@"qa[i]; if ($1 == q) {print;getline;print;getline;print;getline;print;}}}' ${qname_list} ${input_fastq}
            mawk 'FNR == NR && $1 ~ /^@/{qn = $1; fl = $0; getline; sl = $0; getline; tl = $0; getline; ll = $0; total=fl"\n"sl"\n"tl"\n"ll; qa[qn] = total;} NR > FNR{printf "%s\n", qa["@"$1];}' ${input_fastq} ${qname_list} | mawk 'NF > 0 {print;}'
        fi
    fi
}


function strip_chr_vcf () {
    # Help stripping chr prefix for all contigs from vcf file. This function is needed before import vcf to mergeSVvcf
    local input_vcf=${1}
    local tmp_vcf=${input_vcf/.vcf*/.tmp.vcf}
    if [[ -z ${2} ]]; then local output_vcf=${1/.vcf*/.stripchr.vcf}; else local output_vcf=${2}; fi
    if [[ -z ${3} ]]; then local grch37_header=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.vcf.header; fi

    if [[ ${input_vcf} =~ \.[b]*gz$ ]]; then
        zcat ${input_vcf} | awk -F '\t' '{mod_line = gensub(/chr([a-zA-Z_0-9]+)/,"\\1","g",$0); print mod_line;}' > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    else
        awk -F '\t' '{mod_line = gensub(/chr([a-zA-Z_0-9]+)/,"\\1","g",$0); print mod_line;}' ${input_vcf} > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    fi
}


function add_chr_SV_vcf () {
    # Help stripping chr prefix for all contigs from vcf file. This function is needed before import vcf to mergeSVvcf
    local input_vcf=${1}
    local tmp_vcf=${input_vcf/.vcf*/.tmp.vcf}
    if [[ -z ${2} ]]; then local output_vcf=${1/.vcf*/.stripchr.vcf}; else local output_vcf=${2}; fi
    if [[ -z ${3} ]]; then local hg19_header=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.vcf.header; fi

    if [[ ${input_vcf} =~ \.[b]*gz$ ]]; then
        # Upon test, we do not need to escape < in awk regex
        zcat ${input_vcf} | awk -F '\t' '$0 ~ /^##/{ if($0 ~ /^##\<contig/) getline < "'${hg19_header}'"; print $0; } \
                                         $0 ~ /^#CHROM/{print;} \
                                         $0 !~ /^#/{mod_line = gensub(/([\[\]])([a-zA-Z_0-9]+):/, "\\1chr\\2:", "g", $0); printf "chr%s\n", mod_line;}' > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then 
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    else
        awk -F '\t' '$0 ~ /^##/{ if($0 ~ /^##\<contig/) getline < "'${hg19_header}'"; print $0; } \
                     $0 ~ /^#CHROM/{print;} \
                     $0 !~ /^#/{mod_line = gensub(/([\[\]])([a-zA-Z_0-9]+):/, "\\1chr\\2:", "g", $0); printf "chr%s\n", mod_line;}' ${input_vcf} > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then 
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    fi
}


function liftover_from_GRCh37_to_hg19 () {
    local input_vcf=${1}
    local tmp_vcf=${input_vcf/.vcf*/.tmp.vcf}
    if [[ -z ${2} ]]; then local output_vcf=${1/.vcf*/.stripchr.vcf}; else local output_vcf=${2}; fi
    if [[ ${input_vcf} =~ \.[b]*gz$ ]]; then
        # Upon test, we do not need to escape < in awk regex
        zcat ${input_vcf} | awk -F '\t' '$0 ~ /^##/{print;} \
                                         $0 ~ /^#CHROM/{print;} \
                                         $0 !~ /^#/{mod_line = gensub(/([\[\]])([a-zA-Z_0-9]+):/, "\\1chr\\2:", "g", $0); printf "chr%s\n", mod_line;}' > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then 
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    else
        awk -F '\t' '$0 ~ /^##/{print;} \
                     $0 ~ /^#CHROM/{print;} \
                     $0 !~ /^#/{mod_line = gensub(/([\[\]])([a-zA-Z_0-9]+):/, "\\1chr\\2:", "g", $0); printf "chr%s\n", mod_line;}' ${input_vcf} > ${tmp_vcf}
        if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then 
            bgzip -f ${tmp_vcf} && \
            bcftools sort -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
            tabix -f -p vcf ${tmp_vcf}.gz && mv ${tmp_vcf}.gz ${output_vcf} && mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi
        else
            mv ${tmp_vcf} ${output_vcf}
        fi
    fi
}


function check_vcf_ploidy(){
    local input_vcf=${1}
    module load bcftools 2> /dev/null
    bcftools +check-ploidy ${input_vcf} | \
    tail -1 | \
    awk '{print $5}'
    module unload bcftools 2> /dev/null
}


function post_process_pindel_vcf () {
    local pindel_vcf=${1}

    if [[ ${pindel_vcf} =~ \.vcf$ ]]; then
        local plain_input=${pindel_vcf}
        bgzip -f ${pindel_vcf} && tabix -f -p vcf ${pindel_vcf}.gz
        local pindel_vcf=${pindel_vcf}.gz
        local output_vcf=${pindel_vcf::-4}.tmp.vcf
    elif [[ ${pindel_vcf} =~ \.vcf\.gz$ ]]; then
        local plain_input=${pindel_vcf::-3}
        local output_vcf=${pindel_vcf::-7}.tmp.vcf
    fi

    module load bcftools

    bcftools norm -d indels --no-version -c s -f /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
    -Ov -o ${plain_input} ${pindel_vcf}

    ls -lh ${plain_input}

    module unload bcftools

    strip_chr_vcf ${plain_input} ${output_vcf} && mv ${output_vcf} ${plain_input}

    ls -lh ${plain_input}

    if [ -z ${output_vcf} ]; then
        touch ${output_vcf}
    else
        : > ${output_vcf}
    fi
    
    python3 /paedyl01/disk1/yangyxt/ngs_scripts/preprocess_pindel_vcf_before_merging.py \
    ${plain_input} \
    ${output_vcf}

    ls -lh ${output_vcf}
    mv ${output_vcf} ${plain_input} && bgzip -f ${plain_input} && tabix -f -p vcf ${plain_input}.gz
}


function prepare_intrinsic_vcf() {
	local homo_bed=${1} #2 can be the variable name to get the output vcf path
	if [[ -z ${homo_bed} ]];then local homo_bed="/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/hg19/homo_region.expanded.hg19.geneanno.PID.trimmed.bed"; fi
	local ref_assembly="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"
	local whole_region_bed="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.bed"
	
	# make masked genome
    local masked_genome
	bash ${central_scripts}/realign_masked_and_HC_multiploidy.sh \
	prepare_masked_genome \
	-g ${ref_assembly} \
	-s ${homo_bed} \
	-t ${whole_region_bed} \
	-m masked_genome

	# make fastq from ref assembly sequence
	# first extract the sequence from + strand homo regions.
	module load BEDTools
	module load seqtk
	mawk -F '\t' '$10 == "+"{printf "%s\t%s\t%s\n", $4, $5, $6;}' ${homo_bed} > ${homo_bed::-4}.nonsig.forward.bed
	mawk -F '\t' '$10 != "+"{printf "%s\t%s\t%s\n", $4, $5, $6;}' ${homo_bed} > ${homo_bed::-4}.nonsig.reverse.bed

	bedtools getfasta -fi ${ref_assembly} -bed ${homo_bed::-4}.nonsig.forward.bed > ${homo_bed::-4}.nonsig.onlyfor.fq && \
	bedtools getfasta -fi ${ref_assembly} -bed ${homo_bed::-4}.nonsig.reverse.bed | seqtk seq -r - >> ${homo_bed::-4}.nonsig.onlyfor.fq

	module unload seqtk
	module unload BEDTools
	
	# Now map the genome assembly sequence to masked assembly (only significant regions (homologous pair component) are not masked)
	module load minimap2
	minimap2 -k 27 -uf -ax asm5 ${masked_genome} ${homo_bed::-4}.nonsig.onlyfor.fq | samtools sort -O bam -o ${homo_bed::-4}.bam && \
    samtools index ${homo_bed::-4}.bam
	module unload minimap2

    module load bcftools
    bcftools mpileup --annotate INFO/AD -f ${masked_genome} ${homo_bed::-4}.bam | bcftools call -mv -Oz -o ${homo_bed::-4}.vcf.gz && \
    bcftools sort -Oz -o ${homo_bed::-4}.sort.vcf.gz ${homo_bed::-4}.vcf.gz && \
	${gatk} LeftAlignAndTrimVariants \
	-R ${ref_assembly} \
	-V ${homo_bed::-4}.sort.vcf.gz \
	-O ${homo_bed::-4}.vcf.gz && \
	bcftools norm -f /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
    -m -any --no-version \
    -Oz -o ${homo_bed::-4}.bial.vcf.gz ${homo_bed::-4}.vcf.gz && \
    mv ${homo_bed::-4}.bial.vcf.gz ${homo_bed::-4}.vcf.gz && \
	if [[ ! -z ${2} ]]; then eval ${2}="'${homo_bed::-4}.vcf.gz'"; fi
    module unload bcftools
}


function search_by_patientID () {
    local OPTIND s o t
    while getopts s:o:t:: args
    do
        case ${args} in
            s) local search_table=$OPTARG ;; # This table should have one ID per line
            o) local output_table=$OPTARG ;;
            t) local target_file=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done
    
    local wes_total_ped=/paedyl01/disk1/yangyxt/wes/wes_total_ped.ped
    local wgs_total_ped=/paedyl01/disk1/yangyxt/wgs/wgs_total_ped.ped
    local wesplus_total_ped=/paedyl01/disk1/yangyxt/wesplus/wesplus_total_ped.ped

    if [[ -f ${output_table} ]]; then rm ${output_table}; fi

    awk_merge_two_tables ${search_table} ${wes_total_ped} 1 2 ${output_table}.wes
    awk_merge_two_tables ${search_table} ${wgs_total_ped} 1 2 ${output_table}.wgs
    awk_merge_two_tables ${search_table} ${wesplus_total_ped} 1 2 ${output_table}.wesplus

    if [[ -z ${target_file} ]]; then
        awk -F '\t' 'NF > 2{print;}' ${output_table}.wes | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{print;}' ${output_table}.wgs | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{print;}' ${output_table}.wesplus | cut -f 1 --complement >> ${output_table}
    elif [[ ${target_file} == "vcf" ]]; then
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wes/%s/vcfs/%s.vcf.gz\n", $0, $8, $2;}' ${output_table}.wes | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wgs/%s/vcfs/%s.vcf.gz\n", $0, $8, $2;}' ${output_table}.wgs | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wesplus/%s/vcfs/%s.vcf.gz\n", $0, $8, $2;}' ${output_table}.wesplus | cut -f 1 --complement >> ${output_table}
    elif [[ ${target_file} == "bam" ]]; then
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wes/%s/aligned_results/%s.bqsr.bam\n", $0, $8, $3;}' ${output_table}.wes | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wgs/%s/aligned_results/%s.bqsr.bam\n", $0, $8, $3;}' ${output_table}.wgs | cut -f 1 --complement >> ${output_table}
        awk -F '\t' 'NF > 2{printf "%s\t/paedyl01/disk1/yangyxt/wesplus/%s/aligned_results/%s.bqsr.bam\n", $0, $8, $3;}' ${output_table}.wesplus | cut -f 1 --complement >> ${output_table}
    fi

    rm ${output_table}.*
}


function build_index_bwtsw(){
    local tar_fasta=${1}
    local fasta_name=$(basename ${tar_fasta})
    local tmp_dir=${2}
    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    if  [[ ! -f ${tar_fasta}.pac ]] || [[ ${tar_fasta} -nt ${tar_fasta}.pac ]] || \
        [[ ! -f ${tar_fasta}.bwt ]] || [[ ${tar_fasta} -nt ${tar_fasta}.bwt ]] || \
        [[ ! -f ${tar_fasta}.ann ]] || [[ ${tar_fasta} -nt ${tar_fasta}.ann ]] || \
        [[ ! -f ${tar_fasta}.amb ]] || [[ ${tar_fasta} -nt ${tar_fasta}.amb ]] || \
        [[ ! -f ${tar_fasta}.sa ]] || [[ ${tar_fasta} -nt ${tar_fasta}.sa ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Some of the bwa-generated index file is older than fasta file itself, need to get the index files updated." 
        echo "In case some other jobs are also trying to access the fasta file and related index file. We build the fasta under a temp directory."
        rsync -avu --delay-updates --include="${fasta_name}" --exclude="*" $(dirname ${tar_fasta})/ ${tmp_dir}/ && \
        /usr/bin/time bwa index -a bwtsw -b 380000000 ${tmp_dir}/${fasta_name} && \
        ls -lh ${tmp_dir}/${fasta_name}* && \
        rsync -avu --delay-updates --include="${fasta_name}*" --exclude="*" ${tmp_dir}/ $(dirname ${tar_fasta})/ && \
        ls -lh ${tar_fasta}*
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": All the bwa-generated index file is already built."
    fi

    trap "rm -rf ${tmp_dir}" RETURN
}

function build_index_bwamem2(){
    local tar_fasta=${1}
    local fasta_name=$(basename ${tar_fasta})
    local tmp_dir=${2}
    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    if [[ ! -f ${tar_fasta}.0123 ]] || [[ ${tar_fasta} -nt ${tar_fasta}.0123 ]] || \
       [[ ! -f ${tar_fasta}.bwt.2bit.64 ]] || [[ ${tar_fasta} -nt ${tar_fasta}.bwt.2bit.64 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Some of the bwa mem2 generated index file is older than fasta file itself, need to get the index files updated." 
        local bwa2=/home/yangyxt/software/bwa-mem2-2.1_x64-linux/bwa-mem2
        rsync -avu --delay-updates --include="${fasta_name}" --exclude="*" $(dirname ${tar_fasta})/ ${tmp_dir}/ && \
        /usr/bin/time ${bwa2} index ${tmp_dir}/${fasta_name} && \
        ls -lh ${tmp_dir}/${fasta_name}.bwt.2bit.64 && \
        rsync -avu --delay-updates --include="${fasta_name}*" --exclude="*" ${tmp_dir}/ $(dirname ${tar_fasta})/ && \
        ls -lh ${tar_fasta}.bwt.2bit.64
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": All the bwa-mem2-generated index file is already built."
    fi

    trap "rm -rf ${tmp_dir}" RETURN
}

function build_fai_for_genome() {
    local tar_fasta=${1}
    local fasta_name=$(basename ${tar_fasta})
    local genome_dir=$(dirname ${tar_fasta})
    local tmp_dir=${2}
    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    module load samtools
    rsync -avu --delay-updates --include="${fasta_name}" --exclude="*" ${genome_dir}/ ${tmp_dir}/
    samtools faidx ${tmp_dir}/${fasta_name} && ls -lh ${tmp_dir}/${fasta_name}.fai && \
    rsync -avu --delay-updates --include="${fasta_name}.fai" --exclude="*" ${tmp_dir}/ ${genome_dir}/ && \
    ls -lh ${genome_dir}/${fasta_name}.fai
    module unload samtools
}


function build_dict_for_genome() {
    local tar_fasta=${1}
    local fasta_name=$(basename ${tar_fasta})
    local genome_dir=$(dirname ${tar_fasta})
    local tmp_dir=${2}
    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    rsync -avu --delay-updates --include="${fasta_name}" --exclude="*" ${genome_dir}/ ${tmp_dir}/ && \
    ${gatk} CreateSequenceDictionary -R ${tmp_dir}/${fasta_name} && \
    ls -lh ${tmp_dir}/${fasta_name%.f*a*}.dict && \
    rsync -avu --delay-updates --include="${fasta_name%.f*a*}.dict" --exclude="*" ${tmp_dir}/ ${genome_dir}/ && \
    ls -lh ${genome_dir}/${fasta_name%.f*a*}.dict
}


function build_index_for_genome() {
    local genome=${1}
    local fasta_name=$(basename ${genome})
    local genome_dir=$(dirname ${genome})
    local tmp_dir=${2}
    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    build_index_bwtsw ${genome} ${tmp_dir}
    check_return_code
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Let's check the index files $(ls -lh ${genome}.*)"

    if [[ ! -d ${tmp_dir} ]]; then mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}; fi

    # Create index for masked genome
    if [[ ! -f ${genome}.fai ]] || [[ ${genome} -nt ${genome}.fai ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Fai index file for fasta ${genome} is not valid. Recreating now."
        build_fai_for_genome ${genome} ${tmp_dir}
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Fai index file for fasta ${genome} has been already created."
    fi

    if [[ ! -f ${genome%.f*a*}.dict ]] || [[ ${genome} -nt ${genome%.f*a*}.dict ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Dict index file for fasta ${genome} is not valid. Recreating now."
        build_dict_for_genome ${genome} ${tmp_dir}
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Dict file for fasta file ${genome} already created."
    fi

    trap "rm -rf ${tmp_dir}" RETURN
}


function independent_dragmap_alignment(){
    local OPTIND f r s g o t c
    while getopts f::r::s::g::o::t::c:: args
    do
        case ${args} in
            f) local forward_reads=$OPTARG ;;
            r) local reverse_reads=$OPTARG ;;
            c) local cover=$OPTARG ;;
            g) local ref_genome_dir=$OPTARG ;;
            s) local samp_ID=$OPTARG ;;
            o) local output_align=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    # At least input forward read and reverse read and output align

    if [[ -z ${ref_genome_dir} ]]; then local ref_genome_dir=/paedyl01/disk1/yangyxt/indexed_genome/hg19; fi
    if [[ -z ${samp_ID} ]]; then local samp_ID=$(basename ${forward_reads} | awk '{gsub(/_[a-z]*1\.f[ast]*q[\.gz]*$/, "", $0); printf "%s", $0;}'); fi
    if [[ -z ${output_align} ]]; then local output_align=${samp_ID}.bam; fi
    if [[ -z ${cover} ]]; then local cover="Yes"; fi
    
    local cpu=$(get_pbs_cpu)
    local mem=$(get_pbs_mem)
    if [[ ! -z ${mem} ]]; then echo "In function ${FUNCNAME}: $(timestamp): Under current PBS job ${PBS_JOBID}, the total cpu is ${cpu}, and the total mem is ${mem}Gb."; fi
    if [[ -z ${cpu} ]]; then local cpu=8; fi
    if [[ -z ${threads} ]]; then local threads=${cpu}; fi
    if [[ ${threads} -gt 10 ]]; then local threads=10; fi

    rm -f ${output_align}.tmp.*.bam 2> /dev/null && \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are some temporary files left in last alignment, remove them before alignment " || \
    >&2 echo "$(timestamp): In function ${FUNCNAME}, there are no temporary tmp.00x.bam file left in last alignment. "
    
    if check_bam_validity ${output_align} && [[ ${cover} == "No" ]]; then
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: Output alignment file ${output_align} already generated."
    else
        module load boost
        module load dragmap
        
        # Backup parameters
        # --Aligner.sec-phred-delta 0
        # --Aligner.pe-stat-mean-read-len 150
        # --Aligner.pe-stat-stddev-insert 70
        # --Aligner.pe-stat-mean-insert 300
        dragen-os -r ${ref_genome_dir} \
        -1 ${forward_reads} \
        -2 ${reverse_reads} \
        --num-threads ${threads} \
        --verbose \
        --RGID ${samp_ID} \
        --RGSM ${samp_ID} | \
        mawk -F '\t' '{gsub(/\/[1-2]$/, "", $1); print;}' | \
        samtools view -uSh -@ ${threads} - | samtools sort -O bam -@ ${threads} -o ${output_align} - && \
        samtools index ${output_align} && echo "DRAGMAP alignment done"

        module unload dragmap
        module unload boost
    fi
}


function dragmap_build_hashtable() {
    local ref_genome=${1}
    local output_dir=$(dirname ${ref_genome})
    if [[ -z ${2} ]]; then local prefix=$(basename ${ref_genome/.fasta*/}); else local prefix=${2}; fi
    local cpu=$(get_pbs_cpu)
    local mem=$(get_pbs_mem)
    if [[ ! -z ${mem} ]]; then echo "In function ${FUNCNAME}: $(timestamp): Under current PBS job ${PBS_JOBID}, the total cpu is ${cpu}, and the total mem is ${mem}Gb."; fi
    if [[ -z ${cpu} ]]; then local cpu=8; fi

    module load boost
    module load dragmap

    if [[ -f ${output_dir}/hash_table.cfg ]] || \
    [[ -f ${output_dir}/hash_table.cfg.bin ]] || \
    [[ -f ${output_dir}/hash_table.cmp ]] || \
    [[ -f ${output_dir}/hash_table_stats.txt ]] || \
    [[ -f ${output_dir}/reference.bin ]] || \
    [[ -f ${output_dir}/ref_index.bin ]] || \
    [[ -f ${output_dir}/repeat_mask.bin ]] || \
    [[ -f ${output_dir}/str_table.bin ]]; then
        local tmp_dir_name=$(randomID)
        /usr/bin/time dragen-os \
        --num-threads ${cpu} \
        --build-hash-table "true" \
        --ht-reference ${ref_genome} \
        --output-directory ${output_dir}/${tmp_dir_name} \
        --output-file-prefix ${prefix} && \
        rsync --rsync-path=/bin/rsync -avu ${output_dir}/${tmp_dir_name}/ ${output_dir}/
        trap "rm -rf ${output_dir}/${tmp_dir_name}" RETURN
    else
        /usr/bin/time dragen-os \
        --num-threads ${cpu} \
        --build-hash-table "true" \
        --ht-reference ${ref_genome} \
        --output-directory ${output_dir} \
        --output-file-prefix ${prefix}
    fi

    module unload dragmap
    module unload boost
}


function get_pbs_cpu(){
    # echo "In function ${FUNCNAME}: $(timestamp): The current PBS job we are in is $PBS_JOBID"
    /paedyl01/disk1/yangyxt/ngs_scripts/qst | \
    awk -F '\t' '{new = gensub(/ +/, "\t", "g", $0); new = gensub(/^\t/, "", "g", new); print new;}' | \
    mawk -F '\t' '$1 == "'${PBS_JOBID}'"{printf "%s", $5;}' | \
    mawk -F '=' '{printf "%s", $2;}' | \
    mawk -F ':' '{printf "%s", $1;}'
}


function get_pbs_mem(){
    # The unit is Gb
    /paedyl01/disk1/yangyxt/ngs_scripts/qst | \
    awk -F '\t' '{new = gensub(/ +/, "\t", "g", $0); new = gensub(/^\t/, "", "g", new); print new;}' | \
    mawk -F '\t' '$1 == "'${PBS_JOBID}'"{printf "%s", $5;}' | \
    mawk -F '=' '{printf "%s", $3;}' | \
    mawk -F 'Gb' '{printf "%s", $1;}'
}


function change_qual_code(){
    local input_fq=${1}

    if [[ ${input_fq} =~ \.gz$ ]]; then
        local p_fq=${input_fq/.gz/}
        local gz_fq=${input_fq}
    else
        local p_fq=${input_fq}
        local gz_fq="${input_fq}.gz"
    fi

    # echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): plain fastq file is ${p_fq}."
    # echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): gzip fastq file is ${gz_fq}."

    local mawk=$(which mawk 2> /dev/null)
    if [[ -z ${mawk} ]]; then local mawk="/home/yangyxt/anaconda3/bin/mawk"; fi
    if [[ -f ${p_fq} ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): plain fastq file existed: $(ls -lh ${p_fq})"
        if [[ ! -f ${p_fq}.tmp ]]; then touch ${p_fq}.tmp; fi 
        $mawk '{ \
            if((NR+2)/4 == int((NR+2)/4)) {\
                seq_len = length($0); \
                print;} \
            else if(NR/4 == int(NR/4)) {\
                qual_len = length($0); \
                gsub(/[\:\;\<\=\>\?\@]/, "A", $0); \
                gsub(/[^A-Za-z0-9]/, ",", $0); \
                if (seq_len == qual_len) {\
                    print;} \
                else if(seq_len > qual_len) { \
                    delta = seq_len - qual_len; \
                    printf "%s", $0; \
                    for (i=1;i<=delta;i++) printf ","; \
                    printf "\n";} \
                else { \
                    trimmed = substr($0, 1, seq_len); \
                    printf "%s\n", trimmed;}} \                  
            else {print;}}' ${p_fq} > ${p_fq}.tmp && \
        mv ${p_fq}.tmp ${p_fq}
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): New plain fastq file generated: $(ls -lh ${p_fq})"
    elif [[ -f ${gz_fq} ]]; then
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): gzip fastq file existed: $(ls -lh ${gz_fq})"
        touch ${p_fq}
        zcat ${gz_fq} | $mawk '{ \
            if((NR+2)/4 == int((NR+2)/4)) {\
                seq_len = length($0); \
                print;} \
            else if(NR/4 == int(NR/4)) {\
                qual_len = length($0); \
                gsub(/[\:\;\<\=\>\?\@]/, "A", $0); \
                gsub(/[^A-Za-z0-9]/, ",", $0); \
                if (seq_len == qual_len) {\
                    print;} \
                else if(seq_len > qual_len){ \
                    delta = seq_len - qual_len; \
                    printf "%s", $0; \
                    for (i=1;i<=delta;i++) printf ","; \
                    printf "\n";} \
                else { \
                    trimmed = substr($0, 1, seq_len); \
                    printf "%s\n", trimmed;}} \                    
            else {print;}}' > ${p_fq} && \
        gzip -f -c ${p_fq} > ${gz_fq/.f*q*/.tmp.fq.gz} && \
        mv ${gz_fq/.f*q*/.tmp.fq.gz} ${gz_fq} && \
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): New gzip fastq file generated: $(ls -lh ${gz_fq})"
    else
        return 1
    fi
}


function visualize_insert_size_distribution(){
    source "$(which conda | awk -F '/' '{for (i=1;i<NF-1;i++) printf "%s/", $i;}')etc/profile.d/conda.sh"
    conda activate r_env

    local input_bam=${1}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    ${gatk} CollectInsertSizeMetrics \
    -I ${input_bam} \
    -O ${input_bam/.bam/insert_size_metrics.txt} \
    -H ${input_bam/.bam/insert_size_histogram.pdf} \
    -M 0 \
    -W 1500

    source "$(which conda | awk -F '/' '{for (i=1;i<NF-1;i++) printf "%s/", $i;}')etc/profile.d/conda.sh"
    conda deactivate
}


function parallel_visualize_insert_size_distribution(){
    local wkd=${1}
    local threads=${2}
    local -a bams=($(find ${wkd} -regex "${wkd}/[A-Za-z0-9_-]+\.bam"))

    if [[ -z ${threads} ]]; then
        local cpu=$(get_pbs_cpu)
        if [[ ${cpu} -ge 10 ]]; then
            local threads=10
        elif [[ ${cpu} -gt 1 ]]; then
            local threads=$((cpu-1))
        else
            local threads=1
        fi
    fi

    module load parallel
    source $(which env_parallel.bash)
    env_parallel --env visualize_insert_size_distribution -j${threads} -k --dry-run \
    visualize_insert_size_distribution {} ::: ${bams[*]} && \
    env_parallel --env visualize_insert_size_distribution -j${threads} -k \
    visualize_insert_size_distribution {} ::: ${bams[*]}

    module unload parallel
}


function check_bam_total_seq(){
    local input_bam=${1}
    samtools stats ${input_bam} | mawk -F '\t' '$2 == "raw_total_sequences"{printf "%s",$3; exit 0;}'
}


function fade_per_bam(){
    local OPTIND i o r t f
    while getopts i:o::r::t::f:: args
    do 
        case ${args} in
            i) local input_bam=$OPTARG ;; 
            o) local output_bam=$OPTARG ;;
            r) local ref_genome=$OPTARG ;;
            t) local threads=$OPTARG ;;
            f) local force_replace=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
    
    if [[ -z ${ref_genome} ]]; then
        local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta
    fi

    if [[ -z ${output_bam} ]]; then
        local output_bam=${input_bam/.bam/.filtered.bam}
    fi

    if [[ -z ${threads} ]]; then
        local threads=1
    fi

    if [[ ! -f ${input_bam}.bai ]] && [[ ! -f ${input_bam/.bam/.bai} ]]; then
        samtools index ${input_bam}
    fi
    
    if check_bam_validity ${input_bam}; then
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: Start to use FADE to process BAM file ${input_bam}"
        fade annotate -b ${input_bam} ${ref_genome} > ${input_bam/.bam/.anno.bam} && \
        samtools sort -n -O bam -o ${input_bam/.bam/.anno.qsort.bam} -@ {threads} ${input_bam/.bam/.anno.bam} && \
        fade out -b ${input_bam/.bam/.anno.qsort.bam} > ${input_bam/.bam/.tmp.bam} && \
        samtools sort -O bam -o ${output_bam} -@ ${threads} ${input_bam/.bam/.tmp.bam} && \
        samtools index ${output_bam} && \
        ls -lh ${output_bam} || \
        { >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: Failed to execute FADE now. It's currently not stable. If failed, we just skipped this step." && return; }
    fi

    if [[ ! -z ${force_replace} ]] && check_bam_validity ${output_bam}; then
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: ${output_bam} valid and we set up parameter replace: ${force_replace}"
        mv ${output_bam} ${input_bam} && \
        ls -lh ${input_bam} && \
        mv ${output_bam}.bai ${input_bam}.bai && \
        ls -lh ${input_bam}.bai
    elif check_bam_validity ${output_bam}; then
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: ${output_bam} valid but replace argument not set up."
    else
        >&2 echo "Line ${LINENO}: $(timestamp): In function ${FUNCNAME}: ${output_bam} not valid. Maybe due to the failed execution of FADE."
    fi 

    trap "rm ${input_bam/.bam/.anno.bam} ${input_bam/.bam/.anno.qsort.bam} ${input_bam/.bam/.tmp.bam}" 0
}


function parallel_check_bam_total_seq(){
    local wkd=${1}
    local threads=${2}
    local -a bams=($(find ${wkd} -regex "${wkd}/[A-Za-z0-9_-]+\.bam"))

    if [[ -z ${threads} ]]; then
        local cpu=$(get_pbs_cpu)
        if [[ ${cpu} -ge 10 ]]; then
            local threads=10
        elif [[ ${cpu} -gt 1 ]]; then
            local threads=$((cpu-1))
        else
            local threads=1
        fi
    fi

    module load parallel
    source $(which env_parallel.bash)
    env_parallel --env visualize_insert_size_distribution -j${threads} -k --dry-run \
    visualize_insert_size_distribution {} ::: ${bams[*]} && \
    env_parallel --env visualize_insert_size_distribution -j${threads} -k \
    visualize_insert_size_distribution {} ::: ${bams[*]}

    module unload parallel
}


function batch_check_bam_total_seq(){
    local -a wkds=($(echo ${1} | awk 'BEGIN{RS=",";} {printf "%s ", $1;}'))
    local -a dir_names
    local output_tab=${2}

    for wkd in "${wkds[@]}"; do
        dir_names+=( "$(dirname ${wkd})" )
        parallel_check_bam_total_seq ${wkd} > ${output_tab}.tmp
        if [[ $(cat ${output_tab}.tmp | wc -l) -gt $(cat ${output_tab}.tmp | wc -l) ]]; then
            local new_lines=$(cat ${output_tab}.tmp | wc -l)
            local old_lines=$(cat ${output_tab} | wc -l)
            local delta=$((new_lines - old_lines))
            mawk -F '\t' -v delta="${delta}" '\
                {print; nf = NF;} \
                END{ for(i=1;i<=delta;i++) { \
                        for(i=1;i<nf;i++) {printf "\t";} \
                        printf "\n";}}' ${output_tab} > ${output_tab}.new && \
            mv ${output_tab}.new ${output_tab}
        fi
        mawk -F '\t' \
        'NR == FNR{arr[FNR] = $0;} NR > FNR{printf "%s\t%s\n", $0, arr[FNR];}' \
        ${output_tab}.tmp ${output_tab} > ${output_tab}
    done

    echo ${dir_names[*]} | mawk '{for(i=1;i<NF;i++) {printf "%s\t", $i;} printf "%s\n", $NF;}' > ${output_tab}.new && \
    cat ${output_tab}.new ${output_tab} > ${output_tab}.final && \
    mv ${output_tab}.final ${output_tab}

    trap "rm ${output_tab}.*" RETURN
}


function parallel_run_local_func(){
    local OPTIND f a t c
    while getopts f:a:t::c:: args
    do 
        case ${args} in
            f) local func=$OPTARG ;; 
            a) local args=$OPTARG ;;
            c) local constant_args=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${threads} ]]; then
        local cpu=$(get_pbs_cpu)
        if [[ ${cpu} -ge 10 ]]; then
            local threads=10
        elif [[ ${cpu} -gt 1 ]]; then
            local threads=$((cpu-1))
        else
            local threads=1
        fi
    fi

    if [[ ! -z ${constant_args} ]]; then
        local constant_args=$(echo ${constant_args} | awk -F ',' 'BEGIN{RS=",";} {printf "%s ", $1;}')
    fi
    
    module load parallel
    source $(which env_parallel.bash)

    find_involved_vars "${func}" __involved_vars
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": ENV/SHELL variables involved in function ${func} are: "${__involved_vars}
	parallel_args=$(echo ${__involved_vars} | awk 'BEGIN{RS=" ";} { gsub(/\n$/, ""); printf "--env %s ", $1;}')

    env_parallel ${parallel_args} --dry-run -j${threads} -k --link \
    ${func} ${constant_args} {} ::: $(echo ${args} | awk -F ',' 'BEGIN{RS=",";} {printf "%s ", $1;}') && \
    env_parallel ${parallel_args} -j${threads} -k --link \
    ${func} ${constant_args} {} ::: $(echo ${args} | awk -F ',' 'BEGIN{RS=",";} {printf "%s ", $1;}')

    module unload parallel
}


function check_parallel_joblog() {
    local ret_code=$(echo $?)
    local job_log=${1}
    # print failed job commands to stdout
    local job_num=$(tail -n +2 ${job_log} | wc -l)
    local -a fail_job_ids=($(tail -n +2 ${job_log} | awk -F '\t' '$7 != "0"{print $1;}'))
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Parallel job log is ${job_log}"
    if [[ ${#fail_job_ids[@]} -gt 0 ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Totally ${job_num} jobs were running in parallel, ${#fail_job_ids[@]} of the parallel jobs failed, get a return code of ${ret_code}, the commands are:"
        for id in "${fail_job_ids[@]}"; do
            >&2 mawk -F '\t' '$1 == "'${id}'"{printf "%s\n", $NF;}' ${job_log}
            mawk -F '\t' '$1 == "'${id}'"{printf "%s\n", $NF;}' ${job_log}
        done
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": All parellel jobs are finished without an error."
        >&2 cat ${job_log}
    fi
}


function merge_variants_with_priority() {
    local OPTIND r p o f
    while getopts r:p:o:f:: args
    do 
        case ${args} in
            f) local force_replace=$OPTARG ;; 
            p) local prior_vcf=$OPTARG ;;
            r) local raw_vcf=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done
	local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts

	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): Original vcf input is ${raw_vcf}"
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): Priority vcf input is ${prior_vcf}"
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): Output vcf is ${output_vcf}"

	if check_vcf_validity ${output_vcf} && [[ -z ${force_replace} ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): The output vcf file ${output_vcf} is already there. Quit redoing it all over again."
		ls -lh ${output_vcf}
	else
		if [[ ! -z ${force_replace} ]]; then
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): The output vcf file ${output_vcf} is not ready."
		else
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "$(timestamp): The output vcf file ${output_vcf} is already there. But force_replace parameter is set. Re-doing the merging again."
		fi
		python3 ${central_scripts}/merge_variants_with_priority.py \
		-ov ${raw_vcf} -pv ${prior_vcf} -op ${output_vcf}
	fi
}


function quick_get_bam_depth_MQ() {
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
    local mosdepth=/home/yangyxt/software/mosdepth/mosdepth
    local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
    if [[ -z ${MQ_threshold} ]]; then local MQ_threshold=10; fi

    if [[ ! -z ${target_region} ]]; then
        local target_region_ID=$(basename ${target_region} | awk -F '.' '{printf "%s", $1;}')
        local targeted_bam=${input_bam/.bam/}.${target_region_ID}.bam
        bash ${central_scripts}/select_bam_by_regions.sh \
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
        ${mosdepth} -Q ${MQ_threshold} -b ${target_region} ${input_bam/.bam/.highMQ} ${input_bam} && \
        ls -lh ${input_bam/.bam/.highMQ}.per-base.bed.gz || \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": ${mosdepth} still out of memory."
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Command line is: ${mosdepth} -Q ${MQ_threshold} ${input_bam/.bam/.highMQ} ${input_bam}"
        ${mosdepth} -Q ${MQ_threshold} ${input_bam/.bam/.highMQ} ${input_bam} && \
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


function pick_poorcov_region_bam() {
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
    local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
    local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed
    if [[ -z ${MQ_threshold} ]]; then local MQ_threshold=10; fi
    if [[ -z ${depth_threshold} ]]; then local depth_threshold=10; fi

    if [[ -z ${target_region} ]]; then
        quick_get_bam_depth_MQ -i ${input_bam} -o ${tmp_bed} -q ${MQ_threshold}
    else
        quick_get_bam_depth_MQ -i ${input_bam} -t ${target_region} -o ${tmp_bed} -q ${MQ_threshold}
    fi
    ls -lh ${tmp_bed} && head ${tmp_bed}
    pick_poorcov_region_bed -i ${tmp_bed} -o ${output_bed} -d ${depth_threshold}
    trap "rm ${tmp_bed}" RETURN
    ls -lh ${output_bed}    
}


function download() {
    local url=${1}
    local md5=${2}
    local file_name=$(basename ${url})

    if [[ ! -f $(pwd)/${file_name} ]]; then
        touch $(pwd)/${file_name}
    fi

    local actual_md5=$(md5sum $(pwd)/${file_name} | awk '{print $1}')
    if [[ -z ${md5} ]]; then
        while [[ ${new_md5} != ${actual_md5} ]]; do
            local actual_md5=${new_md5}
            wget -c ${url}; \
            local new_md5=$(md5sum $(pwd)/${file_name} | awk '{print $1}')
        done
    else
        while [[ ${actual_md5} != ${md5} ]]; do
            wget -c ${url}; \
            local actual_md5=$(md5sum $(pwd)/${file_name} | awk '{print $1}')
        done
    fi
}


function retrieve_seq_from_fasta(){
    local region_str=${1}
    local ref_fasta=${2}

    if [[ -z ${ref_fasta} ]]; then
        local ref_fasta="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"
    fi

    samtools faidx ${region_str} ${ref_fasta} | awk '$0 ~ /^[A-Za-z]/{gsub(/\n$/, ""); printf "%s", $0;}'
}


function check_var_number(){
    local var_table=${1}
    local id_col=${2}

    if [[ -z ${id_col} ]]; then
        local id_col="uniq_ID"
    fi

    local id_col_ind=$(head -1 ${var_table} | mawk -F '\t' -v idc="${id_col}" '{for(i=1;i<=NF;i++) if($i == idc) {print i; exit 0;}}')
    if [[ -z ${id_col_ind} ]]; then
        first_5cols=$(head -1 ${var_table} | mawk -F '\t' '{printf "%s,%s,%s,%s,%s", $1, $2, $3, $4, $5;}')
        if [[ ${first_5cols} =~ ^[Cc][Hh][Rr],[Ss][Tt][Aa][Rr][Tt],[Ee][Nn][Dd],[Rr][Ee][Ff],[Aa][Ll][Tt]$ ]]; then
            local var_num=$(cut -f 1-5 ${var_table} | mawk 'NR>1{print;}' | sort - | uniq | wc -l)
            >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Table ${var_table} contains ${var_num} variant events. Determined by first 5 columns ${first_5cols}" 
        else
            >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Cannot identify the number of variant events in table ${var_table}. No variant ID nor chr,start,end,ref,alt combined columns are detected."
        fi
    else
        local var_num=$(mawk -F '\t' 'NR > 1{print $'${id_col_ind}'}' ${var_table} | sort - | uniq | wc -l)
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: Table ${var_table} contains ${var_num} variant events, determined by variant ID column ${id_col}" 
    fi
}


function pick_poorcov_region_bed() {
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
    local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
    if [[ -z ${depth_threshold} ]]; then local depth_threshold=10; fi
    # Function used to output a bed file with the 4th column recording depth information

    if [[ -z ${target_region} ]]; then
        if [[ ! ${output_bed} =~ \.bed(\.gz)*$ ]]; then
            local ob_var_name=${output_bed}
            local output_bed=${input_bed/.bed*/.lowDepth.bed}
        fi
        if [[ ${input_bed} =~ \.gz$ ]]; then
            zcat ${input_bed} | \
            mawk -F '\t' '$4 <= '${depth_threshold}' {print;}' > ${output_bed/.gz/}
        else
            mawk -F '\t' '$4 <= '${depth_threshold}' {print;}' ${input_bed} > ${output_bed/.gz/}
        fi
    else
        local target_region_ID=$(basename ${target_region} | awk -F '.' '{printf "%s", $1;}')
        if [[ ! ${output_bed} =~ \.bed(\.gz)*$ ]]; then
            local ob_var_name=${output_bed}
            local output_bed=${input_bed/.bed*/.lowDepth}.${target_region_ID}.bed
        fi
        bedtools intersect -wa -a ${input_bed} -b ${target_region} | \
        mawk -F '\t' '$4 <= '${depth_threshold}'{print;}' > ${output_bed/.gz/}
    fi

    if [[ ${output_bed} =~ \.gz$ ]]; then
        gzip -f -c ${output_bed/.gz/} > ${output_bed}
    fi

    if [[ ! -z ${ob_var_name} ]]; then
        eval ${ob_var_name}="'${output_bed}'"
    fi
}


function  mark_homoseq_recall_shortv(){
    local input_vcf=${1}
    local threads=${2}
    local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

    if [[ -z ${threads} ]]; then
        local threads=1
    fi
    local output_vcf=${input_vcf/.vcf*/.tmp.vcf}

    python3 ${central_scripts}/mark_recalled_shortv_from_homoseq.py \
    -i ${input_vcf} \
    -t ${threads} \
    -o ${output_vcf} && \
    ls -lh ${output_vcf} && \
    bgzip -f -c ${output_vcf} > ${output_vcf}.gz && \
    bcftools sort -Oz -o ${output_vcf}.gz ${output_vcf}.gz && \
    tabix -f -p vcf ${output_vcf}.gz && \
    ${gatk} LeftAlignAndTrimVariants \
    -R /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
    -V ${output_vcf}.gz \
    -O ${output_vcf}.tmp.gz \
    --split-multi-allelics \
    --max-indel-length 500 && \
    bcftools norm -d exact -c ws --no-version \
    -f /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
    -Oz -o ${output_vcf}.gz ${output_vcf}.tmp.gz && \
    bcftools sort -Oz -o ${output_vcf}.gz ${output_vcf}.gz && \
    tabix -f -p vcf ${output_vcf}.gz && \
    if check_vcf_validity ${output_vcf}.gz; then 
        if [[ ${input_vcf} =~ \.gz ]]; then
            mv ${output_vcf}.gz ${input_vcf} && \
            tabix -f -p vcf ${input_vcf} && \
            rm ${output_vcf}.gz.tbi
        else
            gunzip -c ${output_vcf}.gz > ${input_vcf} && \
            rm ${output_vcf}.gz.tbi
        fi
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Output vcf ${output_vcf}.gz is not valid, hence cannot cover input vcf ${input_vcf}"
        return 1
    fi
}


function check_homoseq_recall_variants(){
    local input_vcf=${1}

    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n" ${input_vcf} | \
    mawk -F '\t' '$5 ~ /INTRINSIC/ || $5 ~ /^\.$/{print;}' > ${input_vcf/.vcf*/.vcf.recalled_shortv.tsv}

    trap 'rm -f ${input_vcf/.vcf*/.vcf.recalled_shortv.tsv}' RETURN

    if [[ $(cat ${input_vcf/.vcf*/.vcf.recalled_shortv.tsv} | wc -l) -gt 0 ]]; then
        return 0
    else
        return 1
    fi
}


function determine_job_num(){
    local OPTIND m c r t
    while getopts m::c::r::t:: args
    do 
        case ${args} in
            m) local mem=$OPTARG ;; 
            c) local cpu=$OPTARG ;;
            r) local pbs_mem=$OPTARG ;;
            t) local pbs_cpu=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${cpu} ]]; then local cpu=1; fi

    if [[ -z ${pbs_cpu} ]]; then
        local pbs_cpu=$(get_pbs_cpu)
    fi
    local cpu_jobs=$(awk -v pu="${pbs_cpu}" -v u="${cpu}" 'BEGIN{printf "%i", pu/u;}')

    if [[ -z ${pbs_mem} ]]; then
        local pbs_mem=$(get_pbs_mem)
    fi
    local mem_jobs=$(awk -v pm="${pbs_mem}" -v m="${mem}" 'BEGIN{printf "%i", pm/m;}')

    if [[ ${cpu_jobs} -lt 1 ]]; then local cpu_jobs=1; fi
    if [[ ${mem_jobs} -lt 1 ]]; then local mem_jobs=1; fi

    >&2 echo "$(timestamp): ${FUNCNAME}: ${mem}Gb and ${cpu} vCPU for a job. We can have ${cpu_jobs} parallel jobs given cpu limitation and ${mem_jobs} parallel jobs given memory limitation."

    awk -v cj="${cpu_jobs}" -v mj="${mem_jobs}" 'BEGIN{printf "%s\n%s", cj, mj}' | \
    sort -n | head -1
}


function get_updated_progress_sheet (){
    local progress_tab=${1}
    local seq_type=${2}

    if [[ -z ${seq_type} ]]; then
        local seq_type=$(basename $(dirname ${progress_tab}))
    fi

    >&2 echo "$(timestamp): In function ${FUNCNAME}, Running: python3 ${central_scripts}/fetch_updated_form_response.py ${progress_tab} ${seq_type}"
    python3 ${central_scripts}/fetch_updated_form_response.py \
    ${progress_tab} \
    ${seq_type} && \
    ls -lh ${progress_tab}
}


function check_available_queue() {
    local required_cpu=${1}
    local required_mem=${2}
    local __output_que=${3}

    local tmp_table=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).tsv

    pbsnodes -aSj | \
    tail -n +4 | \
    awk '$1 ~ /^paed[a-z]*y/{gsub(/ +/, "\t"); print;}' | \
    awk -F '\t' '{gsub(/gb\/.+$/, "", $6); gsub(/\/.+$/, "", $7); printf "%s\t%s\t%s\t%s\n", $1, $2, $6, $7;}' > ${tmp_table}

    cat ${tmp_table}
    local -a avail_qs=($(awk -F '\t' -v rc="${required_cpu}" -v rm="${required_mem}" \
    '$2 == "free" && $4/rc >= 1 && $3/rm >= 1{printf "%s ", $1;}' ${tmp_table}))

    if [[ ${#avail_qs[@]} -eq 0 ]]; then
        eval ${__output_que}="large"
    else
        eval ${__output_que}="'${avail_qs}'"
    fi

    trap "rm -f ${tmp_table}" RETURN
}


function run_cmd_qsub() {
    local OPTIND c t m q n a
    while getopts c:t::m::q::n::a:: args
    do 
        case ${args} in
            c) local command_func=$OPTARG ;; 
            t) local threads=$OPTARG ;;
            m) local memory=$OPTARG ;;
            q) local queue=$OPTARG ;;
            n) local jobname=$OPTARG ;;
            a) local arguments=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local email_addr="u3005579@connect.hku.hk"

    if [[ -z ${threads} ]]; then
        local threads=8
    fi
    if [[ -z ${memory} ]]; then
        local memory=80
    fi
    if [[ -z ${queue} ]]; then
        local queue
        check_available_queue ${threads} ${memory} queue
    fi
    if [[ -z ${jobname} ]]; then
        local jobname=$(echo ${command_func} | awk '{gsub(/(bash)|(python) +/, ""); gsub(" ", "_"); print;}')
    fi

    local queue_name=$(echo ${queue} | awk '{gsub(/[0-9]+/, ""); print;}')

    if [[ ! -z ${arguments} ]]; then 
        local -a arg_arr=($(echo ${arguments} | awk -F ',' '{for(i=1;i<=NF;i++) {printf "%s ", $i;}}'))
        >&2 echo "Line ${LINENO}: $(timestamp): Input argument pairs are ${arg_arr[*]}"
        local argument=","
        for arg_pair in "${arg_arr[@]}"; do
            local argument="${argument}$(echo ${arg_pair} | awk -F '=' '{printf "%s=\"%s\",", $1, $2;}')"
        done
        local argument=${argument::-1}
    fi

    if [[ ${queue} == "large" ]]; then
        >&2 echo "Line ${LINENO}: $(timestamp): No queues are available now. Given required cpu is ${threads}, required mem is ${memory}, we assign this job ${jobname} job to node ${queue}."
        local qsub_cmd="qsub -q ${queue_name} -l mem=${memory}g,nodes=1:ppn=${threads},walltime=84:00:00 -N ${jobname} -j eo -e omics:/paedyl01/disk1/yangyxt/batch_run_log/${jobname}.log -v func=\"${command_func}\"${argument} /paedyl01/disk1/yangyxt/ngs_scripts/batch_run_pipeline.sh"
    else
        >&2 echo "Line ${LINENO}: $(timestamp): Given the required cpu is ${threads}, required mem is ${memory}, we assign this job ${jobname} job to node ${queue}."
        local qsub_cmd="qsub -q ${queue_name} -l mem=${memory}g,nodes=${queue}:ppn=${threads},walltime=9999:00:00 -N ${jobname} -j eo -e omics:/paedyl01/disk1/yangyxt/batch_run_log/${jobname}.log -v func=\"${command_func}\"${argument} /paedyl01/disk1/yangyxt/ngs_scripts/batch_run_pipeline.sh"
    fi

    >&2 echo "Line ${LINENO}: $(timestamp): log is /paedyl01/disk1/yangyxt/batch_run_log/${jobname}.log qsub command is ${qsub_cmd}"
    echo "Run job with qsub command: ${qsub_cmd}"$'\n'"Log file is /paedyl01/disk1/yangyxt/batch_run_log/${jobname}.log" | mail -s "Run job ${jobname} at PBS node ${queue_name}" ${email_addr} && \
    ${qsub_cmd}
}



function check_fam_anno_status() {
    local batch_dir=${1}
    local ped_file_name=${2}
    local batch_name=$(basename ${batch_dir})

    if [[ -z ${ped_file_name} ]]; then
        local ped_file_name=$(basename ${batch_dir})
    fi

    local ped_file=${batch_dir}/pedigree/${ped_file_name}.ped

    if [[ ! -f ${ped_file} ]]; then
        >&2 echo "Line ${LINENO}: $(timestamp): ${ped_file} does not exist. Check before acting"
        return 1
    else
        local -a batch_fams=($(mawk -F '\t' 'NR > 1 {print $1;}' ${ped_file} | uniq - | awk '{printf  "%s ", $1;}'))
    fi

    local -a further_need_fams
    for fam in "${batch_fams[@]}"; do
        if [[ -f ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log ]] && [[ -f ${batch_dir}/respective_scripts/${fam}_anno.log ]]; then
            if [[ ${batch_dir}/respective_scripts/${fam}_anno.log -ot ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log ]]; then
                local anno_status=$(mawk -F '\t' '$0 ~ /Traceback/{printf "1"; exit 0;} $0 ~ /Complete all the processing/{printf "0"; exit 0;} END{printf "1";}' ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log)
                if [[ ${anno_status} == "01" ]]; then
                    >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} finished annotation."
                else
                    >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} not finished annotation, by checking ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log"
                    further_need_fams+=( ${fam} )
                fi
            else
                local anno_status=$(mawk -F '\t' '$0 ~ /Traceback/{printf "1"; exit 0;} $0 ~ /Complete all the processing/{printf "0"; exit 0;} END{printf "1";}' ${batch_dir}/respective_scripts/${fam}_anno.log)
                if [[ ${anno_status} == "01" ]]; then
                    >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} finished annotation."
                else
                    >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} not finished annotation, by checking ${batch_dir}/respective_scripts/${fam}_anno.log"
                    further_need_fams+=( ${fam} )
                fi
            fi
        elif [[ -f ${batch_dir}/respective_scripts/${fam}_anno.log ]]; then
            local anno_status=$(mawk -F '\t' '$0 ~ /Traceback/{printf "1"; exit 0;} $0 ~ /Complete all the processing/{printf "0"; exit 0;} END{printf "1";}' ${batch_dir}/respective_scripts/${fam}_anno.log)
            if [[ ${anno_status} == "01" ]]; then
                >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} finished annotation."
            else
                >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} not finished annotation, by checking ${batch_dir}/respective_scripts/${fam}_anno.log"
                further_need_fams+=( ${fam} )
            fi
        else
            local anno_status=$(mawk -F '\t' '$0 ~ /Traceback/{printf "1"; exit 0;} $0 ~ /Complete all the processing/{printf "0"; exit 0;} END{printf "1";}' ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log)
            if [[ ${anno_status} == "01" ]]; then
                >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} finished annotation."
            else
                >&2 echo "Line ${LINENO}: $(timestamp): ${fam} in ${batch_name} not finished annotation, by checking ${batch_dir}/annotations/batch_anno_${batch_name}_fam_${fam}.log"
                further_need_fams+=( ${fam} )
            fi
        fi
    done

    if [[ ${#further_need_fams} -eq 0 ]]; then
        >&2 echo "Line ${LINENO}: $(timestamp): All fams in ${batch_name} finished annotation."
    else
        >&2 echo "Line ${LINENO}: $(timestamp): These fams in ${batch_name} not finished annotation. Now generate the command for batch execution."
        echo "${further_need_fams[@]}" | mawk 'BEGIN{RS=" "; FS="\t";} {gsub(/\n$/, ""); printf "%s\n", $1;}' | \
        mawk '{printf "/usr/bin/time bash '${batch_dir}'/respective_scripts/annovar_filtration_per_family.sh %s > '${batch_dir}'/respective_scripts/%s_anno.log 2>&1 &\n", $1, $1;}'
    fi
}



function check_variants_included_vcf() {
    local OPTIND c p e r a v
    while getopts c:p:r:a:v:e:: args
    do 
        case ${args} in
            c) local chromosome=$OPTARG ;; 
            p) local pos=$OPTARG ;;
            e) local end=$OPTARG ;;
            r) local ref=$OPTARG ;;
            a) local alt=$OPTARG ;;
            v) local vcf_file=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local tmp_tsv=${vcf_file/.vcf*/}.$(randomID).bcft.tsv
    trap 'rm ${tmp_tsv}' RETURN

    if [[ ! -z ${end} ]]; then
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%REF\t%ALT[\t%GT]\n' ${vcf_file} | \
        mawk -F '\t' \
        -v qe="${end}" \
        -v qc="${chromosome}" \
        -v qp="${pos}" \
        -v qr="${ref}" \
        -v qa="${alt}" \
        '$1 == qc && $4 == qr && $5 == qa { f = qe - $2; s = $3 - qp; \
                                            if(f >= s) {ratio = s/f;} else {ratio = f/s;} \
                                            if (ratio >= 0.95) {print;}}' > ${tmp_tsv}
    else
        local end=${pos}
        bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT[\t%GT]\n' ${vcf_file} | \
        mawk -F '\t' \
        -v qc="${chromosome}" \
        -v qp="${pos}" \
        -v qr="${ref}" \
        -v qa="${alt}" \
        '$1 == qc && $2 == qp && $4 == qr && $5 == qa{print;}' > ${tmp_tsv}
    fi

    if [[ $(cat ${tmp_tsv} | wc -l) -gt 0 ]]; then
        >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): ${chromosome}:${pos}-${end}:${ref}->${alt} is found in ${vcf_file}."
        if [[ $(cat ${tmp_tsv} | wc -l) -gt 1 ]]; then
            >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): At least 2 variants found match with ${chromosome}:${pos}-${end}:${ref}->${alt} in ${vcf_file}."
            return 0
        else
            return 0
        fi
    else
        >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): ${chromosome}:${pos}-${end}:${ref}->${alt} is not found in ${vcf_file}."
        return 1
    fi
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
    "$@"
fi