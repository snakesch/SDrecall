function check_bam_validity {
    local input=${1}
    local expected_lines=${2}

    # >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Start to run check_bam_validity function"
    if [[ -z ${expected_lines} ]]; then local expected_lines=1; fi

    if [[ ! -f ${input} ]] && [[ ! -L ${input} ]] ; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, found ${input} not even exist."
        return 1
    fi
    # >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Now try to check qname format of ${input}."

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

    check_bam_index ${input} && \
    check_bam_format ${input}
}

function check_bam_format {
    local input=${1}
    local gatk="bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh gatk_wrapper"

    if quick_check_bam_validity ${input}; then
        >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} format passed quick check by samtools"
        return;
    fi

    local tmp_check_log=${input::-4}.$(randomID).check
    >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} format check log is ${tmp_check_log}"
    ${gatk} --java-options "-Xmx20G" \
    ValidateSamFile \
    -I ${input} \
    --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE \
    -MO 1000000 > ${tmp_check_log} 2>&1 || \
    >&2 echo "In function ${FUNCNAME}: $(timestamp): Seems some error format in ${input}, check log is ${tmp_check_log}"

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
                >&2 mawk '$1 ~ /^ERROR::INVALID_PLATFORM_VALUE/{print $0;}' ${tmp_check_log}
            elif [[ ${err} =~ "MISSING_READ_GROUP" ]]; then
                >&2 echo "In function ${FUNCNAME}: $(timestamp): ${input} does not have the read group line. Just ignore."
            fi
        done
        # IF the for loop above is executed without return 1, then no other error types are identified. return
    fi
    rm ${tmp_check_log} 2> /dev/null || >&2 echo "${tmp_check_log} already deleted"
}

function check_bam_index {
    local bam=${1}
    >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): Input BAM files is ${bam}"

    if [[ -L ${bam} ]]; then
        local source_bam=$(readlink -f ${bam})
        >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): Input BAM files is a symbolic link: ${bam}"
        if [[ -L ${bam}.bai ]] && [[ ${bam}.bai -nt ${bam} ]]; then
            >&2 echo "Line ${LINENO}: In function ${FUNCNAME}: $(timestamp): Input BAM files index is a symbolic link: ${bam}.bai and its updated"
        else
            samtools index ${bam} && \
            ln -f -s ${source_bam}.bai ${bam}.bai || { \
            >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Indexing ${bam} failed. Need to sort it first." && \
            samtools sort -O bam -o ${source_bam/.bam/.tmp.bam} ${source_bam} && \
            samtools index ${source_bam} && \
            ln -f -s ${source_bam}.bai ${bam}.bai || { \
            >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Sorting source bam ${source_bam} failed. Quit with error."; \
            return 1; } }
        fi
    elif [[ ! -f ${bam}.bai ]] && [[ ! -f ${bam::-1}i ]]; then
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

function modify_bam_sq_lines () {
    local input_bam=${1}
    local ref_fasta=${2}
    local output_header=${3}

    local sq_lines=$()

    samtools view -H ${input_bam} | grep -v "@SQ" | grep -v "@PG" > ${output_header} && \
    generate_sq_lines ${ref_fasta} >> ${output_header} && \
    ls -lht ${output_header}
}

function generate_sq_lines () {
    local ref_fasta=${1}

    if ([[ -f ${ref_fasta}.fai ]] && [[ ${ref_fasta}.fai -nt ${ref_fasta} ]]) || \
       ([[ -f ${ref_fasta/.fasta/.fai} ]] && [[ ${ref_fasta/.fasta/.fai} -nt ${ref_fasta} ]]); then
        :
    else
        samtools faidx ${ref_fasta}
    fi

    awk 'BEGIN{OFS="\t"} {printf "@SQ\tSN:%s\tLN:%s\n", $1, $2;}' ${ref_fasta}.fai
}

function modify_masked_genome_coords () {
     local mid_align=${1}
     local pc_index=${3}
     local ref_contig_sizes=${2}

     if [[ -z ${pc_index} ]]; then
         local pc_index="PC1"
     fi

     if [[ -z ${ref_contig_sizes} ]]; then
         local ref_contig_sizes="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.contigsize.genome"
     fi

     mawk 'BEGIN{FS=OFS="\t"; sq_line = 0;} \
           $0 ~ /^@SQ/ && sq_line < 1{print; sq_line+=1;} \
           $0 !~ /^@SQ/{print;}' ${mid_align} | \
     mawk 'BEGIN{FS=OFS="\t";} \
           FNR == NR {size_arr[$1]=$2;} \
           NR > FNR && $0 !~ /^@/{ split($3, arr, ":"); \
                                   split($7, mate_arr, ":"); \
                                   if ($7 == "=") { \
                                     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, arr[1], $4 + arr[2], $5, $6, $7, $8 + arr[2], $9, $10, $11; } \
                                   else { \
                                     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, arr[1], $4 + arr[2], $5, $6, mate_arr[1], $8 + mate_arr[2], $9, $10, $11; } \
                                   for (i=12; i<NF; i++) printf "%s\t", $i; \
                                   printf "%s\n", $NF;} \
           NR > FNR && $0 ~ /^@SQ/{for (i in size_arr) { \
                                   if (i ~ /_/) { printf "@SQ\tSN:%s\tLN:%s\tAH:*\n", i, size_arr[i]; } \
                                   else { printf "@SQ\tSN:%s\tLN:%s\n", i, size_arr[i]; }}} \
           NR > FNR && $0 !~ /^@SQ/ && $0 ~ /^@/{print;}' ${ref_contig_sizes} - | \
     mawk -v pc_tag="${pc_index}" 'BEGIN{FS=OFS="\t";} \
                                   $0 ~ /^@/ {print;} \
                                   $0 !~ /^@/ {printf "%s:%s", $1, pc_tag; \
                                               for (i=2;i<NF;i++) printf "\t%s", $i;
                                               printf "\t%s\n", $NF;}' - | uniq -
 }


function independent_minimap2_masked () {
     local OPTIND f r s g o p t e z m i b a c
     while getopts f::r::p::s::g::o::t::e::z::m::i::b::a::c:: args
     do
         case ${args} in
             f) local forward_reads=$OPTARG ;;
             r) local reverse_reads=$OPTARG ;;
             p) local merged_pair_reads=$OPTARG ;;
             g) local ref_genome=$OPTARG ;;
             a) local masked_genome=$OPTARG ;;
             s) local samp_ID=$OPTARG ;;
             e) local expected_lines=$OPTARG ;;
             o) local output_align=$OPTARG ;;
             t) local threads=$OPTARG ;;
             z) local ref_contig_sizes=$OPTARG ;;
             m) local mode=$OPTARG ;;
             i) local pc_index=$OPTARG ;;
             b) local raw_align=$OPTARG ;;
             c) local nm_cutoff=$OPTARG ;;
             *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
         esac
     done

     if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
     if [[ -z ${samp_ID} ]]; then local samp_ID=$(basename ${forward_reads} | awk '{gsub(/_[a-z]*1\.f[ast]*q[\.gz]*$/, "", $0); printf "%s", $0;}'); fi
     if [[ -z ${output_align} ]]; then local output_align=${samp_ID}.bam; fi
     if [[ -z ${expected_lines} ]]; then local expected_lines=1; fi
     if [[ -z ${ref_contig_sizes} ]]; then local ref_contig_sizes="${ref_genome}.fai"; fi
     if [[ -z ${threads} ]]; then local threads=1; fi
     if [[ -z ${mode} ]]; then local mode="sr"; fi
     if [[ -z ${nm_cutoff} ]]; then local nm_cutoff=5; fi

     if [[ ${forward_reads} =~ \.gz$ ]]; then
         # Since novoalign free version does not support the gz fastq input. we need to prepare plain version fastq
         gunzip -c ${forward_reads} > ${forward_reads/.gz/} && \
         local forward_reads=${forward_reads/.gz/}
     fi

     if [[ ${reverse_reads} =~ \.gz$ ]]; then
         gunzip -c ${reverse_reads} > ${reverse_reads/.gz/} && \
         local reverse_reads=${reverse_reads/.gz/}
     fi

     local mid_align=${output_align/.bam/.sam}

     if [[ ${masked_genome/.fasta/.mmi} -ot ${masked_genome} ]]; then
         minimap2 -x ${mode} -d ${masked_genome/.fasta/.mmi} ${masked_genome}
     fi

     log "Running minimap2 --eqx --MD -F 1000 -ax ${mode} -t ${threads} \
     -R \"@RG\tID:${samp_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}\" \
     ${masked_genome/.fasta/.mmi} ${forward_reads} ${reverse_reads}" && \
     minimap2 -ax ${mode} --eqx --MD -F 1000 -t ${threads} -R "@RG\tID:${samp_ID}\tLB:SureSelectXT\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}" \
     ${masked_genome/.fasta/.mmi} ${forward_reads} ${reverse_reads} > ${mid_align} && \
     modify_masked_genome_coords ${mid_align} ${ref_contig_sizes} ${pc_index} | \
     samtools view -S -b - | samtools sort -O bam -o ${output_align} && \
     samtools index ${output_align}

     if check_bam_validity ${output_align} 1; then
         :
     else
         log "Failed to generate a valid bam file ${output_align} from ${forward_reads} and ${reverse_reads} on ${ref_genome}"
         return 1;
     fi
 }

