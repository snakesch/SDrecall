#!/usr/bin/env bash
self=$(realpath ${BASH_SOURCE[0]})

function log() {
    local msg="$1"
    local script_name="${BASH_SOURCE[1]##*/}"
    local func_name="${FUNCNAME[1]}"
    local line_num="${BASH_LINENO[0]}"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    >&2 echo "[$timestamp] [Script $script_name: Line $line_num] [Func: $func_name] $msg"
}

function silent_remove_tmps {
    local -a targets=($(echo "$@"))

    for target in "${targets[@]}"; do
        rm -rf ${target}
    done
}

function echo_line_no {
    grep -n "$1" $0 | sed "s/echo_line_no//"
}

function randomID {
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}


function display_vcf {
    local input_vcf=${1}
    local head_lines=${2}
    local tmp_tab="/tmp/$(randomID).tsv"

    if [[ -z ${head_lines} ]]; then
        local head_lines=10
    else
        local head_lines=$((head_lines + 1))
    fi

    bcftools view -H ${input_vcf} | tail -n +${head_lines} > ${tmp_tab} && \
    display_table ${tmp_tab} && \
    silent_remove_tmps ${tmp_tab}
}


function display_table {
    # Allow reading from stdin: use '-' as filename or detect piped stdin with no filename
    if [[ "$1" == "-" ]] || { [[ -z "$1" ]] && [ ! -t 0 ]; }; then
        local tmp_tag=$(randomID)
        local stdin_tmp="/tmp/${tmp_tag}.tsv"
        cat - > "${stdin_tmp}"
        display_table "${stdin_tmp}" "$2" "$3"
        silent_remove_tmps "${stdin_tmp}"
        return
    fi
    if [[ -z $2 ]]; then local rows=10; else local rows=${2}; fi
    if [[ -z $3 ]]; then local delimiter="\t"; else local delimiter=${3}; fi

    if [[ ${delimiter} == "\t" ]]; then
        local del_arg=""
    else
        local del_arg="-d ${delimiter}"
    fi

    if [[ ${1} =~ \.vcf$ ]] || \
       [[ ${1} =~ \.vcf\.gz$ ]] || \
       [[ ${1} =~ \.bcf$ ]]; then
        log "Input an VCF file ${1}. Using bcftools to extract the records and view the content of first ${rows} lines."
        display_vcf ${1} ${rows}
        return;
    fi

    ## Another dependency? (Needs to install tsv-utils with conda install bioconda::tsv-utils)
    local row_num=$(tail -n +2 ${1} | wc -l)
    local col_num=$(head -1 ${1} | awk '{print NF;}')

    if [[ ${1} =~ \.gz$ ]]; then
        local tmp_tag=$(randomID)
        zcat ${1} > ${1/.gz/}.${tmp_tag} && \
        local input=${1/.gz/}.${tmp_tag}
    else
        local input=${1}
    fi


    log  "${1} has ${row_num} rows and ${col_num} columns. It looks like:"
    if [[ ${rows} -le 0 ]]; then
        tsv-pretty -u ${del_arg} -m 5000 -l 200 -a ${input} >&2 2> /dev/null
    else
        tsv-pretty -u ${del_arg} -m 5000 -l 200 -a ${input} | \
        head -n ${rows} - >&2 2> /dev/null || >&2 echo ""
    fi

    if [[ ${input} != "${1}" ]]; then
        silent_remove_tmps ${input}
    fi
}


function check_bam_validity {
    local input=${1}
    local expected_lines=${2}

    # log "Start to run check_bam_validity function"
    if [[ -z ${expected_lines} ]]; then local expected_lines=1; fi

    if [[ ! -f ${input} ]] && [[ ! -L ${input} ]] ; then
        log "found ${input} not even exist."
        return 1
    fi
    # log "Now try to check qname format of ${input}."

    local line_num=$( samtools view -c ${input} 2>&1 | awk '{print;}' - )

    if [[ ${line_num} -lt ${expected_lines} ]]; then
        log "found ${input} only has ${line_num} lines where we expect it to have ${expected_lines} lines."
        return 1
    fi

    if check_bam_qname_format ${input}; then
        log "${input} has proper query name format."
    else
        log "${input} has inappropriate query name format."
        return 1
    fi

    check_bam_index ${input} && \
    check_bam_format ${input}
}



function quick_check_bam_validity {
    local input=${1}
    samtools quickcheck -v ${input} && return || return 1
}


function check_bam_format {
    local input=${1}


    if quick_check_bam_validity ${input}; then
        log "${input} format passed quick check by samtools"
        return
    else
        log "${input} format failed quick check by samtools"
        return 1
    fi
}



function check_bam_qname_format {
    local input_bam=${1}
    local threads=${2}

    if [[ -z ${threads} ]]; then local threads=1; fi

    if [[ -L ${input_bam} ]]; then
        local source_bam=$(readlink -f ${input_bam})
        log "${input_bam} is a symbolic link, the source file is ${source_bam}"
    fi

    local issue_lines=$(samtools view ${input_bam} 2>&1 | head -20 | awk -F '\t' '$1 ~ /\/[1-2]$/{print;}' | wc -l)
    if [[ ${issue_lines} -eq 0 ]]; then
        log "No malformed lines found in ${input_bam}."
    else
        log "The bam file ${input_bam} contains problematic query name format: there are hanging /1 or /2 tag at the end of query names."
        return 1
    fi
}


function check_bam_index {
    local bam=${1}
    log "Input BAM files is ${bam}"

    if [[ -L ${bam} ]]; then
        local source_bam=$(readlink -f ${bam})
        log "Input BAM files is a symbolic link: ${bam}"
        if [[ -L ${bam}.bai ]] && [[ ${bam}.bai -nt ${bam} ]]; then
            log "Input BAM files index is a symbolic link: ${bam}.bai and its updated"
        else
            samtools index ${bam} && \
            ln -f -s ${source_bam}.bai ${bam}.bai || { \
            log "Indexing ${bam} failed. Need to sort it first." && \
            samtools sort -O bam -o ${source_bam/.bam/.tmp.bam} ${source_bam} && \
            samtools index ${source_bam} && \
            ln -f -s ${source_bam}.bai ${bam}.bai || { \
            log "Sorting source bam ${source_bam} failed. Quit with error."; \
            return 1; } }
        fi
    elif [[ ! -f ${bam}.bai ]] && [[ ! -f ${bam::-1}i ]]; then
        log "Index file for ${bam} not existed."
        samtools index ${bam} || \
        log "Indexing ${bam} failed. Need to sort it first."
    elif [[ -f ${bam}.bai ]] && [[ ${bam} -nt ${bam}.bai ]]; then
        log "Index file for ${bam} is older than ${bam} itself. Reindexing"
        samtools index ${bam} && ls -lh ${bam}* || \
        log "Indexing ${bam} failed. Need to sort it first."
    elif [[ -f ${bam::-4}.bai ]] && [[ ${bam} -nt ${bam::-4}.bai ]]; then
        log "Index file for ${bam} is older than ${bam} itself. Reindexing"
        samtools index ${bam} || \
        log "Indexing ${bam} failed. Need to sort it first."
    else
        log "Index file for ${bam} is valid."
    fi
}



function check_vcf_validity {
    local input_vcf=${1}
    local expected_lines=${2}
    local expected_samples=${3}

    if [ -z ${expected_lines} ]; then expected_lines=1; fi
    if [ ! -f ${input_vcf} ]; then
        log "${input_vcf} does not even exist."
        return 1
    fi

    local tmp_output=$(mktemp)
    bcftools view -H --threads 1 ${input_vcf} > ${tmp_output} || { \
    log "The input VCF file ${input_vcf} has some format issue detected by bcftools. Quit this function with error."$'\n'; \
    return 1; }

    if [[ $(cat ${tmp_output} | wc -l) -ge ${expected_lines} ]]; then
        log "${input_vcf} has enough records."
    else
        log "${input_vcf} has $(cat ${tmp_output} | wc -l) valid lines while expected to have ${expected_lines} lines."
        return 1
    fi

    silent_remove_tmps ${tmp_output}
    check_vcf_samples ${input_vcf} ${expected_samples}
}


function check_vcf_samples {
    local input_vcf=${1}
    local expected_samples=${2}
    # expected samples should be delimited by comma
    local asks=$(bcftools query -l ${input_vcf} | sort - | awk '{printf "%s\n", $1;}')

    if [[ -z ${expected_samples} ]]; then
        log "No expected samples specified. Quit checking samples"
    else
        local expected_samples=$(echo ${expected_samples} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' | sort - | awk '{printf "%s\n", $1;}')
        if [[ ${expected_samples} == ${asks} ]]; then
            log "Expected samples are identical with actual samples in vcf file ${input_vcf}"
        else
            log "${input_vcf} has these samples:"
            log "${asks}"
            log "${input_vcf} is expected to have these samples:"
            log "${expected_samples}"
            return 1
        fi
    fi
}



function check_gz_file_validity {
    local input=${1}
    if [[ ${input} =~ \.[b]*gz$ ]]; then
        local validity=$(gzip -v -t ${input} 2>&1 | awk 'END{printf $NF;}')
        local size=$(ls -l ${input} | awk '{gsub(/[ ]+/, " "); print;}' | awk '{printf $5}')

        if [[ ${validity} == "OK" ]] && [[ ${size} -gt 30 ]]; then
            log "The GZ format is OK and size is not oddly small."
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
        log "$(timestamp) Input file is not gzipped."
        return 1
    fi
}


function samtools_bam_coverage() {
    # Use getopt to deal the input arguments
    # Only store 3 arguments
    # The input bam
    # The chunksize for splitting
    # The delimiter size (used for bedtools merge -d option)
    local OPTIND i d o b m
    while getopts i:d::o::b::m:: args
    do
        case ${args} in
            i) local input_bam=$OPTARG ;;
            d) local delimiter_size=$OPTARG ;;
            o) local output_bed=$OPTARG ;;
            b) local base_qual=$OPTARG ;;
            m) local map_qual=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -f (input bam) or -c (chunk size)." ;;
        esac
    done

    # Set the default chunk_size to 1Mb and the default delimiter size to 1000bp
    if [[ -z ${delimiter_size} ]]; then
        local delimiter_size=1000
    fi
    if [[ -z ${output_bed} ]]; then
        local output_bed=${input_bam/.bam/.coverage.bed}
    fi
    if [[ -z ${base_qual} ]]; then
        local base_qual=0
    fi
    if [[ -z ${map_qual} ]]; then
        local map_qual=0
    fi


    if [[ -f ${output_bed} ]]; then
        local ori_md5=$(md5sum ${output_bed} | awk '{print $1;}')
    fi

    # Use samtools depth to calculate the coverage for each base in the bam file
    # Use mawk to convert the result into a BED file format
    # Then perform bedtools merge
    # Finally, split the merged BED file into chunks
    samtools depth -q ${base_qual} -Q ${map_qual} ${input_bam} | \
    mawk -F '\t' '{printf "%s\t%s\t%s\n", $1, $2 -1, $2;}' | \
    bedtools merge -d ${delimiter_size} -i - | \
    bedtools sort -i - > ${output_bed}.tmp && \
    local new_md5=$(md5sum ${output_bed}.tmp | awk '{print $1;}')

    if [[ ${ori_md5} != ${new_md5} ]]; then
        mv ${output_bed}.tmp ${output_bed}
    else
        rm ${output_bed}.tmp
    fi
    
    ls -lhtr ${output_bed}
}



function modify_bam_sq_lines () {
    local input_bam=${1}
    local ref_fasta=${2}
    local output_header=${3}

    samtools view -H ${input_bam} | grep -v "@SQ" | grep -v "@PG" > ${output_header} && \
    generate_sq_lines ${ref_fasta} >> ${output_header}
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
    local ref_contig_sizes=${2}
    local rg_index=${3}

    if [[ -z ${rg_index} ]]; then
        local rg_index="PC1"
    fi

    if [[ -z ${ref_contig_sizes} ]]; then
        log "No reference contig sizes provided. Return with error"
        return 1
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
    mawk -v rg_tag="${rg_index}" 'BEGIN{FS=OFS="\t";} \
                                   $0 ~ /^@/ {print;} \
                                   $0 !~ /^@/ {printf "%s:%s", $1, rg_tag; \
                                               for (i=2;i<NF;i++) printf "\t%s", $i;
                                               printf "\t%s\n", $NF;}' - | uniq - | \
    mawk -F '\t' '$0 !~ /^@/ {print;}' - | \
	mawk -F '\t' '$2 < 256 {print;}' -
}


function independent_minimap2_masked () {
    local OPTIND f r s g o p t z m i a k
    while getopts f::r::s::g::o::t::z::m::i::a::k:: args
    do
        case ${args} in
            f) local forward_reads=$OPTARG ;;
            r) local reverse_reads=$OPTARG ;;
            g) local ref_genome=$OPTARG ;;
            a) local masked_genome=$OPTARG ;;
            s) local samp_ID=$OPTARG ;;
            o) local output_align=$OPTARG ;;
            t) local threads=$OPTARG ;;
            z) local ref_contig_sizes=$OPTARG ;;
            m) local mode=$OPTARG ;;
            i) local rg_index=$OPTARG ;;
            k) local kwargs=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    if [[ -z ${ref_genome} ]]; then { echo >&2 "No reference genome provided. Exit. " && exit 1; } fi
    if [[ -z ${samp_ID} ]]; then local samp_ID=$(basename ${forward_reads} | awk '{gsub(/_[a-z]*1\.f[ast]*q[\.gz]*$/, "", $0); printf "%s", $0;}'); fi
    if [[ -z ${output_align} ]]; then local output_align=${samp_ID}.bam; fi
    if [[ -z ${ref_contig_sizes} ]]; then local ref_contig_sizes="${ref_genome}.fai"; fi
    if [[ -z ${threads} ]]; then local threads=1; fi
    if [[ -z ${mode} ]]; then local mode="sr"; fi

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

    log "Minimap2 index file ${masked_genome/.fasta/.mmi} should be generated with minimap2 -x ${mode} -d ${masked_genome/.fasta/.mmi} ${masked_genome}"
    mmi=${masked_genome/.fasta/.mmi}
    tmp=$(mktemp "${mmi}.tmp.XXXXXX")
    minimap2 -x "$mode" -d "$tmp" "$masked_genome" && \
    mv -f "$tmp" "$mmi"

    if [[ ! -f ${ref_genome}.fai ]]; then
        log "Reference genome index file ${ref_genome}.fai not existed. Generate it."
        samtools faidx ${ref_genome} && \
        ls -lh ${ref_genome}.fai
    fi

    log "Running minimap2 --eqx --MD -F 1000 -ax ${mode} --end-bonus 10 --no-end-flt -t ${threads} \
    -R \"@RG\tID:${samp_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}\" \
    ${masked_genome/.fasta/.mmi} ${forward_reads} ${reverse_reads}" && \
    minimap2 -ax ${mode} --eqx --MD -F 1000 --end-bonus 10 --no-end-flt -t ${threads} ${kwargs} -R "@RG\tID:${samp_ID}\tLB:SureSelectXT\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}" \
    ${masked_genome/.fasta/.mmi} ${forward_reads} ${reverse_reads} > ${mid_align} && \
    log "Modify the bam header with modify_bam_sq_lines ${mid_align} ${ref_genome} ${output_align/.bam/.header}" && \
    modify_bam_sq_lines ${mid_align} ${ref_genome} ${output_align/.bam/.header} && \
    log "Modify the masked genome coordinates with modify_masked_genome_coords ${mid_align} ${ref_contig_sizes} ${rg_index}" && \
    modify_masked_genome_coords ${mid_align} ${ref_contig_sizes} ${rg_index} >> ${output_align/.bam/.header} && \
    cat ${output_align/.bam/.header} | \
    samtools view -u -S -@ ${threads} -t ${ref_genome}.fai - | \
    samtools sort -O bam -@ ${threads} -o ${output_align} && \
    samtools index ${output_align} && \
    silent_remove_tmps ${mid_align} ${output_align/.bam/.header} && \
    ls -lh ${output_align}

    if check_bam_validity ${output_align} 1; then
        :
    else
        log "Failed to generate a valid bam file ${output_align} from ${forward_reads} and ${reverse_reads} on ${ref_genome}"
        return 1;
    fi
}


function bcftools_call_per_RG {
    local OPTIND m o c b p
    while getopts m:o::c::b:p:: args
    do
        case ${args} in
            m) local masked_genome=$OPTARG ;;
            b) local masked_bam=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            c) local cpu_threads=$OPTARG ;;
            p) local rg_tag=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -s (priority component bed files) or -t (whole region bed path), -a (input align file), -m (masked genomes paths)." ;;
        esac
    done

    if [[ -z ${cpu_threads} ]]; then local cpu_threads=1; fi

    log "In this iteration, the masked genome path is ${masked_genome}"
    # Skip the record if output_vcf is valid and updated.
    if check_vcf_validity ${output_vcf} 1 && \
       [[ ${output_vcf} -nt ${masked_genome} ]] && \
       [[ ${output_vcf} -nt ${masked_bam} ]]; then
        log "${output_vcf} already valid and updated. Skip the rest"
        return 0
    fi

    export OPENBLAS_NUM_THREADS=${cpu_threads}
    log "Running bcftools mpileup --indels-2.0 --threads ${cpu_threads} -A -a FORMAT/AD,FORMAT/DP -q 10 -Q 15 -f ${masked_genome} ${masked_bam} | bcftools call --threads ${cpu_threads} -mv -P "4e-2" -f GQ -Ou | bcftools norm --threads ${cpu_threads} -m -both -f ${masked_genome} --multi-overlaps 0 -a -Ou - "
    bcftools mpileup --indels-2.0 --threads ${cpu_threads} -A -a FORMAT/AD,FORMAT/DP -q 10 -Q 15 -f ${masked_genome} ${masked_bam} | \
    bcftools call --threads ${cpu_threads} -mv -P "4e-2" -f GQ -Ou | \
    bcftools norm --threads ${cpu_threads} -m -both -f ${masked_genome} --multi-overlaps 0 -a -Ou - | \
    bcftools norm --threads ${cpu_threads} -d exact - | \
    bcftools view --threads ${cpu_threads} -i 'ALT!="*"' -Ov | \
    mawk 'BEGIN{FS=OFS="\t";} {printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, toupper($4), toupper($5); \
                            for(i=6;i<=NF;i++) {printf "\t%s", $i; } \
                            printf "\n"; }' - | \
    bcftools filter --threads ${cpu_threads} -e 'GT != "mis"' -s "${rg_tag}" - | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    display_table ${output_vcf}
    check_vcf_validity ${output_vcf} 1 || \
    { >&2 echo "Generated ${output_vcf} is not valid." && return 1; }
}



function bcftools_concatvcfs {
    local OPTIND v o e s a t c
    while getopts v:o::e::s::a::t::c:: args
    do
        case ${args} in
            v) local input_vcfs=$OPTARG ;;
            o) local merged_vcf=$OPTARG ;;
            e) local ignore_error=$OPTARG ;;
            a) local other_args=$OPTARG ;;
            c) local cpu_threads=$OPTARG ;;
            t) local tmpdir=$OPTARG ;;
            s) local samples=$OPTARG ;;  # Should be delimited by comma
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ ${input_vcfs}  =~ \/ ]] && [[ ! ${input_vcfs} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${input_vcfs} =~ , ]]; then
        local -a vcfs=($(awk '{printf "%s ", $1;}' < ${input_vcfs}))
    else
        local -a vcfs=($(echo ${input_vcfs} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    fi

    if [[ -z ${cpu_threads} ]]; then
        local cpu_threads=1
    fi

    if [[ -z ${tmpdir} ]]; then
        local tmpdir=$TMPDIR
    fi

    log "the vcfs are ${vcfs[*]}"
    log "The merged vcf is ${merged_vcf}"
    # Check file existence.
    if [[ -z ${ignore_error} ]]; then
        local -a invalid_vcfs
        for vcf in "${vcfs[@]}"; do
            if check_vcf_validity ${vcf} 1 2> /dev/null; then
                log "To be merged ${vcf} is valid."
            else
                log "${vcf} not existed or corrupted. Run bash ${central_scripts}/common_bash_utils.sh check_vcf_validity ${vcf} to see for yourself"
                invalid_vcfs+=( ${vcf} )
            fi
            if [[ ${vcf} =~ \.vcf$ ]]; then
                log "Since ${vcf} is plain text format and bcftools -f need to use bgzipped vcfs, we compress the vcf and index it with tabix."
                bcftools sort --temp-dir ${tmpdir} -Oz -o ${vcf}.gz ${vcf} && tabix -f -p vcf ${vcf}.gz && \
                ls -lh ${vcf}.gz
            fi
            if [[ ! -z ${samples} ]]; then
                bcftools view -s "${samples}" -Oz -o ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
                ls -lht ${vcf/.vcf*/.samp.vcf.gz} && \
                check_vcf_validity ${vcf/.vcf*/.samp.vcf.gz} 1 && \
                mv ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
                tabix -f -p vcf ${vcf/.vcf*/.vcf.gz}
            fi
        done

        if [[ ${#invalid_vcfs[@]} -gt 0 ]]; then
            log "The following vcfs are invalid: ${invalid_vcfs[*]}"
            return 1
        fi
    fi

    local tmp_file_list=${tmpdir}/$(randomID).lst
    echo "${vcfs[*]}" | awk -F '\t' 'BEGIN{RS=" ";} length($1) > 0{gsub(/\n$/, ""); if ($1 ~ /\.gz$/) printf "%s\n", $1; else printf "%s.gz\n", $1;}' > ${tmp_file_list}
    log "Here is the temp list file storing the paths of to be concat vcfs:"
    ls -lh ${tmp_file_list}
    cat ${tmp_file_list}

    if [[ ${merged_vcf} =~ \.vcf$ ]]; then
        local plain_merged=${merged_vcf}
        local merged_vcf=${merged_vcf}.gz
    fi

    if [[ $(cat ${tmp_file_list} | wc -l) -eq 1 ]]; then
        cp -f $(cat ${tmp_file_list}) ${merged_vcf}
    else
        # If using file list for input a list of vcfs, each one of them need to be bgzipped and tabix indexed
        bcftools concat -o ${merged_vcf} ${other_args} -a --threads ${cpu_threads} -d exact -Oz -f ${tmp_file_list} && \
        tabix -f -p vcf ${merged_vcf} && \
        bcftools sort --temp-dir ${tmpdir} -o ${merged_vcf/.vcf.gz/.sorted.vcf.gz} -Oz ${merged_vcf} && \
        mv ${merged_vcf/.vcf.gz/.sorted.vcf.gz} ${merged_vcf}
    fi

    if [[ ! -z ${plain_merged} ]]; then
        gunzip -c ${merged_vcf} > ${plain_merged}
    else
        tabix -f -p vcf ${merged_vcf}
    fi

    display_table ${merged_vcf} 20
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    $@
fi
