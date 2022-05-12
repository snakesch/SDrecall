#!/bin/bash

central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
. ${central_scripts}/common_bash_utils.sh

# Create a function to collect read depth from any bam without considering the MAPQ.
function extract_ave_depth_without_MQ(){
    # This function directly output the number of the average depth
    local align_file=${1}
    local target_region=${2}

    if [[ -f ${align_file} ]] && [[ -f ${target_region} ]] && [[ $(wc -l ${target_region} | awk '{print $1;}') -gt 0 ]]; then
        module load samtools 2> /dev/null
        samtools depth -a -b ${target_region} -g DUP,UNMAP,QCFAIL,SECONDARY ${align_file} | awk '{sum+=$3} END{ printf "%s",sum/NR;}' || \
        { >&2 echo "${align_file} at region ${target_region} seems to have no reads covered." && echo "0"; }
        module unload samtools 2> /dev/null
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Not implementing proper args, quit now."
        return 127
    fi
}

# Create a function to collect read depth from any bam without considering the MAPQ.
function extract_ave_depth_with_highqual(){
    # This function directly output the number of the average depth
    local align_file=${1}
    local target_region=${2}

    if [[ -f ${align_file} ]] && [[ -f ${target_region} ]] && [[ $(wc -l ${target_region} | awk '{print $1;}') -gt 0 ]]; then
        module load samtools 2> /dev/null
        samtools view -h ${align_file} | mawk -F '\t' '$0 !~ /XA:Z:/{print;}' | samtools view -b -o ${align_file/.bam/.tmp.bam} && \
        samtools index ${align_file/.bam/.tmp.bam} 2> /dev/null && \
        samtools depth -a -b ${target_region} -Q 30 -g DUP,UNMAP,QCFAIL,SECONDARY ${align_file/.bam/.tmp.bam} | awk '{sum+=$3} END{ printf "%s",sum/NR;}' && \
        rm -f ${align_file/.bam/.tmp.bam} 2> /dev/null
        module unload samtools 2> /dev/null
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Not implementing proper args, quit now."
        return 127
    fi
}


function remove_entire_contigs_with_all_N(){
    local OPTIND f o b t
    while getopts f:o::b:t:: args
    do
        case ${args} in
            f) local fasta_file=$OPTARG ;;
            o) local output_fasta=$OPTARG ;;
            b) local bed_file=$OPTARG ;;
            t) local tmp_dir=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -f fasta file and -b ${bed_file}." ;;
        esac
    done

    if [[ -z ${tmp_dir} ]]; then
        ls /staging 2> /dev/null && local tmp_dir=/staging/test_tmp/$(randomID) || \
        local tmp_dir=/paedyl01/disk1/yangyxt/test_tmp/$(randomID)
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    elif [[ ! -d ${tmp_dir} ]]; then
        mkdir -p ${tmp_dir} && chmod -R ug+rw ${tmp_dir}
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input fasta is $(ls -lh ${fasta_file}) And it looks like:"
    head -5 ${fasta_file}

    if [[ -z ${output_fasta} ]]; then local output_fasta="${fasta_file}.tmp"; fi
    if check_masked_genome ${fasta_file} ${bed_file}; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input fasta is $(ls -lh ${fasta_file}) And it's been already removed redundant all N contigs."
    else
        local fasta_name=$(basename ${fasta_file})
        mawk 'BEGIN {currc = 0;}{
                if ($1 ~ /^>/) { start = FNR + 1; alln = 1; currc = $0; }
                else if ($0 !~ /^[N]+$/) { arr[FNR] = $0"\n"; alln = 0; start = FNR; }
                else { start = FNR; arr[FNR] = $0"\n"; alln = 1;}
                while (FNR > 0) {
                    ret = getline;
                    if (ret <= 0 ) { end = FNR; break; }
                    else if ($0 ~ /^>/) { end = FNR - 1; nextc = $0; break; }
                    else if ($0 !~ /^[N]+$/) { arr[FNR] = $0"\n"; alln *= 0; }
                    else { arr[FNR] = $0"\n"; alln *= 1; }
                    }
                if ( alln == 0 ) { print currc; for (i=start;i<=end;i++) { printf "%s", arr[i];}}
                delete arr; currc = nextc; }' < ${fasta_file} > ${tmp_dir}/${fasta_name} && \
        mv ${tmp_dir}/${fasta_name} ${output_fasta} && \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Output fasta is $(ls -lh ${output_fasta}) And it looks like: "
        head ${output_fasta}
        if [[ "${output_fasta}" == "${fasta_file}.tmp" ]]; then mv ${output_fasta} ${fasta_file}; fi
    fi

    trap "rm -rf ${tmp_dir}" RETURN
}


# Create a function to get the complementary region of a bed file
function complementing_bed(){
    local OPTIND i c w
    while getopts i:c::w:: args
    do 
        case ${args} in
            i) local input_bed=$OPTARG ;;
            w) local whole_region=$OPTARG ;;
			c) local complement_bed=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${whole_region} ]]; then local whole_region=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.bed; fi
    if [[ -z ${complement_bed} ]]; then
        local complement_bed=$(dirname ${input_bed})/$(basename ${whole_region} | awk -F '.' '{printf $1"."$2;}').sub.$(basename ${input_bed} | awk -F '.' '{printf $1;}').bed
    fi 

    module load BEDTools
    bedtools subtract -a ${whole_region} -b ${input_bed} > ${complement_bed}
    module unload BEDTools
}


# Create a function to mask the ref genome at the interval specified by a bed file
function mask_genome(){
    local OPTIND r b o
    while getopts r::b:o:: args
    do
        case ${args} in
            r) local ref_genome=$OPTARG ;;
            b) local bed_file=$OPTARG ;;
            o) local output_genome=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${output_genome} ]]; then local output_genome=$(dirname ${bed_file})/$(basename ${ref_genome} | awk -F '.' '{printf "%s.%s", $1, $2;}').$(basename ${bed_file} | awk -F '.' '{printf $1;}').masked.fasta; fi


    local tmp_output_genome=${output_genome::-6}.$(randomID).fasta 
    module load seqtk
    seqtk seq -l60 -M ${bed_file} -n N ${ref_genome} > ${tmp_output_genome} && \
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Masked fasta is $(ls -lh ${output_genome})"
    module unload seqtk

    remove_entire_contigs_with_all_N -f ${tmp_output_genome} -o ${output_genome} -b ${bed_file}
    trap "rm -f ${tmp_output_genome}" RETURN
}


# Create a function to call GVCF file with multiple diploidy
function gvcf_calling_polyploidy(){
    local OPTIND r p i g b s
    while   getopts r::p:i:b::g::s: args
    do
        case ${args} in
            r) local ref_genome=$OPTARG ;;
            p) local ploidy=$OPTARG ;;
            i) local input_bam=$OPTARG ;;
            s) local special_interval=$OPTARG ;;
            b) local output_bam=$OPTARG ;;
            g) local output_gvcf=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac
    done

    local interval_file_name=$(basename ${special_interval} | awk -F '.' '{for (i=1;i<NF-1;i++) printf "%s.", $i; printf $(NF-1);}')
    local input_bam_ID=$(basename ${input_bam} | awk -F '.' '{printf $1;}')
    local bamf=$(dirname ${input_bam})
    local gvcf=$(echo ${bamf} | awk -F '/' '{for (i=1;i<NF;i++) printf "%s/", $i;}')/gvcfs

    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${output_gvcf} ]]; then 
        local output_gvcf=${gvcf}/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz
    elif [[ ! ${output_gvcf} =~ (/|\.) ]]; then
        local __output_gvcf=${output_gvcf}
        local output_gvcf=${gvcf}/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz
    fi
    if [[ -z ${output_bam} ]]; then
        local output_bam=${bamf}/${input_bam_ID}.only_${interval_file_name}.realigned.bam
    elif [[ ! ${output_bam} =~ (/|\.) ]]; then
        local __output_bam=${output_bam}
        local output_bam=${bamf}/${input_bam_ID}.only_${interval_file_name}.realigned.bam
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Trying to call GVCFs in region specified by ${interval_file_name}, with BAM file ${input_bam}, ploidy ${ploidy}. Output to ${output_gvcf}"
    
    if [[ -z ${gatk} ]]; then local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk; fi

    time $gatk --java-options "-Xmx75G" HaplotypeCaller \
    --emit-ref-confidence GVCF \
    -R ${ref_genome} \
    --ploidy ${ploidy} \
    -I ${input_bam} \
    -O ${output_gvcf} \
    -bamout ${output_bam} \
    --force-active \
    --allow-non-unique-kmers-in-ref \
    --debug-assembly \
    --kmer-size 21 \
    --max-num-haplotypes-in-population 256 \
    --bam-writer-type CALLED_HAPLOTYPES \
    --linked-de-bruijn-graph \
    -L ${special_interval} \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -G StandardHCAnnotation \
    -imr OVERLAPPING_ONLY \
    -OBI true \
    --verbosity WARNING

    if [[ ! -z ${__output_gvcf} ]]; then eval $__output_gvcf="'${output_gvcf}'"; fi
    if [[ ! -z ${__output_bam} ]]; then eval $__output_bam="'${output_bam}'"; fi
}


function check_masked_genome() {
    local masked_genome=${1}
    local target_region=${2}

    if [[ ! -f ${masked_genome} ]]; then false; fi

    module load BEDTools/2.27.1
    local validate=$(bedtools getfasta -fi ${masked_genome} -bed ${target_region} | \
    mawk -F '\t' '$1 !~ /^>/{printf "%s", $0;}' - | \
    mawk '$0 ~ /N/{printf "False";} $0 !~ /N/{printf "True";}' - )
    module unload BEDTools/2.27.1

    if [[ ${validate} == "True" ]]; then
        return 0
    elif [[ ${validate} == "False" ]]; then
        return 1
    fi
}


function validate_genome(){
    local tar_fasta=${1}
    local expected_lines=${2}
    if [[ -z ${expected_lines} ]]; then local expected_lines=100; fi 

    if [[ -f ${tar_fasta} ]] && [[ $(wc -l ${masked_genome} | awk '{printf $1;}') -ge ${expected_lines} ]]; then
        build_index_for_genome ${tar_fasta} && \
        true
    else
        false
    fi
}


function prepare_masked_genome(){
    local OPTIND g t s m r
    while   getopts g::t:s:m::r:: args
    do
        case ${args} in
            g) local ref_genome=$OPTARG ;;
            s) local priority_components=$OPTARG ;;
            t) local whole_target_region=$OPTARG ;;
            m) local __masked_genomes_varname=$OPTARG ;;
            r) local replace=$OPTARG;; 
            *) echo "No argument passed. Pls at least specify -s (priority_components) or -t (whole target region bed file path)." ;;
        esac
    done

    # if [[ -z ${replace} ]]; then local replace="True"; fi

    local -a priority_component_arr=($(echo ${priority_components} | awk 'BEGIN{RS=",";} {printf "%s ",$1;}'))
    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${central_scripts} ]]; then local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts; fi

    local region_name=$(basename ${whole_target_region} | awk -F '.' '{printf $1;}')
    local -a masked_genomes

    for target_bed in "${priority_component_arr[@]}"; do
        local interval_file_name=$(basename ${target_bed} | awk -F '.' '{for (i=1;i<NF-1;i++) printf "%s.", $i; printf $(NF-1);}')
        
        # First create a complementary bed file.
        local complement_bed=$(dirname ${target_bed})/$(basename ${ref_genome} | awk -F '.' '{printf $1"."$2;}').sub.$(basename ${target_bed} | awk -F '.' '{printf $1;}').bed
        complementing_bed -i ${target_bed} -c ${complement_bed}
        
        # Then create a masked genome.
        local masked_genome=$(dirname ${complement_bed})/$(basename ${ref_genome} | awk -F '.' '{printf "%s.%s", $1, $2;}').$(basename ${complement_bed} | awk -F '.' '{for (i=1;i<NF;i++) printf "%s.", $i;}')masked.fasta
        if [[ ${replace} == "True" ]]; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": We set up to replace the current masked genome."
            mask_genome -b ${complement_bed} -r ${ref_genome} -o ${masked_genome} && \
            build_index_for_genome ${masked_genome}
        elif validate_genome ${masked_genome} && check_masked_genome ${masked_genome} ${target_bed}; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": We already have a masked genome: ${masked_genome} and it has $(cat ${masked_genome} | wc -l) lines"
        elif validate_genome ${masked_genome}; then
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The masked genome may not have been removed contigs with all N."
            remove_entire_contigs_with_all_N -f ${masked_genome} -b ${target_bed}
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Masked genome $(ls -lh ${masked_genome}) is corrupted or not existed. Need to regenerate it."
            mask_genome -b ${complement_bed} -r ${ref_genome} -o ${masked_genome} && \
            build_index_for_genome ${masked_genome}
        fi
        check_return_code
        masked_genomes+=( "${masked_genome}" )
    done

    local masked_genome_list=$(echo "${masked_genomes[@]}" | awk 'BEGIN{RS=" ";} { gsub(/\n$/, ""); printf "%s,", $1;}' | awk '{gsub(/,$/, ""); print; }')
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Masked genomes are ${masked_genome_list}. Output masked genome var name is $__masked_genomes_varname"
    if [[ ${#__masked_genomes_varname} -gt 0 ]]; then eval $__masked_genomes_varname="'${masked_genome_list}'"; fi
}


function prepare_fastq_files(){
    local OPTIND t a f r
    while   getopts t:a:f::r:: args
    do
        case ${args} in
            t) local whole_target_region=$OPTARG ;;
            a) local input_align_file=$OPTARG ;;
            f) local __forward_read=$OPTARG ;;
            r) local __reverse_read=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -s (priority component bed files) or -t (whole region bed path), -a (input align file), -m (masked genomes paths)." ;;
        esac
    done

    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${central_scripts} ]]; then local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts; fi

    local sample_ID=$(basename ${input_align_file} | awk -F '.' '{printf $1;}')
    local region_name=$(basename ${whole_target_region} | awk -F '.' '{printf $1;}')

    local extract_fastq=$(dirname ${input_align_file})/${sample_ID}.${region_name}-XA.fastq

    if [[ -f ${__forward_read} ]] && [[ -f ${__reverse_read} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input forward reads file and reverse reads file exist."
        ls -lh ${__forward_read} && \
        ls -lh ${__reverse_read} && \
        bash ${central_scripts}/select_bam_by_regions.sh \
        -i ${input_align_file} \
        -b ${whole_target_region} \
        -o ${extract_fastq} \
        -f ${__forward_read} \
        -r ${__reverse_read} \
        -x yes
    else
        bash ${central_scripts}/select_bam_by_regions.sh \
        -i ${input_align_file} \
        -b ${whole_target_region} \
        -o ${extract_fastq} \
        -x yes
    fi

    check_return_code

    local ext_forward_read=$(gnws ${extract_fastq})_1.fastq
    local ext_reverse_read=$(gnws ${extract_fastq})_2.fastq

    ls -lh ${ext_forward_read}
    ls -lh ${ext_reverse_read}

    if [[ ! -f ${__forward_read} ]]; then eval $__forward_read="'${ext_forward_read}'"; fi
    if [[ ! -f ${__reverse_read} ]]; then eval $__reverse_read="'${ext_reverse_read}'"; fi
}


function call_polyploidy_per_priority_region(){
    local OPTIND g t a s f r m
    while getopts g::t:a:s:f::r::m: args
    do
        case ${args} in
            g) local ref_genome=$OPTARG ;;
            s) local priority_component=$OPTARG ;;
            t) local whole_target_region=$OPTARG ;;
            a) local input_align_file=$OPTARG ;;
            f) local forward_read=$OPTARG ;;
            r) local reverse_read=$OPTARG ;;
            m) local masked_genome=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -s (priority component bed files) or -t (whole region bed path), -a (input align file), -m (masked genomes paths)." ;;
        esac
    done

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input forward read is $(ls -lh ${forward_read})"
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Input reverse read is $(ls -lh ${reverse_read})"
    
    if [[ -z ${ref_genome} ]]; then local ref_genome=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ -z ${central_scripts} ]]; then local central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts; fi

    local sample_ID=$(basename ${input_align_file} | awk -F '.' '{printf $1;}')
    local region_name=$(basename ${whole_target_region} | awk -F '.' '{printf $1;}')

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": In this iteration, the masked genome path is ${masked_genome}"
    local input_bam_ID=$(basename ${input_align_file} | awk -F '.' '{printf $1;}')
    local interval_file_name=$(basename ${priority_component} | awk -F '.' '{for (i=1;i<NF-1;i++) printf "%s.", $i; printf $(NF-1);}')

    # Calculate the raw depth in raw bam file at this area.
    local ori_depth=$(extract_ave_depth_without_MQ ${input_align_file} ${priority_component})
    local ori_highqual_depth=$(extract_ave_depth_with_highqual ${input_align_file} ${priority_component})

    >&2 echo $(timestamp)": The original average depth in ${input_align_file} across interval ${priority_component} is ${ori_depth}(without considerinng MQ)."
    
    # Align generated reads to masked genome.
    local masked_align=$(gnws ${input_align_file}).only_$(get_name ${priority_component}).bam
    if [[ -f ${masked_genome}.pac ]] && [[ ${masked_genome}.pac -nt ${masked_genome} ]] && \
    [[ -f ${masked_genome}.bwt ]] && [[ ${masked_genome}.bwt -nt ${masked_genome} ]] && \
    [[ -f ${masked_genome}.ann ]] && [[ ${masked_genome}.ann -nt ${masked_genome} ]] && \
    [[ -f ${masked_genome}.amb ]] && [[ ${masked_genome}.amb -nt ${masked_genome} ]] && \
    [[ -f ${masked_genome}.sa ]] && [[ ${masked_genome}.sa -nt ${masked_genome} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Seems that index files are built by bwa index -a bwtsw command:"
        ls -lh ${masked_genome}.*
        independent_bwa_alignment -g ${masked_genome} -s ${sample_ID} -f ${forward_read} -r ${reverse_read} -o ${masked_align} && ls -lh ${masked_align} || { \
        independent_bwamem2_alignment -g ${masked_genome} -s ${sample_ID} -f ${forward_read} -r ${reverse_read} -o ${masked_align} && ls -lh ${masked_align} || \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Mapping failed with both bwa-mem and bwa-mem2."; }
    elif [[ -f ${masked_genome}.bwt.2bit.64 ]] && [[ -f ${masked_genome}.0123 ]] && [[ ${masked_genome}.0123 -nt ${masked_genome} ]] && [[ ${masked_genome}.bwt.2bit.64 -nt ${masked_genome} ]]; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Seems the index files are built by bwamem2: $(ls -lh ${masked_genome}.*) And we try to align reads using bwa-mem2" && \
        independent_bwamem2_alignment -g ${masked_genome} -s ${sample_ID} -f ${forward_read} -r ${reverse_read} -o ${masked_align} && ls -lh ${masked_align} || \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Mapping failed with bwa-mem2."
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Seems that index files are not ready."
        ls -lh ${masked_genome}.*
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Now we try to index that genome ${masked_genome}."
        build_index_for_genome ${masked_genome}
    fi

    # Calculate the current depth in bam aligned to masked genome.
    local multi_aligned_depth=$(extract_ave_depth_without_MQ ${masked_align} ${priority_component})
    local post_depth=$(awk -v md="${multi_aligned_depth}" -v od="${ori_highqual_depth}" 'BEGIN{printf "%s", md + od;}')

    # Calculate ploidy
    local ploidy=$(awk -v od="${ori_depth}" -v pd="${post_depth}" 'BEGIN{pl=2*pd/od; printf "%.0f", pl; }')
    >&2 echo $(timestamp)": The depth in this region specified by ${priority_component} get increased by "$(awk -v od="${ori_depth}" -v pd="${post_depth}" 'BEGIN{pl=pd/od; printf "%.1f", pl;}')" fold."
    >&2 echo $(timestamp)": After multiplying increased fold with diploidy, we get ploidy of "$(awk -v od="${ori_depth}" -v pd="${post_depth}" 'BEGIN{pl=2*pd/od; printf "%.1f", pl;}')
    echo "$(timestamp): We round that ploidy number to integer "$(awk -v od="${ori_depth}" -v pd="${post_depth}" 'BEGIN{pl=2*pd/od; printf "%.0f", pl;}')

    if [[ ${ploidy} -ge 2 ]]; then
        # Generate gvcf file with polyploidy
        gvcf_calling_polyploidy \
        -p ${ploidy} \
        -i ${masked_align} \
        -s ${priority_component} \
        -b masked_reassembly_align \
        -g $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz
        
        local gvcf_hom="$(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz"
        if [[ $(zcat ${gvcf_hom} | awk '$1 !~ /^#/{print;}' | wc -l | awk '{print $1;}') -eq 0 ]]; then
            touch $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz && \
            touch $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.dip.vcf.gz && \
            return 0
        fi

        if [[ -z ${gatk} ]]; then local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk; fi

        ${gatk} GenotypeGVCFs \
        -R /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
        -V $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz \
        -O $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz \
        -G StandardAnnotation \
        -G AS_StandardAnnotation && \
        check_vcf_validity $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz 1 || \
        { rm -f $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz && \
        >&2 echo "Generated $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz is not valid. Deleted." && return 1; }

        local vcf_hom=$(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz
        if [[ $(zcat ${vcf_hom} | awk '$1 !~ /^#/{print;}' | wc -l | awk '{print $1;}') -eq 0 ]]; then
            >&2 echo "$(timestamp): In function ${FUNCNAME}: Generated ${vcf_hom} contains 0 records. Abandon this interval region."
            # touch $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.dip.vcf.gz
            return 1
        fi

        bash ${central_scripts}/deal_with_multi_allelic_records_at_vcf_level.sh $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}
        bash ${central_scripts}/common_bash_utils.sh replace_contig_lines_in_vcfhead $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz

        python3 ${central_scripts}/convert_vcf_to_diploidy.py -vp $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz || \
        { >&2 echo  "$(timestamp): In function ${FUNCNAME}: Failed to convert vcf to diploid vcf. Quit with error." && return 1; }
        local poor_cov_bed=${input_align_file/.bam/}.${interval_file_name}.txt
        pick_poorcov_region_bam -i ${input_align_file} -t ${priority_component} -o ${poor_cov_bed} -d 15 && \
        >&2 echo "$(timestamp): In function ${FUNCNAME}: The poor coverage region in ${input_bam} overlapping with ${priority_component} is stored in ${poor_cov_bed}: (Show first 5 rows)"
        ls -lh ${poor_cov_bed} && \
        head -5 ${poor_cov_bed}
        bcftools view -R ${poor_cov_bed} -Oz -o $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz
        if check_vcf_validity $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz; then
            >&2 echo "$(timestamp): In function ${FUNCNAME}: Pickout variants in poor coverage regions to store the variants in $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz"
            ls -lh $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz
        else
            >&2 echo "$(timestamp): In function ${FUNCNAME}: Pickout variants in poor coverage regions, the result vcf seems problematic."
            ls -lh $(dirname ${input_align_file})/${input_bam_ID}.only_${interval_file_name}.vcf.gz
            return 1
        fi
    fi
}

function main_call_polyploidy(){
    local OPTIND g t a s f r m
    while getopts g::t:a:s:f::r::m: args
    do
        case ${args} in
            g) local ref_genome=$OPTARG ;;
            s) local priority_components=$OPTARG ;;
            t) local whole_target_region=$OPTARG ;;
            a) local input_align_file=$OPTARG ;;
            f) local forward_read=$OPTARG ;;
            r) local reverse_read=$OPTARG ;;
            m) local masked_genomes=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -s (priority component bed files) or -t (whole region bed path), -a (input align file), -m (masked genomes paths)." ;;
        esac
    done

    prepare_fastq_files -t ${whole_target_region} -a ${input_align_file} -f ${forward_read} -r ${reverse_read}

    local sample_ID=$(basename ${input_align_file} | awk -F '.' '{printf $1;}')
    local region_name=$(basename ${whole_target_region} | awk -F '.' '{printf $1;}')
    local extract_fastq=$(dirname ${input_align_file})/${sample_ID}.${region_name}-XA.fastq
    local ext_forward_read=$(gnws ${extract_fastq})_1.fastq
    local ext_reverse_read=$(gnws ${extract_fastq})_2.fastq

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Extracted forward read is $(ls -lh ${ext_forward_read})"
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Extracted reverse read is $(ls -lh ${ext_reverse_read})"

    local -a priority_component_arr=($(echo ${priority_components} | awk 'BEGIN{RS=",";} {printf "%s ",$1;}'))
    local -a masked_genome_arr=($(echo ${masked_genomes} | awk 'BEGIN{RS=",";} {printf "%s ",$1;}'))

    local -A region_genome_map
    
    for i in "${!priority_component_arr[@]}"; do
        region_genome_map["${priority_component_arr[$i]}"]="${masked_genome_arr[$i]}"
    done

    for target_bed in "${priority_component_arr[@]}"; do
        local masked_genome=${region_genome_map["${target_bed}"]}

        call_polyploidy_per_priority_region \
        -s ${target_bed} \
        -t ${whole_target_region} \
        -a ${input_align_file} \
        -f ${ext_forward_read} \
        -r ${ext_reverse_read} \
        -m ${masked_genome}
    done
}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    # main_call_polyploidy -s $(pwd)/NCF1.bed -t $(pwd)/NCF1_NCF1C_NCF1B.geneonnly.bed -a $(pwd)/WGS_simu_chr7.bam
    # main_call_polyploidy "$@"
    "$@"
fi