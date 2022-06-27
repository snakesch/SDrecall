#!/home/yangyxt/anaconda3/bin/bash
# This script is designed to process one bam file only
central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
gene_anno_dir=/paedyl01/disk1/yangyxt/public_data/gene_annotation
gatk=/home/yangyxt/software/gatk-4.2.6.1/gatk

module load java
module load samtools

. ${central_scripts}/common_bash_utils.sh

function select_bam {
    local OPTIND g c b i o a t
    while getopts g::c::b::i::o::a::t:: args
    do 
        case ${args} in
            g) local genelist=$OPTARG ;;
            c) local coordinates=$OPTARG ;;
            b) local region_bed=$OPTARG ;;
            i) local input_bam=$OPTARG ;;
            o) local output_bam=$OPTARG ;;
            a) local assembly=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bam sample path." ;;
        esac
    done

    if [[ -z ${genelist} ]]; then
        if [[ -z ${coordinates} ]]; then
            if [[ -z ${region_bed} ]]; then
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Must specify target regions either by gene names or coordinates. Does not receive any of the two info, exiting now."
                exit 1
            else
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region specified by the region based bed file, the first three columns must be coordinates"
                display_table ${region_bed}
                local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).txt
                samtools view -L ${region_bed} ${input_bam} | awk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
                ${gatk} FilterSamReads -I ${input_bam} -O ${output_bam} --FILTER includeReadList -RLF ${tmp_read_ref_names} -SO coordinate && \
                rm ${tmp_read_ref_names}
            fi
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region sepcified by the coordinates"
            local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).txt
            samtools view ${input_bam} "${coordinates}" | awk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
            ${gatk} FilterSamReads -I ${input_bam} -O ${output_bam} --FILTER includeReadList -RLF ${tmp_read_ref_names} -SO coordinate && \
            rm ${tmp_read_ref_names}
        fi
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region sepcified by the genelist, separated by comma "${genelist}
        local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).txt
        local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed
        bash ${central_scripts}/generate_bed_by_gene_names.sh -g ${genelist} -o ${tmp_bed}
        display_table ${tmp_bed}
        samtools view -L ${tmp_bed} ${input_bam} | awk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names} && \
        display_table ${tmp_read_ref_names} && \
        rm ${tmp_bed}
        ${gatk} FilterSamReads -I ${input_bam} -O ${output_bam} --FILTER includeReadList -RLF ${tmp_read_ref_names} -SO coordinate && \
        rm ${tmp_read_ref_names}
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "extracted bam from target region is "${output_bam}
}


function select_fastq {
    # The input path of files should be abs path
    local OPTIND g c b i o a x f r t
    while getopts g::c::b::i::o::a::x::f::r::t:: args
    do 
        case ${args} in
            g) local genelist=$OPTARG ;;
            c) local coordinates=$OPTARG ;;
            b) local region_bed=$OPTARG ;;
            i) local input_bam=$OPTARG ;;
            o) local output_fastq=$OPTARG ;;
            a) local assembly=$OPTARG ;;
            f) local input_forward_read=$OPTARG ;;
            r) local input_reverse_read=$OPTARG ;;
            x) local xa_tag=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bam sample path." ;;
        esac
    done

    module load seqtk/1.3
    
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "We are going to directly output reads from fastq files"

    local bam_dir=$(echo ${input_bam} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}')
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The directory where input bam locating is "${bam_dir}
    local bam_ID=$(get_RG_SM ${input_bam})
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Input bam is ${input_bam} while the sample ID is ${bam_ID}."

    if [[ ${#bam_dir} -gt 1 ]]; then
        local trim_fastq_r1=$(echo ${bam_dir} | awk -F '/' '{for(i=1;i<NF-1;i++) printf $i"/";}' | awk '{printf $0"trimmed_sequences/'${bam_ID}'_1_val_1.fq.gz"}')
        local trim_fastq_r2=$(echo ${bam_dir} | awk -F '/' '{for(i=1;i<NF-1;i++) printf $i"/";}' | awk '{printf $0"trimmed_sequences/'${bam_ID}'_2_val_2.fq.gz"}')
    else
        local trim_fastq_r1=$(pwd | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{printf $0"trimmed_sequences/'${bam_ID}'_1_val_1.fq.gz"}')
        local trim_fastq_r2=$(pwd | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{printf $0"trimmed_sequences/'${bam_ID}'_2_val_2.fq.gz"}')
    fi

    local output_fastq_r1=${output_fastq/.fastq*/_1.fastq}
    local output_fastq_r2=${output_fastq/.fastq*/_2.fastq}

    if [[ ! -f ${trim_fastq_r1} ]] && [[ ! -f ${trim_fastq_r2} ]] && [[ -f ${input_forward_read} ]] && [[ -f ${input_reverse_read} ]]; then
        # If input forward read is like a file path, then use the path to extract reads
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Trimmed fastq r1 file "${trim_fastq_r1}" does not exist, use the input one."
        ls -lh ${input_forward_read}
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Trimmed fastq r2 file "${trim_fastq_r2}" does not exist, use the input one."
        ls -lh ${input_reverse_read}

        local trim_fastq_r1=${input_forward_read}
        local trim_fastq_r2=${input_reverse_read}
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "To be extracted forward fastq file is: "$(ls -lh ${trim_fastq_r1})
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "To be extracted reverse fastq file is: "$(ls -lh ${trim_fastq_r2})
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Output forward reads should be ${output_fastq_r1}"
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Output reverse reads should be ${output_fastq_r2}"

    if [[ -z ${genelist} ]]; then
        if [[ -z ${coordinates} ]]; then
            if [[ -z ${region_bed} ]]; then
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Must specify target regions either by gene names or coordinates. Does not receive any of the two info, exiting now."
                exit 1
            else
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region specified by the region based bed file, the first three columns must be coordinates"
                display_table ${region_bed}
                >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "First sort, then merge bed intervals to a tmp bed file." 
                local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed
                local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).lst
                if [[ $(wc -l ${region_bed} | awk '{printf $1;}') -gt 1 ]]; then
                    bedtools sort -i ${region_bed} | bedtools merge -i - > ${tmp_bed} && \
                    ls -lh ${tmp_bed}
                else
                    cat ${region_bed} > ${tmp_bed}
                    ls -lh ${tmp_bed}
                fi
                >&2 echo $(timestamp)": Temporary bed file ${tmp_bed} looks like:"
                head ${tmp_bed}
                if [[ -z ${xa_tag} ]]; then
                    samtools view -L ${tmp_bed} ${input_bam} | mawk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
                else 
                    samtools view -L ${tmp_bed} ${input_bam} | mawk -F '\t' '$0 ~ /XA:Z:/ || $5 < 30 {print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
                fi
            fi
        else
            >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region sepcified by the coordinates"
            local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).txt
            if [[ -z ${xa_tag} ]]; then
                samtools view ${input_bam} "${coordinates}" | mawk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
            else
                samtools view ${input_bam} "${coordinates}" | mawk -F '\t' '$0 ~ /XA:Z:/{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
            fi
        fi
    else
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "Target region sepcified by the genelist, separated by comma "${genelist}
        local tmp_read_ref_names=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).txt
        local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed
        local tmp_merged_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed
        bash ${central_scripts}/generate_bed_by_gene_names.sh -g ${genelist} -o ${tmp_bed}
        display_table ${tmp_bed}
        bedtools sort -i ${tmp_bed} | bedtools merge -s -i - > ${tmp_merged_bed}
        if [[ -z ${xa_tag} ]]; then
            samtools view -L ${tmp_merged_bed} ${input_bam} | awk -F '\t' '{print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
        else
            samtools view -L ${tmp_merged_bed} ${input_bam} | awk -F '\t' '$0 ~ /XA:Z:/ || $5 < 30 {print $1;}' | sort - | uniq - > ${tmp_read_ref_names}
        fi   
    fi

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Now the query name list overlapping with target region extracted and recorded in ${tmp_read_ref_names}. Let's inspect"
    ls -lh ${tmp_read_ref_names}

    if [[ $(wc -l ${tmp_read_ref_names} | awk '{print $1;}') -eq 0 ]]; then
        echo "WARNING !!!Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): Failed to fetch any reads from ${input_bam}. Try run: "$'\n'"samtools view -L ${tmp_bed} ${input_bam} | mawk -F '\t' '\$0 ~ /XA:Z:/ || \$5 < 30 {print \$1;}' | sort - | uniq - "
        local bed_len=$(awk -F '\t' 'BEGIN{s=0;} {s+=$3-$2;} END{printf "%s",s;}' ${tmp_bed})
        local exon_regions=$(intersectBed -a ${tmp_bed} -b /paedyl01/disk1/yangyxt/public_data/gene_annotation/refgene_latest_anno_exon.bed -wo | awk 'BEGIN{FS=OFS="\t"} {printf "%s:%s-%s_overlap_with_%s(%s,strand_%s);", $1, $2, $3, $7, $8, $9;}' | awk '{gsub(/;$/, ""); print;}')
        >&2 echo "Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): There are $(samtools view -L ${tmp_bed} ${input_bam} | wc -l | awk '{print $1;}') reads align record in ${input_bam} within region of ${tmp_bed}, (None of them are Multi-aligned) the total length of the region is ${bed_len}bp, these regions overlap with exons: ${exon_regions}"
        samtools view -L ${tmp_bed} ${input_bam}
        touch ${output_fastq_r1}
        touch ${output_fastq_r2}
        return
    elif check_fq_file_validity ${trim_fastq_r1} 80 ${threads} && check_fq_file_validity ${trim_fastq_r2} 80 ${threads}; then
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Now since we have some read overlaping with target regions lets check fastq pairs."
        >&2 echo "$(timestamp): In function ${FUNCNAME}: Both ${trim_fastq_r1} and ${trim_fastq_r2} are valid."
        >&2 echo "$(timestamp): In function ${FUNCNAME}: Check whether ${trim_fastq_r1} and ${trim_fastq_r2} have unwanted suffix."

        # Remove suffix if there is one
        local qname_sample=$(seqtk seq ${trim_fastq_r1} | head -1 | awk '{printf $1;}')
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": The query name in fastq file looks like: "$(seqtk seq ${trim_fastq_r1} | head -1 | awk '{printf $1;}')
        if [[ ${qname_sample} =~ \/[1-2]$ ]]; then
            if [[ ${trim_fastq_r1} = \.gz$ ]]; then
                seqtk seq ${trim_fastq_r1} | \
                mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r1::-3} && \
                gzip -f ${trim_fastq_r1::-3}
                # seqtk seq ${trim_fastq_r2} | mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r2::-3} && gzip -f ${trim_fastq_r2::-3}
            else
                seqtk seq ${trim_fastq_r1} | \
                mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r1}.tmp && \
                mv ${trim_fastq_r1}.tmp ${trim_fastq_r1}
                # seqtk seq ${trim_fastq_r2} | mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r2}.tmp && mv ${trim_fastq_r2}.tmp ${trim_fastq_r2}
            fi
        fi
        local qname_sample=$(seqtk seq ${trim_fastq_r2} | head -1 | awk '{printf $1;}')
        if [[ ${qname_sample} =~ \/[1-2]$ ]]; then
            if [[ ${trim_fastq_r2} = \.gz$ ]]; then
                # seqtk seq ${trim_fastq_r1} | mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r1::-3} && gzip -f ${trim_fastq_r1::-3}
                seqtk seq ${trim_fastq_r2} | \
                mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r2::-3} && \
                gzip -f ${trim_fastq_r2::-3}
            else
                # seqtk seq ${trim_fastq_r1} | mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r1}.tmp && mv ${trim_fastq_r1}.tmp ${trim_fastq_r1}
                seqtk seq ${trim_fastq_r2} | \
                mawk '$1 ~ /^@/{gsub(/\/[1-2]$/,"", $1); print;} $1 !~ /^@/{print;}' > ${trim_fastq_r2}.tmp && \
                mv ${trim_fastq_r2}.tmp ${trim_fastq_r2}
            fi
        fi
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}": Now the query name in fastq file looks like: "$(seqtk seq ${trim_fastq_r1} | head -1 )
        rm ${tmp_bed} && \
        >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The query names extracted from ${input_bam} at region ${region_bed} are : "
        echo $(mawk '{printf $1" ";}' < ${tmp_read_ref_names})


        >&2 echo "Now need to perform function extract_read_from_fastq_by_qname on ${trim_fastq_r1} and ${trim_fastq_r2}"
        ls -lh ${tmp_read_ref_names} && \
        /usr/bin/time bash ${central_scripts}/common_bash_utils.sh extract_read_from_fastq_by_qname ${trim_fastq_r1} ${tmp_read_ref_names} ${output_fastq_r1} && \
        /usr/bin/time bash ${central_scripts}/common_bash_utils.sh extract_read_from_fastq_by_qname ${trim_fastq_r2} ${tmp_read_ref_names} ${output_fastq_r2}
        # seqtk subseq ${trim_fastq_r1} ${tmp_read_ref_names} > ${output_fastq_r1} && \
        # seqtk subseq ${trim_fastq_r2} ${tmp_read_ref_names} > ${output_fastq_r2}
        # rm ${tmp_read_ref_names}
    else
        >&2 echo "$(timestamp): In function ${FUNCNAME}: ${trim_fastq_r1} and ${trim_fastq_r2} are not both valid."
        >&2 echo "$(timestamp): In function ${FUNCNAME}: Try to direct extract qname from ${input_bam} using qname list file ${tmp_read_ref_names}"
        /usr/bin/time bash ${central_scripts}/common_bash_utils.sh \
        select_qname_to_paired_fastq_from_bam \
        ${input_bam} \
        ${tmp_read_ref_names} \
        ${output_fastq_r1} \
        ${output_fastq_r2}
    fi  

    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "extracted fastq files from target region are ${output_fastq_r1},${output_fastq_r2}"
    ls -lh ${output_fastq_r1}
    ls -lh ${output_fastq_r2}
    module unload seqtk/1.3

    trap "rm ${tmp_read_ref_names} ${tmp_bed} ${tmp_merged_bed}" 0
}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    while getopts g::c::b::i:o::a::x::f::r::t:: args; do 
        case ${args} in
            g) genelist=$OPTARG ;;
            c) coordinates=$OPTARG ;;
            b) region_bed=$OPTARG ;;
            i) input_bam=$OPTARG ;;
            o) output_file=$OPTARG ;;
            a) assembly=$OPTARG ;;
            x) xa_tag=$OPTARG ;;
            f) forward_read=$OPTARG ;;
            r) reverse_read=$OPTARG ;;
            t) threads=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bam sample path." ;;
        esac
    done
    current_conda=$(return_to_conda_base)
    if [[ -z ${assembly} ]]; then assembly=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta; fi
    if [[ ${output_file} =~ \.bam$ ]]; then select_bam "$@"; fi
    if [[ ${output_file} =~ \.fastq(\.gz)*$ ]]; then select_fastq "$@"; fi
    if [[ ${current_conda} != "base" ]]; then
        source "$(which conda | awk -F '/' '{for (i=1;i<NF-1;i++) printf "%s/", $i;}')etc/profile.d/conda.sh" && \
        conda activate ${current_conda} || \
        /home/yangyxt/anaconda3/bin/conda activate ${current_conda}
    fi
fi



            
    
