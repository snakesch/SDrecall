# This is the script to showcase the arguments we used to perform GATK Best Practice calling on GIAB samples

function BQSR {
    local input_bam=${1}
    local ref_genome=${2}

    if [[ -z ${ref_genome} ]]; then
        local ref_genome="<path>/indexed_genome/ucsc.hg19.fasta"
    fi

    local ref_gen=$(dirname ${ref_genome})
    local gatk="bash <path>/ngs_scripts/common_bash_utils.sh gatk_wrapper"

    if [[ ${ref_genome} =~ hg19 ]]; then
        time ${gatk} BaseRecalibrator \
        -R ${ref_genome} \
        -I ${input_bam} \
        --known-sites $ref_gen/1000G_phase1.indels.hg19.sites.vcf \
        --known-sites $ref_gen/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
        --known-sites $ref_gen/dbsnp_138.hg19.vcf \
        -O ${input_bam/.bam/}.bqsr.table \
        -LE TRUE && \
        >&2 echo "***** Generating BQSR table for hg19 done ******"
    else
        time ${gatk} BaseRecalibrator \
        -R ${ref_genome} \
        -I ${input_bam} \
        --known-sites $ref_gen/Homo_sapiens_assembly38.known_indels.vcf.gz \
        --known-sites $ref_gen/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --known-sites $ref_gen/Homo_sapiens_assembly38.dbsnp138.vcf \
        -O ${input_bam/.bam/}.bqsr.table \
        -LE TRUE && \
        >&2 echo "***** Generating BQSR table done ******"
    fi

    # Gold_standard small indels 
    # Gold_standard indels from another project named 1000G
    # Gold standard SNV
    # Designates the output directory

    # The second step is to use ApplyBQSR to change the current quality scores

    time ${gatk} ApplyBQSR \
    -R ${ref_genome} \
    -I ${input_bam} \
    -bqsr-recal-file ${input_bam/.bam/}.bqsr.table \
    -O ${input_bam/.bam/.bqsr.bam} && \
    >&2 echo "***BQSR DONE for ${ref_genome}***"
}


function run_BQSR_and_HC_and_genotype_unit(){
    local input_bam=${1}
    local ref_genome=${2}

    if [[ -z ${ref_genome} ]]; then
        local ref_genome="<path>/indexed_genome/ucsc.hg19.fasta"
    fi

    local tag="gatk"
    local input_bam_ID=$(basename ${input_bam/.SD.deduped.bam/})
    local bamf=$(dirname ${input_bam})
    local gatk="bash ${central_scripts}/common_bash_utils.sh gatk_wrapper"
    local working_bam=${input_bam/.bam/.bqsr.bam}
    local output_vcf=${input_bam/.SD.deduped.bam/.${tag}.vcf.gz}  # Now the output name should be HG002.gatk.vcf.gz
    
    if [[ ! -f ${input_bam}.bai ]] || [[ ${input_bam} -nt ${input_bam}.bai ]]; then
        samtools index ${input_bam}
    fi

    if [[ ${working_bam} -nt ${input_bam} ]] && check_bam_validity ${working_bam}; then
        >&2 echo $'\n\n'"Line "${LINENO}": In function ${FUNCNAME}: $(timestamp): SKIP BQSR bam"
    else
        BQSR \
        ${input_bam} \
        ${ref_genome} && \
        samtools index ${working_bam}
    fi

    if check_vcf_validity ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz && [[ ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz -nt ${working_bam} ]]; then
        >&2 echo "$(timestamp): In function ${FUNCNAME}, ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz already generated. Check it out."
        ls -lh ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz
    else
        time ${gatk} --java-options "-Xmx40G -XX:ConcGCThreads=20 -XX:ParallelGCThreads=20" \
            HaplotypeCaller \
            --emit-ref-confidence GVCF \
            -R ${ref_genome} \
            -I ${working_bam} \
            -O ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz \
            --assembly-region-out ${bamf}/${input_bam_ID}.HC.active_region.tsv \
            --bam-writer-type CALLED_HAPLOTYPES \
            --native-pair-hmm-threads 60 \
            -bamout ${bamf}/${input_bam_ID}.SD.HC.realigned.bam \
            --linked-de-bruijn-graph \
            -G StandardAnnotation \
            -G AS_StandardAnnotation \
            -G StandardHCAnnotation \
            -OBI "true" && \
        >&2 echo "** GVCF ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz done **"
    fi

    if check_vcf_validity ${output_vcf} && [[ ${output_vcf} -nt ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz ]]; then
         >&2 echo "$(timestamp): In function ${FUNCNAME}, ${output_vcf} already generated and validated"
    else
        /usr/bin/time ${gatk} --java-options "-Xmx50g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
            -R ${ref_genome} \
            --tmp-dir "<path>/test_tmp" \
            -V ${bamf}/${input_bam_ID}.SD.HC.g.vcf.gz \
            -G StandardAnnotation \
            -G AS_StandardAnnotation \
            -O ${output_vcf} && \
        bcftools sort -Oz -o ${output_vcf} ${output_vcf} && \
        tabix -f -p vcf ${output_vcf}
    fi
}






############################################################################################################################################################
# Below we paste the commands to showcase the arguments we used for DeepVariant calling

function run_deepvariant() {
    local bamID=${1}
    local reference_genome=${2}
    local output_dir=${3}

    local genome_tag=$(basename ${reference_genome} | cut -f 2 -d ".")
    local sd_working_dir="<path>/wgs/GIAB_samples/aligned_results/${genome_tag}/${bamID}_refSD_priority_component_pairs"
    local working_bam="$(dirname ${sd_working_dir})/${bamID}.SD.deduped.bam"
    local output_vcf="${output_dir}/${bamID}.deepvariant.vcf.gz"
    local output_gvcf="${output_dir}/${bamID}.deepvariant.g.vcf.gz"
    # Remember that this BAM is now mapped to standard ucsc.hg19.fasta

    bash ${central_scripts}/common_bash_utils.sh \
    ind_deepvariant \
    -r ${reference_genome} \
    -b ${working_bam} \
    -o ${output_vcf} \
    -g ${output_gvcf} \
    -q "wgs" \
    -u 4 && \
    display_vcf ${output_vcf}
}

function ind_deepvariant {
    local OPTIND r b o p g t s c q u
    while getopts r::b:o:p::g::t::s::c::q::u:: args
    do
        case ${args} in
            r) local ref_genome=$OPTARG ;;
            b) local bam_file=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            g) local output_gvcf=$OPTARG ;;
            t) local target_region=$OPTARG;;
            p) local ped_file=$OPTARG ;;
            s) local samples=$OPTARG ;;
            c) local dv_container=$OPTARG ;;
            q) local seq_type=$OPTARG ;;
            u) local total_cpu=$OPTARG ;;
            *) echo "No argument passed. Pls at least specify -r (ref fasta path) or -b (bed_file path)." ;;
        esac 
    done
    module load singularity/3.7.2

    if [[ ! -z ${target_region} ]]; then
        local region_arg="--regions=${target_region}"
    else
        local region_arg=""
    fi

    if [[ ${seq_type} == "wgs" ]]; then
        local model_type="WGS"
    else
        local model_type="WES"
    fi

    if [[ -z ${dv_container} ]]; then
        local dv_container="<path>/Tools/DeepVariant/deepvariant_1.5.0.sif"
    fi

    if [[ -z ${ref_genome} ]]; then
        local ref_genome=<path>/indexed_genome/ucsc.hg19.fasta
    fi

    # model type options: exactly one of the following [WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]
    if [[ -z ${total_cpu} ]]; then
        local total_cpu=$(get_pbs_cpu)
        if [[ -z ${total_cpu} ]]; then
            local total_cpu=1
        fi
    fi

    local threads=$(determine_job_num -m 1 -c 1 -t ${total_cpu})
    local tmp_tag=$(randomID)
    local tmp_files=<path>/test_tmp/${tmp_tag}.txt
    local singularity_run=singularity_run_${tmp_tag}
    local singularity_tmp=singularity_tmp_${tmp_tag}
    local singularity_cache=singularity_cache_${tmp_tag}
    local singularity_wkd=singularity_wkd_${tmp_tag}
    local singularity_inter=singularity_inter_${tmp_tag}

    mkdir -p \
    <path>/test_tmp/${singularity_run} \
    <path>/test_tmp/${singularity_tmp} \
    <path>/test_tmp/${singularity_cache} \
    <path>/test_tmp/${singularity_wkd} \
    <path>/test_tmp/${singularity_inter} 2> /dev/null || :

    find $HOME -maxdepth 1 -type f -name ".*" | awk '{printf "basename %s\n", $0;}' | bash - > ${tmp_files} && \
    rsync -avu --delay-updates --dry-run --files-from=${tmp_files} $HOME/ <path>/home/ && \
    rsync -avu --files-from=${tmp_files} $HOME/ <path>/home/ && \
    rm -f ${tmp_files} || ls -lh ${tmp_files}

    trap 'silent_remove_tmps \
    <path>/test_tmp/${singularity_run} \
    <path>/test_tmp/${singularity_tmp} \
    <path>/test_tmp/${singularity_cache} \
    <path>/test_tmp/${singularity_wkd} \
    <path>/test_tmp/${singularity_inter} \
    ${dv_container/.sif/}.${tmp_tag}.sif' SIGTERM

    log "The DeepVariant running command is: \
    singularity exec \
    -B "<path>,/usr/lib/locale,<path>/test_tmp/${singularity_run}:/tmp" \
    --env LANG="en_US.UTF-8" \
    --env LC_ALL="C" \
    --env LANGUAGE="en_US.UTF-8" \
    --env LC_CTYPE="UTF-8" \
    --env OMP_NUM_THREADS="${threads}" \
    --env TF_NUM_INTRAOP_THREADS="${threads}" \
    --env TF_NUM_INTEROP_THREADS="${threads}" \
    --env TMPDIR="<path>/test_tmp/${singularity_run}" \
    --env SINGULARITY_TMPIDR="<path>/test_tmp/${singularity_tmp}" \
    --env SINGULARITY_CACHEDIR="<path>/test_tmp/${singularity_cache}" \
    --home "<path>/home:/home" \
    --workdir "<path>/test_tmp/${singularity_wkd}" \
    --contain \
    ${dv_container/.sif/}.${tmp_tag}.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model_type} \
    --ref="${ref_genome}" \
    --reads="${bam_file}" \
    --make_examples_extra_args="normalize_reads=true" \
    ${region_arg} \
    --output_vcf="${output_vcf}" \
    --output_gvcf="${output_gvcf}" \
    --intermediate_results_dir="<path>/test_tmp/${singularity_inter}" \
    --num_shards=${threads}"

    cp -f ${dv_container} ${dv_container/.sif/}.${tmp_tag}.sif && \
    ls -lh ${dv_container/.sif/}.${tmp_tag}.sif && \
    singularity exec \
    --no-home \
    --cleanenv \
    -B "<path>,/usr/lib/locale,<path>/test_tmp/${singularity_run}:/tmp" \
    --env LANG="en_US.UTF-8" \
    --env LC_ALL="C" \
    --env LANGUAGE="en_US.UTF-8" \
    --env LC_CTYPE="UTF-8" \
    --env OMP_NUM_THREADS="${threads}" \
    --env TF_NUM_INTRAOP_THREADS="${threads}" \
    --env TF_NUM_INTEROP_THREADS="${threads}" \
    --env TMPDIR="<path>/test_tmp/${singularity_run}" \
    --env SINGULARITY_TMPIDR="<path>/test_tmp/${singularity_tmp}" \
    --env SINGULARITY_CACHEDIR="<path>/test_tmp/${singularity_cache}" \
    --workdir "<path>/test_tmp/${singularity_wkd}" \
    --containall \
    ${dv_container/.sif/}.${tmp_tag}.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model_type} \
    --ref="${ref_genome}" \
    --reads="${bam_file}" \
    --make_examples_extra_args="normalize_reads=true" \
    ${region_arg} \
    --output_vcf="${output_vcf}" \
    --output_gvcf="${output_gvcf}" \
    --intermediate_results_dir="<path>/test_tmp/${singularity_inter}" \
    --num_shards=${threads} && \
    ls -lh ${output_vcf} && \
    ls -lh ${output_gvcf} && \
    silent_remove_tmps \
    <path>/test_tmp/${singularity_run} \
    <path>/test_tmp/${singularity_tmp} \
    <path>/test_tmp/${singularity_cache} \
    <path>/test_tmp/${singularity_wkd} \
    ${dv_container/.sif/}.${tmp_tag}.sif \
    <path>/test_tmp/${singularity_inter} || { \
    silent_remove_tmps \
    <path>/test_tmp/${singularity_run} \
    <path>/test_tmp/${singularity_tmp} \
    <path>/test_tmp/${singularity_cache} \
    <path>/test_tmp/${singularity_wkd} \
    <path>/test_tmp/${singularity_inter} \
    ${dv_container/.sif/}.${tmp_tag}.sif; \
    return 1; }

    module unload singularity
}