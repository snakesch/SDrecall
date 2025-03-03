self_script=$(realpath $0)
project_dir=$(dirname $(dirname ${self_script}))
central_scripts=${project_dir}/shell_utils.sh

source ${central_scripts}

# This script is a showcase of how we performed the benchmarking for the SDrecall results.
# Due to extended requirements of environments and dependency files, we cannot offer a fully portable benchmark workflow to reproduce the entire benchmarking process
# This script is the main entry point of the benchmarking workflow, while the dependent scripts and functions are either stored in the same directory or defined in the same script. 

function annotate_cadd() {
    local called_vcf=${1}
    local output_vcf=${2}
    local genome_tag=${3}

    if [[ -z ${genome_tag} ]]; then
        local genome_tag="GRCh37"
    fi

    source /home/yangyxt/miniforge3/etc/profile.d/conda.sh && source /home/yangyxt/miniforge3/etc/profile.d/mamba.sh
    mamba activate ngs_pipeline

    bash ${central_scripts}/annotation_per_family.sh \
    Calculate_CADD \
    -i ${called_vcf} \
    -o ${called_vcf/.vcf*/.cadd.txt} \
    -t 10 \
    -g ${genome_tag} && \
    mawk 'BEGIN{OFS=FS="\t";} NR > 1 {printf "chr%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $NF, $21;}' ${called_vcf/.vcf*/.cadd.txt} > ${called_vcf/.vcf*/.cadd.tsv} && \
    bgzip -f -c ${called_vcf/.vcf*/.cadd.tsv} > ${called_vcf/.vcf*/.cadd.tsv}.gz && \
    tabix -f -s 1 -b 2 -e 2 ${called_vcf/.vcf*/.cadd.tsv}.gz && \
    echo "##INFO=<ID=PHRED,Number=1,Type=Float,Description=\"PHRED scaled score of deleterious effect prediction of the variant\">" > ${called_vcf/.vcf*/.cadd.header} && \
    echo "##INFO=<ID=Gene.refGene,Number=1,Type=String,Description=\"Symbol of the Gene the variant overlaps with\">" >> ${called_vcf/.vcf*/.cadd.header} && \
    bcftools annotate -a ${called_vcf/.vcf*/.cadd.tsv}.gz -h ${called_vcf/.vcf*/.cadd.header} -c CHROM,POS,REF,ALT,PHRED,Gene.refGene -Ov ${called_vcf} | \
    bcftools filter -e 'INFO/PHRED >= 20' -s 'CADD_deleterious' -m + | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf}

    mamba deactivate
}


function annotate_gnomAD_common () {
    local input_vcf=${1}
    local output_vcf=${2}
    local assembly_version=${3}

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="hg19"
    fi

    if [[ ${assembly_version} == "hg19" ]] || [[ ${assembly_version} == "GRCh37" ]]; then
        local genome_gnomad_header="gnomAD_genome_ALL"
        local exome_gnomad_header="gnomAD_exome_ALL"
    elif [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        local genome_gnomad_header="gnomad41_genome_AF"
        local exome_gnomad_header="gnomad41_exome_AF"
    fi

    local annovar_vcf
    basic_annovar ${input_vcf} annovar_vcf ${assembly_version} && \
    display_table ${annovar_vcf} && \
    bgzip -c ${annovar_vcf} > ${annovar_vcf}.gz && \
    tabix -f -p vcf ${annovar_vcf}.gz && \
    bcftools annotate -a ${annovar_vcf}.gz -c CHROM,POS,REF,ALT,INFO/${genome_gnomad_header},INFO/${exome_gnomad_header},INFO/ExonicFunc.refGene,INFO/Gene.refGene -Ov ${input_vcf} | \
    bcftools filter -e "INFO/${genome_gnomad_header} >= 0.05 || INFO/${exome_gnomad_header} >= 0.05" -s 'gnomAD_BA' -m + -Ov - | \
    bcftools filter -e "INFO/${genome_gnomad_header} >= 0.01 || INFO/${exome_gnomad_header} >= 0.01" -s 'gnomAD_common' -m + -Ov - | \
    bcftools sort -Oz -o ${output_vcf} - && \
    tabix -f -p vcf ${output_vcf} && \
    display_table ${output_vcf}
}


function annotate_golden_record () {
    local input_vcf=${1}
    local gold_vcf=${2}
    local output_vcf=${3}

    source /home/yangyxt/miniforge3/etc/profile.d/conda.sh && source /home/yangyxt/miniforge3/etc/profile.d/mamba.sh
    mamba activate SDrecall

    log "Annotate golden records on ${input_vcf} from ${gold_vcf} and output to ${output_vcf}."

    bcftools view -h ${gold_vcf} > ${gold_vcf/.vcf*/.header}
    if [[ $(awk '$0 ~ /platforms/{print}' ${gold_vcf/.vcf*/.header} | wc -l) -eq 0 ]]; then
        awk '$0 ~ /^#CHROM/{printf "##INFO=<ID=platforms,Number=1,Type=Integer,Description=\"Number of platforms the variant is detected\">\n";} \
             $0 !~ /^#CHROM/{print;}' ${gold_vcf/.vcf*/.header} > ${gold_vcf/.vcf*/.header.tmp} && \
        mv ${gold_vcf/.vcf*/.header.tmp} ${gold_vcf/.vcf*/.header}
    else
        awk '$0 !~ /^#CHROM/{print;}' ${gold_vcf/.vcf*/.header} > ${gold_vcf/.vcf*/.header.tmp} && \
        mv ${gold_vcf/.vcf*/.header.tmp} ${gold_vcf/.vcf*/.header}
    fi

    bcftools query -i 'FILTER == "."' -f '%CHROM\t%POS\t%REF\t%ALT\n' ${gold_vcf} | \
    mawk 'BEGIN{OFS=FS="\t";} {printf "%s\t1\n", $0;}' > ${gold_vcf/.vcf*/.MHC.var.tsv} && \
    bgzip -f ${gold_vcf/.vcf*/.MHC.var.tsv} && \
    tabix -f -s 1 -b 2 -e 2 ${gold_vcf/.vcf*/.MHC.var.tsv}.gz && \
    log "The header is ${gold_vcf/.vcf*/.header}" && \
    bcftools annotate -a ${gold_vcf/.vcf*/.MHC.var.tsv}.gz -h ${gold_vcf/.vcf*/.header} -c CHROM,POS,REF,ALT,INFO/platforms -Oz -o ${gold_vcf/.vcf/.anno.vcf} ${gold_vcf} && \
    bcftools index -t ${gold_vcf/.vcf/.anno.vcf} && \

    bcftools view -i 'GT=="hom"' -Oz -o ${gold_vcf/.vcf/.anno.hom.vcf} ${gold_vcf/.vcf/.anno.vcf} && \
    tabix -f -p vcf ${gold_vcf/.vcf/.anno.hom.vcf} && \
    bcftools view -i 'GT=="het"' -Oz -o ${gold_vcf/.vcf/.anno.het.vcf} ${gold_vcf/.vcf/.anno.vcf} && \
    tabix -f -p vcf ${gold_vcf/.vcf/.anno.het.vcf} && \
    bcftools annotate -a ${gold_vcf/.vcf/.anno.hom.vcf} -c CHROM,POS,REF,ALT,INFO/platforms -Ov ${input_vcf} | \
    bcftools filter -e 'INFO/platforms >= 1' -s 'golden_variant_hom' -m + | \
    bcftools sort -Oz -o ${output_vcf/.vcf/.hom.vcf} && \
    tabix -f -p vcf ${output_vcf/.vcf/.hom.vcf} && \
    display_vcf ${output_vcf/.vcf/.hom.vcf} && \
    bcftools annotate -a ${gold_vcf/.vcf/.anno.het.vcf} -c CHROM,POS,REF,ALT,INFO/platforms -Ov ${output_vcf/.vcf/.hom.vcf} | \
    bcftools filter -e 'INFO/platforms >= 1 && FILTER !~ "golden_variant_hom"' -s 'golden_variant_het' -m + | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    display_table ${output_vcf}
}



function bench_callset_per_sample() {
    # Default values
    local sampID=""
    local caller_tag=""
    local golden_vcf=""
    local bench_region=""
    local called_vcf=""
    local called_bam=""
    local target_region=""
    local recall_region=""
    local original_bam=""
    local result_meta_file=""

    # Parse command line arguments
    while getopts "s:c:g:b:v:a:t:r:o:m:" opt; do
        case ${opt} in
            s) sampID=${OPTARG} ;;
            c) caller_tag=${OPTARG} ;;
            g) golden_vcf=${OPTARG} ;;
            b) bench_region=${OPTARG} ;;
            v) called_vcf=${OPTARG} ;;
            a) called_bam=${OPTARG} ;;
            t) target_region=${OPTARG} ;;
            r) recall_region=${OPTARG} ;;
            o) original_bam=${OPTARG} ;;
            m) result_meta_file=${OPTARG} ;;
            \?) echo "Invalid option: -${OPTARG}" >&2; return 1 ;;
        esac
    done
    shift $((OPTIND -1))

    # Check required parameters
    if [[ -z "${sampID}" ]]; then
        echo "Sample ID (-s) is required"
        return 1
    fi
    
    if [[ -z "${caller_tag}" ]]; then
        echo "Caller tag (-c) is required"
        return 1
    fi

    # Set default paths if not provided
    if [[ -z "${golden_vcf}" ]]; then
        golden_vcf="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/${sampID}_GRCh38_1_22_v4.2.1_benchmark.norm.vcf.gz"
    fi
    
    if [[ -z "${bench_region}" ]]; then
        bench_region="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/${sampID}_GRCh38_1_22_v4.2.1_benchmark.bed"
    fi
    
    if [[ -z "${called_vcf}" ]]; then
        called_vcf="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/${sampID}.sdrecall.ref.${caller_tag}.merged.vcf.gz"
    fi
    
    if [[ -z "${called_bam}" ]]; then
        called_bam="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/${sampID}.SD.deduped.bam"
    fi
    
    if [[ -z "${target_region}" ]]; then
        target_region="/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.40_GRCh38.p14_genomic.func.coding.sorted.pad20.bed"
    fi
    
    if [[ -z "${recall_region}" ]]; then
        recall_region="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/${sampID}_refSD_priority_component_pairs/all_PC_regions.bed"
    fi
    
    if [[ -z "${result_meta_file}" ]]; then
        result_meta_file="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/hg38/${sampID}.bench.${caller_tag}.meta.tsv"
    fi
    
    if [[ -z "${original_bam}" ]]; then
        original_bam="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/${sampID}.GRCh38.300x.bam"
        if [[ ${sampID} == "HG006" ]] || [[ ${sampID} == "HG007" ]]; then
            original_bam=${original_bam/.300x/.100x}
        fi
    fi


    if [[ ! -f ${output_poor_bed/.bed/.bench.bed} ]] && \
       [[ ${called_bam} -nt ${output_poor_bed/.bed/.bench.bed} ]] && \
       [[ ${bench_region} -nt ${output_poor_bed/.bed/.bench.bed} ]]; then
        local targeted_golden=${golden_vcf/.vcf*/.bench.vcf.gz}
        local targeted_called=${called_vcf/.vcf*/.bench.vcf.gz}
        local output_poor_bed=${called_bam/.bam/.multi_aligned.bed}
        python3 ${central_scripts}/pick_multialign_region.py \
        -f main_func_pick_region \
        -a ${called_bam} \
        -k "output_bed=${output_poor_bed};inferred_coverage=False;target_region=${target_region};target_tag=FCRs;MQ_threshold=41;depth_threshold=10;minimum_depth=3;multialign_frac=0.5;threads=4" && \
        ls -lhtr ${output_poor_bed} && \
        bedtools intersect -a ${output_poor_bed} -b ${recall_region} | \
        bedtools intersect -a stdin -b ${bench_region} > ${output_poor_bed/.bed/.bench.bed} || \
        { log "Failed to generate benchmark region. Exit."; return 1; }
    else
        log "The benchmark region ${output_poor_bed/.bed/.bench.bed} is already generated."
    fi

    log "Though specify the target region ${target_region}, we use it to overlap with poor cov region (Depth < 10 when MQ > 10) in ${called_bam}. Stored in ${output_poor_bed}" && \
    display_table ${output_poor_bed/.bed/.bench.bed}
    if [[ $(cat ${output_poor_bed/.bed/.bench.bed} | wc -l) -gt 0 ]] && \
       [[ ${output_poor_bed/.bed/.bench.bed} -nt ${targeted_golden} ]] && \
       [[ ${output_poor_bed/.bed/.bench.bed} -nt ${targeted_called} ]]; then
        bcftools view -R ${output_poor_bed/.bed/.bench.bed} -Oz -o ${targeted_golden} ${golden_vcf} && \
        tabix -f -p vcf ${targeted_golden} && \
        bcftools view -R ${output_poor_bed/.bed/.bench.bed} -Oz -o ${targeted_called} ${called_vcf} && \
        tabix -f -p vcf ${targeted_called} && \
        display_table ${targeted_golden} && \
        display_table ${targeted_called}
    else
        log "The poor covered region in ${called_bam} overlapping with ${target_region} is empty."
    fi

    local original_bam=/paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/${sampID}.GRCh38.300x.bam
    if [[ ${sampID} == "HG006" ]] || [[ ${sampID} == "HG007" ]]; then
        local original_bam=${original_bam/.300x/.100x}
    fi

    python3 /paedyl01/disk1/yangyxt/ngs_scripts/identify_allele_bias_missing_TP.py \
    -gv ${targeted_golden} \
    -tv ${targeted_called} \
    --bam_300x ${original_bam} \
    --bam_30x /paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/${sampID}.SD.deduped.bam \
    -rg /paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta && \
    ls -lht ${targeted_golden/.vcf.gz/.no_lost_alt.vcf.gz} && \
    local targeted_golden=${targeted_golden/.vcf.gz/.no_lost_alt.vcf.gz}

    log "Now, after removing the TPs with no alt alleles in input BAM files, the to be compared golden vcf is ${targeted_golden} and the called vcf is ${targeted_called}"

    # Tag out the variants with gnomAD common and CADD_deleterious (gnomAD_common and CADD deleterious filter tags are determined based on ANNOVAR annotation results)
    annotate_gnomAD_common ${targeted_called} ${targeted_called/.vcf/.gnomad.vcf} "hg38" && \
    annotate_cadd ${targeted_called/.vcf/.gnomad.vcf} ${targeted_called/.vcf/.cadd.vcf} "GRCh38" && \
    annotate_golden_record ${targeted_called/.vcf/.cadd.vcf} ${targeted_golden} ${targeted_called/.vcf/.final.vcf} || \
    { log "Failed to add annotations"; return 1; }

    bcftools query \
    -i '(FILTER~"golden_variant_het" || FILTER~"golden_variant_hom")' \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%INFO/Gene.refGene\n' \
    ${targeted_called/.vcf/.final.vcf} > ${targeted_called/.vcf*/.final.tsv}

    : > ${result_meta_file}
    # Perform calculation sample-wise and record the result to a sample-wise meta file.
    if check_vcf_validity ${targeted_called} && check_vcf_validity ${targeted_golden}; then
        # awk -v gv="${targeted_golden}" -v cv="${targeted_called}" 'BEGIN{printf "golden_vcf\tcalled_vcf\n%s\t%s", gv, cv;}' > ${result_meta_file}
        log " The meta file to record the precision and recall rate is ${result_meta_file}."
        python3 /paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/calculate_precision_recall_vcf.py \
        -mf "${result_meta_file}" \
        -ot ${caller_tag} \
        -d 10 \
        -b ${output_poor_bed} \
        -gc "${targeted_golden}" \
        -cc "${targeted_called/.vcf/.final.vcf}" && \
        display_table ${result_meta_file}
    else
        log " Either ${targeted_called} is not valid or ${targeted_golden} is not valid. Drop calculating the precision and recall."
    fi
}


function batch_perform_benchmarking () {
	local -a sample_IDs=( "HG002" "HG003" "HG004" "HG005" "HG006" "HG007" )
	local -a caller_tags=( "GATK" "DeepVariant" )
	local -a assemblies=( "hg38" "hg19" )
	local bench_region_v421="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/{2}/golden_vcfs/{1}_{2}_v4.2.1_benchmark.bed"
	local target_region="/paedyl01/disk1/yangyxt/public_data/gene_annotation/{2}_genomic.func.coding.sorted.pad20.bed"
	local result_meta_file="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/{1}/{2}.bench.{3}.meta.tsv"
	local called_vcf="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/{2}/{1}.sdrecall.ref.{3}.merged.vcf.gz"
	local called_bam="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{2}/{1}.SD.deduped.bam"
	local original_bam="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/{1}.{2}.300x.bam"
	local golden_vcf_v421="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/{2}/golden_vcfs/{1}_{2}_v4.2.1_benchmark.norm.vcf.gz"
	
	
	parallel -j6 \
	bench_callset_per_sample \
	-s {1} \
	-c {3} \
	-g ${golden_vcf_v421} \
	-b ${bench_region_v421} \
	-v ${called_vcf} \
	-a ${called_bam} \
	-t ${target_region} \
	-r ${recall_region} \
	-o ${original_bam} \
	-m ${result_meta_file} ::: ${sample_IDs[@]} ::: ${assemblies[@]} ::: ${caller_tags[@]}


	local golden_vcf_cmrg="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/{1}/golden_vcfs/HG002_{1}_CMRG_smallvar_v1.00.norm.small.vcf.gz"
	local bench_region_cmrg="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/{1}/golden_vcfs/HG002_{1}_CMRG_smallvar_v1.00.norm.small.bed"

	parallel -j4 \
	bench_callset_per_sample \
	-s {1} \
	-c {3} \
	-g ${golden_vcf_cmrg} \
	-b ${bench_region_cmrg} \
	-v ${called_vcf} \
	-a ${called_bam} \
	-t ${target_region} \
	-r ${recall_region} \
	-o ${original_bam} \
	-m ${result_meta_file} ::: "HG002" ::: ${assemblies[@]} ::: ${caller_tags[@]}

}
