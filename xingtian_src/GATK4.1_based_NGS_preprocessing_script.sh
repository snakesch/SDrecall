#!/bin/bash
# This script is based on GATK4.1 which is released in Jan 2019. 
# Define a var for iteration of processing which is $1-- the first section of text cut by awk -F
# Define vars for addresses of required files
. /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh
. /paedyl01/disk1/yangyxt/ngs_scripts/realign_masked_and_HC_multiploidy.sh

seq_type=wes
sample_batch=new_14_samples_2019CNY
wkd=/paedyl01/disk1/yangyxt/${seq_type}/${sample_batch}
rawf=${wkd}/raw_data
trimf=${wkd}/trimmed_sequences
fastqf=${wkd}/fastqcresults
bamf=${wkd}/aligned_results
gvcf=${wkd}/gvcfs
ref_gen=/paedyl01/disk1/yangyxt/indexed_genome # Where I have indexed g1k_v37 and ucsc_hg19 ref hum genome, as well as the standard vcf files.
gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
ped=${wkd}/pedigree
scrf=${wkd}/respective_scripts
batchf=${wkd}/job_qsub_scripts
vcf=${wkd}/vcfs
central_scripts=/paedyl01/disk1/yangyxt/ngs_scripts
# save_py=/paedyl01/disk1/yangyxt/ngs_scripts/pysam_save_mq.py

sample_ID=$1
if [[ -z $2 ]]; then special_interval="no_low_complex_interval"; else special_interval=$2; fi

du -h --max-depth=0 /staging 2> /dev/null && tmp_dir="/staging/test_tmp" || tmp_dir="/paedyl01/disk1/yangyxt/test_tmp"

#software loading
module load bwa
module load samtools
module load java
PICARD=/home/yangyxt/software/picard/picard.jar

function migrate_to_staging_per_sample(){
	trap 'find /staging/${seq_type}/${sample_batch}/${rawf} -name \"${sample_ID}*q.gz\" -type f -exec rm -f {} \; return 0' ERR
	local usage=$(du -h --max-depth=0 /staging | awk -F '\t' '{prinf $1;}')
	du -h --max-depth=0 /staging 2> /dev/null || local ssd_avail="not_available"
	if [[ ${ssd_avail} == "not_available" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This job is not running on paedyl02 node, abandon migrate the process to staging."
		return 0 
	elif [[ ${usage} =~ ^([4-5]\.[5-9]+T|6\.[0-9]+T) ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": /staging folder is nearly full. Quit copy the job to it."
		return 0 
	elif [[ ${seq_type} == "wgs" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": /staging path is too small for WGS data."
		return 0
	elif [[ ${#sample_ID} -gt 3 ]]; then
		mkdir -p /staging/${seq_type}/${sample_batch}/CNV_calling 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/CNV_calling"
		mkdir -p /staging/${seq_type}/${sample_batch}/aligned_results 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/aligned_results"
		mkdir -p /staging/${seq_type}/${sample_batch}/annotations 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/annotations"
		mkdir -p /staging/${seq_type}/${sample_batch}/fastqcresults 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/fastqcresults"
		mkdir -p /staging/${seq_type}/${sample_batch}/gvcfs 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/gvcfs"
		mkdir -p /staging/${seq_type}/${sample_batch}/job_qsub_scripts 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/job_qsub_scripts"
		mkdir -p /staging/${seq_type}/${sample_batch}/respective_scripts 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/respective_scripts"
		mkdir -p /staging/${seq_type}/${sample_batch}/trimmed_sequences 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/trimmed_sequences"
		mkdir -p /staging/${seq_type}/${sample_batch}/vcfs 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/vcfs"
		mkdir -p /staging/${seq_type}/${sample_batch}/pedigree 2> /dev/null || echo "Already created /staging/${seq_type}/${sample_batch}/pedigree"
		
		rsync --dry-run --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}_*.f*q.gz" --exclude="*" ${rawf}/ /staging/${seq_type}/${sample_batch}/raw_data/
		rsync --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}_*.f*q.gz" --exclude="*" ${rawf}/ /staging/${seq_type}/${sample_batch}/raw_data/

		ori_wkd_parent_dir=/paedyl01/disk1/yangyxt/${seq_type}
		wkd=/staging/${seq_type}/${sample_batch}
		rawf=${wkd}/raw_data
		trimf=${wkd}/trimmed_sequences
		fastqf=${wkd}/fastqcresults
		bamf=${wkd}/aligned_results
		gvcf=${wkd}/gvcfs
		ped=${wkd}/pedigree
		scrf=${wkd}/respective_scripts
		batchf=${wkd}/job_qsub_scripts
		vcf=${wkd}/vcfs
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": sample ID not correctly specified: sample_ID=${sample_ID}"
	fi
}


function migrate_from_staging_per_sample(){
	local usage=$(du -h --max-depth=0 /staging | awk -F '\t' '{prinf $1;}')
	du -h --max-depth=0 /staging 2> /dev/null || local ssd_avail="not_available"
	if [[ ${ssd_avail} == "not_available" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This job is not running on paedyl02 node, abandon migrate the process to staging."
		return 0
	elif [[ ${seq_type} == "wgs" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": /staging path is too small for WGS data."
		return 0
	elif [[ ${#sample_ID} -gt 3 ]]; then
		# First migrate back the trimmed sequences and FastQC report
		rsync --dry-run --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}_*" --exclude="*" /staging/${seq_type}/${sample_batch}/trimmed_sequences/ ${ori_wkd_parent_dir}/${sample_batch}/trimmed_sequences
		rsync --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}_*" --exclude="*" /staging/${seq_type}/${sample_batch}/trimmed_sequences/ ${ori_wkd_parent_dir}/${sample_batch}/trimmed_sequences
		# Then migrate back the alignments
		rsync --dry-run --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/aligned_results/ ${ori_wkd_parent_dir}/${sample_batch}/aligned_results
		rsync --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/aligned_results/ ${ori_wkd_parent_dir}/${sample_batch}/aligned_results
		# Finally migrate back the gvcf files
		rsync --dry-run --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/gvcfs/ ${ori_wkd_parent_dir}/${sample_batch}/gvcfs
		rsync --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/gvcfs/ ${ori_wkd_parent_dir}/${sample_batch}/gvcfs
		# Migrate back the homo_region vcfs
		rsync --dry-run --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/vcfs/ ${ori_wkd_parent_dir}/${sample_batch}/vcfs
		rsync --rsync-path=/bin/rsync -avu --temp-dir="${tmp_dir}" --include="${sample_ID}*" --exclude="*" /staging/${seq_type}/${sample_batch}/vcfs/ ${ori_wkd_parent_dir}/${sample_batch}/vcfs
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": sample ID not correctly specified: sample_ID=${sample_ID}"
	fi
}

function check_raw_fastq(){
	local input_fq1=${1}
	local input_fq2=${2}
	local expected_lines=${3}

	if [[ -z ${expected_lines} ]]; then local expected_lines=80000000; fi
	local samp_ID=$(basename ${input_fq1} | awk -F '_1' '{printf "%s", $1;}')
	local rawf=$(dirname ${input_fq1})

	if check_fq_file_validity ${input_fq1} ${expected_lines} && check_fq_file_validity ${input_fq2} ${expected_lines}; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Raw FASTQ files ${input_fq1} and ${input_fq2} all valid."
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": ${input_fq1} or ${input_fq2} corrupted. Try to restore fastq file from old bams" && \
		if check_bam_validity ${rawf}/old_bam_incase_corrupt_fastq/${samp_ID}.bam 40021430; then
			module load BEDTools
			bedtools bamtofastq \
			-i ${rawf}/old_bam_incase_corrupt_fastq/${samp_ID}.bam \
			-fq ${rawf}/${samp_ID}_1.tmp.fastq \
			-fq2 ${rawf}/${samp_ID}_2.tmp.fastq && \
			gzip -f ${rawf}/${samp_ID}_1.tmp.fastq && \
			gzip -f ${rawf}/${samp_ID}_2.tmp.fastq && \
			mv ${rawf}/${samp_ID}_1.tmp.fastq.gz ${rawf}/${samp_ID}_1.fastq.gz && \
			mv ${rawf}/${samp_ID}_2.tmp.fastq.gz ${rawf}/${samp_ID}_2.fastq.gz
			module unload BEDTools
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": New pair of fastq file generated."; \
			ls -lh ${rawf}/${samp_ID}_1.fastq.gz
			ls -lh ${rawf}/${samp_ID}_2.fastq.gz
		else
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": No solid bam file available to replace the corrupted fastq files."
			return 1
		fi
	fi
}


function trim_adapter_and_fastqc(){
	module load TrimGalore/0.4.5
	module load cutadapt/2.3
	local samp_ID=${1}
	time trim_galore --paired --fastqc_args '--outdir=$fastqf' -o $trimf $rawf/"$1"_1.fastq.gz $rawf/"$1"_2.fastq.gz && echo "********** Trim + FastQC done ********"
	
	module unload TrimGalore/0.4.5
	module unload cutadapt/2.3

	local new_folderID
	local sub_folderID
	local uploaded_new_ID

	upload_gdrive_per_file_or_dir \
	-p ${sample_batch} \
	-i 1TiavZWvojroz7Oh2JibJTjuXJXGLy14B \
	-n new_folderID \
	-d dir
	upload_gdrive_per_file_or_dir \
	-p ${sample_batch}_fastq_QC \
	-i ${new_folderID} \
	-n sub_folderID \
	-d dir
	upload_gdrive_per_file_or_dir \
	-p ${trimf}/${sample_ID}_1_val_1_fastqc.html \
	-i ${sub_folderID} \
	-n uploaded_new_ID
	upload_gdrive_per_file_or_dir \
	-p ${trimf}/${sample_ID}_2_val_2_fastqc.html \
	-i ${sub_folderID} \
	-n uploaded_new_ID
}


function bwa_alignment(){
	if [[ -z ${ref_gen} ]]; then local ref_gen=/paedyl01/disk1/yangyxt/indexed_genome; fi
	if [[ ! -z ${2} ]]; then local trimf=${2}; fi
	if [[ ! -z ${3} ]]; then local bamf=${3}; fi
	if [[ -z ${4} ]]; then local ref_genome=${ref_gen}/ucsc.hg19.fasta; else local ref_genome=${4}; fi

	rm -f $bamf/$1.bam.tmp.*.bam || >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": No temp bam remained."
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Aligning paired end reads $trimf/$1_1_val_1.fq.gz $trimf/$1_2_val_2.fq.gz to ${ref_genome} now."
	if [[ ${seq_type} != wgs ]]
	then
		set -x 
		time bwa mem -t 10 -M -R "@RG\tID:$1\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:$1" ${ref_genome} \
		$trimf/"$1"_1_val_1.fq.gz $trimf/"$1"_2_val_2.fq.gz | samtools view -uSh -@ 8 - | samtools sort -O bam -@ 10 -o $bamf/$1.bam - && \
		samtools index $bamf/$1.bam && echo "BWA-MEM alignment done"
		set +x
	else
		set -x
		time bwa mem -t 10 -M -R "@RG\tID:$1\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:$1" ${ref_genome} \
		$trimf/"$1"_1_val_1.fq.gz $trimf/"$1"_2_val_2.fq.gz > $bamf/$1.sam && \
		time samtools sort -O bam -@ 8 -o $bamf/$1.bam < $bamf/$1.sam && \
		time samtools index $bamf/$1.bam && echo "BWA-MEM alignment done" && rm -f $bamf/$1.sam
		set +x
	fi
}


function bwa_align_contigs(){
	module load bwa
	
	local contig_fasta=${1}
	local output_sam=${2}
	local ID=$(echo ${contig_fasta} | awk -F '/' '{printf $NF}' | awk -F '.' '{for (i=1;i<NF;i++) printf $1"."}' | awk '{gsub(/\.$/, ""); print}')
	local ref_gen_ass=$ref_gen/ucsc.hg19.fasta
	time bwa mem -t 10 -M -R "@RG\tID:${ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${ID}" ${ref_gen_ass} \
	${contig_fasta} > ${output_sam}

	module unload bwa
}


function picard_mark_duplicates(){
	time java -XX:ConcGCThreads=1 -jar $PICARD MarkDuplicates \
	INPUT=$bamf/$1.bam \
	OUTPUT=$bamf/$1.deduped.bam \
	CREATE_INDEX=TRUE \
	REMOVE_DUPLICATES=TRUE \
	ASSUME_SORTED=TRUE \
	VALIDATION_STRINGENCY=LENIENT \
	METRICS_FILE=$bamf/$1.duplicatesindex && \
	samtools stats ${bamf}/${1}.deduped.bam > ${bamf}/${1}.deduped.stats && \
	echo "***** MarkDuplicates Done *****"
}


function spark_mark_duplicates(){
	mkdir ${bamf}/tmp_spark

	tmp_path=spark.local.dir=${bamf}/tmp_spark
	# Use new GATK tools(Picard implemented into GATK tool kit)
	time $gatk MarkDuplicatesSpark \
		-I $bamf/$1.bam \
		-O $bamf/$1.deduped.bam \
		-OBI true \
		--remove-all-duplicates true \
		-VS LENIENT \
		--num-executors 1 \
		--executor-cores 2 \
		--conf ${tmp_path} \
		-M $bamf/$1.duplicatesindex && \
	samtools stats ${bamf}/${1}.deduped.bam > ${bamf}/${1}.deduped.stats && \
	echo "***** MarkDuplicates Done *****"
}


function setNmMdUqTags(){
	time ${gatk} SetNmMdAndUqTags \
    -R ${ref_gen}/ucsc.hg19.fasta \
    -I ${bamf}/${1}.deduped.bam \
    -O ${bamf}/${1}.fixed.bam && \
	samtools stats ${bamf}/${1}.fixed.bam > ${bamf}/${1}.fixed.stats
}


function BQSR(){
	local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk

	time ${gatk} BaseRecalibrator \
	-R $ref_gen/ucsc.hg19.fasta \
	-I $bamf/"$1".fixed.bam \
	--known-sites $ref_gen/1000G_phase1.indels.hg19.sites.vcf \
	--known-sites $ref_gen/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	--known-sites $ref_gen/dbsnp_138.hg19.vcf \
	-O $bamf/"$1".bqsr.table \
	-LE TRUE && echo "***** Generating BQSR table done ******"

	# Gold_standard small indels 
	# Gold_standard indels from another project named 1000G
	# Gold standard SNV
	# Designates the output directory

	# The second step is to use ApplyBQSR to change the current quality scores

	time ${gatk} ApplyBQSR \
	-R $ref_gen/ucsc.hg19.fasta \
	-I $bamf/"$1".fixed.bam \
	-bqsr-recal-file $bamf/"$1".bqsr.table \
	-O $bamf/"$1".bqsr.bam && \
	samtools stats ${bamf}/${1}.bqsr.bam > ${bamf}/${1}.bqsr.stats && \
	echo "***BQSR DONE***"
}


function identify_gender_by_ploidy_on_Y(){
	local OPTIND c a p s e
    while getopts c::a::p:s:e: args
    do 
        case ${args} in
            a) local align_file=$OPTARG ;; 
            c) local count_file=$OPTARG ;;
            p) local probe=$OPTARG ;;
			s) local seq_type=$OPTARG ;;
			e) local ped_file=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

	local interval_file=/paedyl01/disk1/yangyxt/${seq_type}/healthy_bams_for_CNV/using_${probe}_probe/${probe}.hg19.preprocessed.interval_list
	local collect_read_script=/paedyl01/disk1/yangyxt/ngs_scripts/collect_read_counts.sh
	local ploidy_cohort_dir=/paedyl01/disk1/yangyxt/${seq_type}/healthy_bams_for_CNV/using_${probe}_probe/${probe}_ploidy_model
	local ploidy_subdir_prefix=${probe}_ploidy_normal_cohort
	local ploidy_cohort_model=${ploidy_cohort_dir}/${ploidy_subdir_prefix}-model

	if [[ -z ${count_file} ]] && [[ -f ${align_file} ]]; then
		bash ${collect_read_script} -p ${probe} -s ${align_file} -o $(dirname ${align_file}) -f hdf5
		local sample_ID=$(basename ${align_file} | awk -F '.' '{printf $1;}')
		local count_file=$(dirname ${align_file})/${sample_ID}.counts.hdf5
		local ploidy_case_dir=$(dirname ${align_file})/../CNV_calling/${probe}_ploidy_case_${sample_ID}
	elif [[ -f ${count_file} ]]; then
		local sample_ID=$(basename ${count_file} | awk -F '.' '{printf $1;}')
		local ploidy_case_dir=$(dirname ${count_file})/../CNV_calling/${probe}_ploidy_case_${sample_ID}
	fi

	local ploidy_case_subdir_prefix=${sample_ID}_ploidy_probe_${probe}

	source "$(which conda | awk -F '/' '{for (i=1;i<NF-1;i++) printf "%s/", $i;}')etc/profile.d/conda.sh"
	conda activate gcnv
	
	${gatk} DetermineGermlineContigPloidy \
	--model ${ploidy_cohort_model} \
	-I ${count_file} \
	-imr OVERLAPPING_ONLY \
	--output ${ploidy_case_dir} \
	--output-prefix ${ploidy_case_subdir_prefix} \
	--verbosity DEBUG

	conda deactivate

	# Extract the gender of this sample
	local Y_copy_number=$(awk -F '\t' '$1 == "chrY" || $1 == "Y"{printf $2;}' < ${ploidy_case_dir}/${ploidy_case_subdir_prefix}-calls/SAMPLE_O/contig_ploidy.tsv)
	if [[ ${Y_copy_number} -eq 1 ]]; then
		local gender=1	
	elif [[ ${Y_copy_number} -eq 0 ]]; then
		local gender=2
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): This sample does not have a normal chromosome Y copy number. Quit processing.
		cat ${ploidy_case_dir}/${ploidy_case_subdir_prefix}-calls/SAMPLE_O/contig_ploidy.tsv
		exit 1
	fi

	# Compare the gender with the registered gender in ped.
	local regist_gender=$(awk -F '\t' 'NR>1 && $2 == "'${sample_ID}'"{printf $6;}' < ${ped_file})
	if [[ ! -z ${regist_gender} ]] && [[ ${gender} == ${regist_gender} ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): This sample does have the gender prediction the same with ped table.
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): This sample does not have its gender prediction the same with ped table. Quit.
		exit 1
	fi
}


function parallel_run_haplotypecaller_unit(){
	# Here $1 should be the input file.
	local input_bam=${1}
	local input_bam_ID=$(echo ${1} | awk -F '/' '{printf $NF}' | awk -F '.' '{printf $1}')
	local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
	# Here $2 should be a path leading to a interval file specifying a region in the genome. 
	local special_region_interval=$2
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: The current calling chromosome is chr${3}
	if [[ ${special_region_interval} == "no_low_complex_interval" ]]
	then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "No special region, Directly run HC"
		time $gatk --java-options "-Xmx8G" HaplotypeCaller \
		--emit-ref-confidence GVCF \
		-R ${ref_gen}/ucsc.hg19.fasta \
		-I ${input_bam} \
		-O ${gvcf}/${input_bam_ID}.chr${3}.HC.g.vcf.gz \
		-bamout $bamf/${input_bam_ID}.chr${3}.realigned.bam \
		--assembly-region-out $gvcf/${input_bam_ID}.chr${3}.HC.active_region.tsv \
		--active-probability-threshold 0.0012 \
		--assembly-region-padding 100 \
		--bam-writer-type CALLED_HAPLOTYPES \
		--linked-de-bruijn-graph \
		-L chr${3} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-OBI true && echo '** GVCF $1.raw.g.vcf.gz done **'
	else
		local interval_file_name=$(echo ${special_region_interval} | awk -F '/' '{printf $NF}' | awk -F '.' '{printf $1}')
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The interval_file_name is "${interval_file_name}

		# Run the HC without the special region using -XL argument in normal mode. Then run the 
		time $gatk --java-options "-Xmx8G" HaplotypeCaller \
		--emit-ref-confidence GVCF \
		-R $ref_gen/ucsc.hg19.fasta \
		-I ${input_bam} \
		-O ${gvcf}/${input_bam_ID}.chr${3}.without_${interval_file_name}.HC.g.vcf.gz \
		-bamout ${bamf}/${input_bam_ID}.chr${3}.without_${interval_file_name}.realigned.bam \
		--active-probability-threshold 0.0012 \
		--assembly-region-out ${gvcf}/${input_bam_ID}.chr${3}.HC.active_region.tsv \
		--assembly-region-padding 100 \
		--bam-writer-type CALLED_HAPLOTYPES \
		--linked-de-bruijn-graph \
		-XL ${special_region_interval} \
		-L chr${3} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-OBI true
	fi
}


function parallel_run_HC(){
	module load parallel/20210122
	source /home/yangyxt/software/parallel_20210122/env_parallel.bash

	local input_bam=${1}
	local input_bam_ID=$(echo ${1} | awk -F '/' '{printf $NF}' | awk -F '.' '{printf $1}')
	local gatk=/home/yangyxt/software/gatk-4.2.5.0/gatk
	local special_interval=${2}

	trap 'catch_error_report ${?} ${LINENO} ${BASH_COMMAND}' ERR

	check_bam_index ${input_bam}
	# export -f parallel_run_haplotypecaller_unit

	local -a contigs=(M 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
	if [[ ${seq_type} != "wgs" ]]; then
		local threads=$(determine_job_num 15 1)
	else
		local threads=$(determine_job_num 50 1)
	fi

	find_involved_vars "parallel_run_haplotypecaller_unit" __involved_vars
    >&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": ENV/SHELL variables involved in function parallel_run_haplotypecaller_unit are: "${__involved_vars}
    local parallel_args=$(echo ${__involved_vars} | awk 'BEGIN{RS=" ";} length($1) > 1{printf "--env %s ", $1;}')

	# backup code for redirect output log content '>' ${gvcf}/${input_bam_ID}_HC_calling_on_chr{}.log '2>&1' '&&' chmod 664 ${gvcf}/${input_bam_ID}_HC_calling_on_chr{}.log 
	env_parallel ${parallel_args}--dry-run --joblog ${gvcf}/${input_bam_ID}.HC.allchr.parallel.log --keep-order --jobs ${threads} parallel_run_haplotypecaller_unit ${input_bam} ${special_interval} {} '>' ${gvcf}/${input_bam_ID}.HC.{}.log '2>&1' ::: $(echo ${contigs[*]}) && \
	env_parallel ${parallel_args}--keep-order --joblog ${gvcf}/${input_bam_ID}.HC.allchr.parallel.log --jobs ${threads} parallel_run_haplotypecaller_unit ${input_bam} ${special_interval} {} '>' ${gvcf}/${input_bam_ID}.HC.{}.log '2>&1' ::: $(echo ${contigs[*]})
	check_parallel_joblog ${gvcf}/${input_bam_ID}.HC.allchr.parallel.log
	# merge_logs $(echo "${contigs[*]}" | awk '{for (i=1;i<NF;i++) {printf "'${gvcf}'/'${input_bam_ID}'_HC_calling_on_chr"$i".log,";} printf "'${gvcf}'/'${input_bam_ID}'_HC_calling_on_chr"$NF".log"}') ${gvcf}/${input_bam_ID}_HC_calling_on_all_chrs.log && \
	# chmod 664 ${gvcf}/${input_bam_ID}_HC_calling_on_all_chrs.log

	if [[ ${special_interval} == "no_low_complex_interval" ]]; then
		local tobe_merged=$(echo "${contigs[*]}" | awk '{for (i=1;i<NF;i++) {printf "'${gvcf}'/'${input_bam_ID}'.chr"$i".HC.g.vcf.gz,";} printf "'${gvcf}'/'${input_bam_ID}'.chr"$NF".HC.g.vcf.gz";}')
		local merged=${gvcf}/${input_bam_ID}.HC.g.vcf
		local gather_script=${batchf}/awkgen_gathervcfs_GVCF.sh
		local tobe_merged_bams=$(echo "${contigs[*]}" | awk '{for (i=1;i<NF;i++) {printf "'${bamf}'/'${input_bam_ID}'.chr"$i".realigned.bam,";} printf "'${bamf}'/'${input_bam_ID}'.chr"$NF".realigned.bam";}')
		
		mergebams ${tobe_merged_bams} ${bamf}/${input_bam_ID}.realigned.bam yes

		# Check file existence
		local -a test_existed=($(echo ${tobe_merged} | awk 'BEGIN{RS=",";} {printf "%s ", $1;}'))
		for item in "${test_existed[@]}"; do
			ls -lh ${item} || echo "No sample existed naming ${item}"
		done

		bcftools_concatvcfs ${tobe_merged} ${merged} ${gather_script} && \
		bgzip -f ${gvcf}/${input_bam_ID}.HC.g.vcf && tabix -f -p vcf ${gvcf}/${input_bam_ID}.HC.g.vcf.gz && \
		rm -f ${gvcf}/${input_bam_ID}.chr*.HC.g.vcf*
	else
		local interval_file_name=$(echo ${special_interval} | awk -F '/' '{printf $NF}' | awk -F '.' '{printf $1}')
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: "The interval_file_name is "${interval_file_name}

		time $gatk --java-options "-Xmx75G" HaplotypeCaller \
		--emit-ref-confidence GVCF \
		-R $ref_gen/ucsc.hg19.fasta \
		-I ${input_bam} \
		-O $gvcf/${input_bam_ID}.only_${interval_file_name}.HC.g.vcf.gz \
		-bamout $bamf/${input_bam_ID}.only_${interval_file_name}.realigned.bam \
		--force-active \
		--allow-non-unique-kmers-in-ref \
		--debug-assembly \
		--kmer-size 21 \
		--max-num-haplotypes-in-population 256 \
		--assembly-region-padding 100 \
		--bam-writer-type CALLED_HAPLOTYPES \
		--linked-de-bruijn-graph \
		-L ${special_interval} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-imr OVERLAPPING_ONLY \
		-OBI true

		local tobe_merged=$(echo "${contigs[*]}" | awk '{for (i=1;i<=NF;i++) {printf "'${gvcf}'/'${input_bam_ID}'.chr"$i".without_'${interval_file_name}'.HC.g.vcf.gz,";}} END{printf "'${gvcf}'/'${input_bam_ID}'.only_'${interval_file_name}'.HC.g.vcf.gz";}')
		local merged=${gvcf}/${input_bam_ID}.HC.g.vcf
		local gather_script=${batchf}/awkgen_gathervcfs_GVCF_${interval_file_name}.sh
		local tobe_merged_bams=$(echo "${contigs[*]}" | awk '{for (i=1;i<=NF;i++) {printf "'${bamf}'/'${input_bam_ID}'.chr"$i".without_'${interval_file_name}'.realigned.bam,";}} END{printf "'${bamf}'/'${input_bam_ID}'.only_'${interval_file_name}'.realigned.bam";}')
		
		mergebams ${tobe_merged_bams} ${bamf}/${input_bam_ID}.realigned.bam yes
		
		bcftools_concatvcfs ${tobe_merged} ${merged} ${gather_script} && \
		bgzip -f ${gvcf}/${input_bam_ID}.HC.g.vcf && tabix -f -p vcf ${gvcf}/${input_bam_ID}.HC.g.vcf.gz && \
		rm -f ${gvcf}/${input_bam_ID}.chr*.HC.g.vcf*
	fi

	module unload parallel/20210122
}


function call_shortV_on_homo_regions() {
	local input_bam=${1}
	local bamf=$(dirname ${input_bam})
	local sample_ID=$(basename ${input_bam} | awk -F '.' '{printf "%s", $1;}')
	local -a homologous_regions=($(ls -d /paedyl01/disk1/yangyxt/indexed_genome/homologous_regions/*/*.bed))
	local -a hr_names
	for hr in "${homologous_regions[@]}"; do hr_names+=( $(basename ${hr} | awk '{name = gensub(/(.*)_related.*/, "\\1", "g", $0); printf "%s", name;}') ); done
	local tmp_map_file=${gvcf/staging/paedyl01\/disk1\/yangyxt}/${sample_ID}_homo_region_map.tsv
	if [[ ! -f ${tmp_map_file} ]]; then touch ${tmp_map_file}; else : > ${tmp_map_file}; fi
	
	# Prepare argument map file for following parallel 
	module load parallel

	if [[ ${seq_type} != "wgs" ]]; then
		local threads=$(determine_job_num 20 1)
	else
		local threads=$(determine_job_num 120 1)
	fi

	parallel --joblog ${bamf}/parallel_prepare_map_file.${sample_ID}.log --link --dry-run -j${threads} bash ${scrf/staging/paedyl01\/disk1\/yangyxt}/GATK4.1_based_NGS_preprocessing_script.sh prepare_map_file_per_region {1} ${input_bam} ${tmp_map_file} '>' ${bamf}/${sample_ID}.bqsr.{2}.genmap.log '2>&1' ::: "${homologous_regions[@]}" ::: "${hr_names[@]}" && \
	parallel --joblog ${bamf}/parallel_prepare_map_file.${sample_ID}.log --link -j${threads} bash ${scrf/staging/paedyl01\/disk1\/yangyxt}/GATK4.1_based_NGS_preprocessing_script.sh prepare_map_file_per_region {1} ${input_bam} ${tmp_map_file} '>' ${bamf}/${sample_ID}.bqsr.{2}.genmap.log '2>&1' ::: "${homologous_regions[@]}" ::: "${hr_names[@]}" && \
	cat ${bamf}/parallel_prepare_map_file.${sample_ID}.log
	cat ${tmp_map_file::-4}.*.tsv > ${tmp_map_file}.tmp && \
	mv ${tmp_map_file}.tmp ${tmp_map_file} && \
	rm -f ${tmp_map_file::-4}.*.tsv

	>&2 echo $'\n'"Here is the tmp_map_file after prepare_map_file_per_region: ${tmp_map_file}"
	ls -lh ${tmp_map_file} && \
	cat ${tmp_map_file} && \
	>&2 echo $'\n'

	if [[ $(cat ${tmp_map_file} | wc -l) -eq 0 ]]; then
		>&2 echo "The ${tmp_map_file} seems empty. Abandon calling homologous gene for sample ${sample_ID}, the ${input_bam} seems not overlapping with any homologous regions."
		return 1
	else
		if [[ ${seq_type} != "wgs" ]]; then
			local threads=$(determine_job_num 22 1)
		else
			local threads=$(determine_job_num 60 1)
		fi

		parallel --dry-run --link -j${threads} \
		--joblog ${bamf}/call_homo_region_vcf.${sample_ID}.log \
		bash ${central_scripts}/realign_masked_and_HC_multiploidy.sh \
		call_polyploidy_per_priority_region \
		-s {2} -t {1} -f {4} -r {5} \
		-a ${input_bam} \
		-m {3} '>' ${bamf}/${sample_ID}.bqsr.{6}.recall.log '2>&1' ::: \
		$(awk -F '\t' '{printf "%s ", $1;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $2;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $3;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $5;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $6;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s\n", $2;}' < ${tmp_map_file} | awk -F '/' '{printf "%s ", $NF;}') && \
		parallel --link -j${threads} --group \
		--joblog ${bamf}/call_homo_region_vcf.${sample_ID}.log \
		bash ${central_scripts}/realign_masked_and_HC_multiploidy.sh \
		call_polyploidy_per_priority_region \
		-s {2} -t {1} -f {4} -r {5} \
		-a ${input_bam} \
		-m {3} '>' ${bamf}/${sample_ID}.bqsr.{6}.recall.log '2>&1' ::: \
		$(awk -F '\t' '{printf "%s ", $1;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $2;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $3;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $5;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s ", $6;}' < ${tmp_map_file}) ::: \
		$(awk -F '\t' '{printf "%s\n", $2;}' < ${tmp_map_file} | awk -F '/' '{printf "%s ", $NF;}'); \
		check_parallel_joblog ${bamf}/call_homo_region_vcf.${sample_ID}.log | \
		awk 'NR == FNR {arr[$5] = $5;} NR > FNR{if ($2 in arr) print;}' - ${tmp_map_file} > ${tmp_map_file}.tmp

		if [[ $(cat ${tmp_map_file}.tmp | wc -l) -gt 0 ]]; then
			>&2 echo $'\n'"Remove the failed records in ${tmp_map_file}.tmp"
			mawk -F '\t' 'NR == FNR {arr[$0] = $0;} NR > FNR {if ($0 in arr) {next;} else print;}' ${tmp_map_file}.tmp ${tmp_map_file} > ${tmp_map_file}.new && \
			mv ${tmp_map_file}.new ${tmp_map_file} && rm -f ${tmp_map_file}.tmp
		else
			>&2 echo "All parallel jobs finished properly, remove ${tmp_map_file}.tmp, the covered regions are:"
			>&2 awk -F '\t' '{printf "%s ", $2;}' < ${tmp_map_file}
			rm -f ${tmp_map_file}.tmp
		fi
			
		>&2 echo $'\n'"Here is the tmp_map_file after variants calling: ${tmp_map_file}"
		ls -lh ${tmp_map_file} && \
		cat ${tmp_map_file}
	fi
	check_return_code
	module unload parallel

	# concat variant records on all homo regions (per sample)
	local -a vcf_candidates=($(cut -f 4 ${tmp_map_file} | uniq - | awk '{printf "%s ", $1;}'))
	local valid_vcfs=${tmp_map_file/.tsv/.valid.dip.vcf.lst}
	local invalid_vcfs=${tmp_map_file/.tsv/.invalid.dip.vcf.lst}
	: > ${valid_vcfs}
	: > ${invalid_vcfs}
	for vc in "${vcf_candidates[@]}"; do  # First, these to be concated vcf files should be gzipped
		if check_vcf_validity ${vc} 1 && [[ -f ${vc} ]]; then
			awk 'BEGIN{print "'${vc}'";}' >> ${valid_vcfs}
		else
			awk 'BEGIN{print "'${vc}'";}' >> ${invalid_vcfs}
		fi
	done
	cat ${valid_vcfs} | sort - | uniq - > ${valid_vcfs}.tmp && mv ${valid_vcfs}.tmp ${valid_vcfs}
	cat ${invalid_vcfs} | sort - | uniq - > ${invalid_vcfs}.tmp && mv ${invalid_vcfs}.tmp ${invalid_vcfs}
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Before concating vcfs, $(wc -l ${valid_vcfs} | awk '{print $1;}') of them are valid, $(wc -l ${invalid_vcfs} | awk '{print $1;}') of them are invalid."
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The invalid vcfs are: $(awk '{printf "%s ", $1;}' ${invalid_vcfs})"
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The valid vcfs are: $(awk '{printf "%s ", $1;}' ${valid_vcfs})"

	local tobe_concat=$(cat ${valid_vcfs} | uniq - | awk -F '\t' '{printf "%s,", $1;}')
	bash ${central_scripts}/common_bash_utils.sh bcftools_concatvcfs \
	${tobe_concat::-1} \
	${vcf}/${sample_ID}.homo_region.vcf.gz && \
	cut -f 1,2,4 ${tmp_map_file} > ${vcf}/${sample_ID}.homo_regions.vcf.lst && \
	rm -f $(echo ${tobe_concat::-1} | awk 'BEGIN{RS=",";} {printf "%s ", $1;}')

	# remove false positive on intrinsic variants
	# if [[ -f ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz ]]; then
	# 	rm -f ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz && \
	# 	touch ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz && \
	# 	chmod 664 ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz
	# else
	# 	touch ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz && \
	# 	chmod 664 ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz
	# fi

	python3 ${central_scripts}/compare_with_intrinsic_vcf.py -qv ${vcf}/${sample_ID}.homo_region.vcf.gz -ov ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz && \
	tabix -f -p vcf ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz
}


function prepare_map_file_per_region() {
	local region_bed=${1}
	local input_bam=${2}
	local wkd=$(dirname ${input_bam})
	local sample_ID=$(basename ${input_bam} | awk -F '.' '{printf "%s", $1;}')
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The sample ID is ${sample_ID}"
	local region_name=$(basename ${region_bed} | awk -F '.' '{printf $1;}')
	local tmp_map_file=${gvcf/staging/paedyl01\/disk1\/yangyxt}/${sample_ID}_homo_region_map.${region_name}.tsv
	if [[ ! -f ${tmp_map_file} ]]; then touch ${tmp_map_file}; else : > ${tmp_map_file}; fi
	local ref_assembly=/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta

	check_bam_validity ${input_bam}
	check_return_code

	local region_bed_dir=$(dirname ${region_bed})
	local priority_components=$(ls -d ${region_bed_dir}/*/ | awk -F '/' '{printf "%s%s.bed,", $0, $(NF-1);}' | awk '{gsub(/,$/, ""); print;}')
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Priority components for this set of homo regions are ${priority_components}"
	local -a pc_arr=($(echo ${priority_components} | awk 'BEGIN{RS=",";} {printf "%s ", $1;}'))
	
	local forward_read_varname=${wkd}/${sample_ID}.${region_name}-XA_1.fastq
	local reverse_read_varname=${wkd}/${sample_ID}.${region_name}-XA_2.fastq

	# IF I run this function by running the script, this will cause the function ran in a subshell
	# Which will fail to pass any value back by eval command
	bash ${central_scripts}/realign_masked_and_HC_multiploidy.sh prepare_fastq_files \
	-t ${region_bed} \
	-a ${input_bam} \
	-f ${wkd}/${sample_ID}_read1.fq \
	-r ${wkd}/${sample_ID}_read2.fq

	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Here is the extracted forward reads ${forward_read_varname} for ${region_bed}"
	>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Here is the extracted reverse reads ${reverse_read_varname} for ${region_bed}"
	ls -lh ${forward_read_varname}
	ls -lh ${reverse_read_varname}

	if [[ $(wc -l ${forward_read_varname} | awk '{print $1;}') -lt 4 ]]; then
		echo "WARNING!!!Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Failed to fetch any multi-aligned reads from ${input_bam} in region ${region_bed}. Skip this region and now deleting these extracted fastq files." && \
		rm ${forward_read_varname} && rm ${reverse_read_varname} && true
	else
		for pc in "${pc_arr[@]}"; do
			local ori_depth=$(extract_ave_depth_without_MQ ${input_bam} ${pc} | awk '{printf "%.0f", $1;}')
			if [[ ${ori_depth} -le 0 ]]; then 
				>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This priority region ${pc} does not have reads covered. Skip this region for preparing map file."
				continue
			fi
			local pc_name=$(basename ${pc} | awk -F '.' '{printf "%s", $1;}')
			local complement_bed=$(dirname ${pc})/$(basename ${ref_assembly} | awk -F '.' '{printf $1"."$2;}').sub.$(basename ${pc} | awk -F '.' '{printf $1;}').bed
			local masked_genome=$(dirname ${complement_bed})/$(basename ${ref_assembly} | awk -F '.' '{printf "%s.%s", $1, $2;}').$(basename ${complement_bed} | awk -F '.' '{for (i=1;i<NF;i++) printf "%s.", $i;}')masked.fasta
			mkdir /paedyl01/disk1/yangyxt/indexed_genome/masked_genomes/${pc_name} || echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": /paedyl01/disk1/yangyxt/indexed_genome/masked_genomes/${pc_name} folder already created."
			# Prepare masked genome.
			bash ${central_scripts}/realign_masked_and_HC_multiploidy.sh prepare_masked_genome \
			-g ${ref_assembly} \
			-s ${pc} \
			-t ${region_bed} && \
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": This is the masked genome: $(ls -lh ${masked_genome})"
			local output_gvcf=$(dirname ${input_bam})/${sample_ID}.only_${pc_name}.dip.vcf.gz
			awk -v wr="${region_bed}" \
			-v pc="${pc}" \
			-v mg="${masked_genome}" \
			-v og="${output_gvcf}" \
			-v fr="${forward_read_varname}" \
			-v rr="${reverse_read_varname}" \
			'BEGIN{printf "%s\t%s\t%s\t%s\t%s\t%s\n", wr, pc, mg, og, fr, rr;}' >> ${tmp_map_file} && \
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Now the ${tmp_map_file} looks like:"
			cat ${tmp_map_file}
		done
	fi
}


function check_homo_region_called() {
	local map_file=${1}
	local -a homo_regions=($(ls -d /paedyl01/disk1/yangyxt/indexed_genome/homologous_regions/*/*.bed))
	
	local pcs  # pcs stands for principle components
	for hr in "${homo_regions[@]}"; do
		local pcs="${pcs}$(ls -d $(dirname ${hr})/*/ | awk -F '/' '{printf "%s%s.bed ", $0, $(NF-1);}')"
	done

	local -a pc_arr=($(echo ${pcs} | awk 'BEGIN{RS=" ";} { gsub(/\n$/, ""); print $1;}' | sort - | uniq -))
	local -a map_w=($(cut -f 1 ${map_file} | sort - | uniq - | awk '{printf "%s ", $1;}'))
	local -a map_p=($(cut -f 2 ${map_file} | sort - | uniq - | awk '{printf "%s ", $1;}'))

	local -a left_regions=($(bash ${central_scripts}/return_array_substraction "${map_w[*]}" "${homo_regions[*]}"))
	local -a left_pcs=($(bash ${central_scripts}/return_array_substraction "${map_p[*]}" "${pc_arr[*]}"))

	if [[ ${#left_regions[@]} -eq 0 ]] && [[ ${#left_pcs[@]} -eq 0 ]]; then
		true
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Here are the total pcs: ${map_p[*]}"$'\n'"Total whole_regions: ${map_w[*]}"$'\n'"pcs from map ${pc_arr[*]}"$'\n'"whole regions from map ${homo_regions[*]}"
		return 1
	fi
}


function main_variant_calling(){
	set -E
	trap 'catch_error_report ${?} ${LINENO} ${BASH_COMMAND}' ERR
	
	migrate_to_staging_per_sample

	# Check whether raw fastq file is solid
	check_raw_fastq ${rawf}/${sample_ID}_1.fastq.gz ${rawf}/${sample_ID}_2.fastq.gz

	# migrate_to_staging_per_sample
	# Before starting the whole process, we need to trim the adapter and use FastQC to inspect the results. While doing this, using option'--fastqc' will have another round of validation using fastqc after the trimming.
	if [[ -f ${trimf}/${sample_ID}_1_val_1.fq.gz ]] && [[ -f ${trimf}/${sample_ID}_2_val_2.fq.gz ]]; then
		if check_fq_file_validity ${trimf}/${sample_ID}_1_val_1.fq.gz 88533928 && check_fq_file_validity ${trimf}/${sample_ID}_2_val_2.fq.gz 88533928; then
			>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Trimmed fastq files already generated. Skip the TrimGalore execution."$'\n\n'
		else
			trim_adapter_and_fastqc ${sample_ID}
			check_return_code
		fi
	else
		trim_adapter_and_fastqc ${sample_ID}
		check_return_code
	fi

	# First, we need to do alignment (read-paris against the reference genome you picked), which uses BWA-- the most widely used software for alignment
	# To use BWA, you need to first index the genome with `bwa index'. Ususally we pre-indexed reference genome, and then use them here directly. 
	# There are three alignment algorithms in BWA: 1. `mem', 2. `bwasw', and 3. `aln/samse/sampe'. Here "mem" means we chose BWA-MEM algorithm

	# -t，用于设定排序时的线程数，我们可以设置为10
	# -M Mark shorter split hits as secondary (for Picard compatibility).
	# -R 用于设定Read Group的信息（每项之间必须用制表符\t分开）：ID设置为测序的lane的ID， PL设置为测序平台信息， SM样本ID信息， LB测序文库的名字（重要性稍低）。 不同组之间的测序过程被认为是相互独立的，这个信息对于后续MarkDuplicates非常重要。

	# The subsequent command align read pairs to ref first, then transfer the SAM to BAM and sort the BAM in the order of genome coordinates
	# samtools is one tool developed by the developer of BWA, which is used to help process the data before variant Calling and after alignment.
	# samtools view is used to transform SAM to BAM. While -b means output BAM files, -u means output uncompressed BAM, -S regulates the input file is automatically detected. -@ 10 regulates the threads number. "-" in the end indicates the position of input file from previous pipeline.
	# samtools sort is used to resort the reads in the order of genome coordinates. -o is used to regulates the passway of output files. -@ is used to regulates threads number.
	if [[ -f ${bamf}/${sample_ID}.bam ]] && check_bam_validity ${bamf}/${sample_ID}.bam; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Raw BAM file already generated. Skip the BWA_MEM execution now."$'\n\n'
	elif [[ ${seq_type} != "wgs" ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Raw BAM file not existed or corrupted. Redo BWA_MEM execution now."$'\n\n'
		bwa_alignment ${sample_ID}
		check_return_code
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Raw BAM file not existed or corrupted. And this is WGS sample Redo mapping with DragMap"$'\n\n'
		independent_dragmap_alignment \
		-f ${trimf}/${sample_ID}_1_val_1.fq.gz \
		-r ${trimf}/${sample_ID}_2_val_2.fq.gz \
		-o ${bamf}/${sample_ID}.bqsr.bam \
		-t 10
	fi

	# Then wat we need to do is marking duplicates
	# This line regulates that the after duplicates marking, we should build an index for this bam file for further local realignment, it's necessary.
	# This line regulates during duplicates marking, we also remove those duplicated reads 
	# This is set to LENIENT incase SAM validation error leading to program failure
	# This line designates the file to store all the marked duplicates
	if [[ -f ${bamf}/${sample_ID}.bqsr.bam ]] && check_bam_validity ${bamf}/${sample_ID}.bqsr.bam; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": BQSR BAM file already generated. Skip the MarkDuplicate execution now."$'\n\n'
	else
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": BQSR BAM not existed or corrupted. Do MarkDuplicate execution now."$'\n\n'
		picard_mark_duplicates ${sample_ID}
		check_return_code
	fi
	# Then we need to do the local realignment because in BWA algorithm, the scoring matrix ultilized does support multiple mismatches and short gaps rather than a long gap. This leads to the false negative detection of many relatively longer indels.
	# The realignment algorithm we used is the well-known Smith-Waterman local realignment algorithm
	# In GATK4, since using HaplotypeCaller for Raw-variant calling is universal, the Broad Institute no longer provide Local Realignment tool in the stage of Data preprocessing. 

	# MarkDuplicatesSpark processing can replace both the MarkDuplicates and SortSam steps of the Best Practices single sample pipeline. 
	# After flagging duplicate sets, the tool automatically coordinate-sorts the records. 
	# It is still necessary to subsequently run SetNmMdAndUqTags before running BQSR.
	if [[ -f ${bamf}/${sample_ID}.bqsr.bam ]] && check_bam_validity ${bamf}/${sample_ID}.bqsr.bam; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": BQSR BAM file already generated. Skip the GATK SetNmMdAndUqTags execution now."$'\n\n'
	else
		setNmMdUqTags ${sample_ID}
		check_return_code
	fi

	# Then we should perform Base Quality Recalibration
	# 在WGS分析中，变异检测是一个极度依赖测序碱基质量值的步骤。因为这个质量值是衡量我们测序出来的这个碱基到底有多正确的重要（甚至是唯一）指标。
	# 如果我们在看到某一个碱基报告的质量值是20时，那么它的预期错误率是1%，反过来想，就等于是说如果有100个质量值都是20的碱基，那么从统计上讲它们中将只有1个是错的！做了这个等效变换之后，我们的问题就可以转变成为寻找错误碱基的数量了。
	# 首先排除掉所有的已知变异位点，然后计算每个（报告出来的）质量值下面有多少个碱基在比对之后与参考基因组上的碱基是不同的，这些不同碱基就被我们认为是错误的碱基，它们的数目比例反映的就是真实的碱基错误率，换算成Phred score。
	if [[ -f ${bamf}/${sample_ID}.bqsr.bam ]] && check_bam_validity ${bamf}/${sample_ID}.bqsr.bam; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": BQSR BAM file already generated. Skip the GATK BQSR execution now."$'\n\n'
	else
		BQSR ${sample_ID}
		check_return_code
	fi

	if [[ -f ${bamf}/${sample_ID}.bqsr.bam ]] && check_bam_validity ${bamf}/${sample_ID}.bqsr.bam && [[ -L ${bamf}/${sample_ID}.bqsr.fade.bam ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": BQSR BAM file already generated. And its been processed by FADE to remove sequencing artifacts."
	else
		:
		# fade_per_bam -i ${bamf}/${sample_ID}.bqsr.bam -f replace
		# check_return_code
		# ln -s ${bamf}/${sample_ID}.bqsr.bam ${bamf}/${sample_ID}.bqsr.fade.bam
	fi

	# We need to call the variants on homologous sequences.
	if check_vcf_validity ${vcf}/${sample_ID}.homo_region.filtered.vcf.gz && [[ -f ${vcf}/${sample_ID}.homo_regions.vcf.lst ]] && check_homo_region_called ${vcf}/${sample_ID}.homo_regions.vcf.lst ; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Variants on homologous regions have already been called without a problem. Skip this step now."$'\n\n'
	else
		call_shortV_on_homo_regions ${bamf}/${1}.bqsr.bam
		check_return_code
	fi
	
	# Generated ${vcf}/${1}.homo_region.vcf.gz
	# We are going to save the Mapping qualities of the bqsr.bam file
	if [[ -L ${bamf}/${1}.saved.bam ]]; then
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Sample ${1} bam file has already been processed with saveMQ script."$'\n\n'
	else
		# /usr/bin/time bash ${central_scripts}/grab_multisplit_aligned_reads.sh ${bamf}/${1}.bqsr.bam || echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Failed saving MQ for this bam."$'\n\n'
		# if check_bam_validity ${bamf}/${1}.bqsr.saved.bam; then
		# 	# You need to clear old index files to make sure downstream tools use the latest one.
		# 	{ samtools stats ${bamf}/${1}.bqsr.saved.bam > ${bamf}/${1}.bqsr.saved.stats && \
		# 	mv ${bamf}/${1}.bqsr.saved.bam ${bamf}/${1}.bqsr.bam && \
		# 	{ rm -f ${bamf}/${1}.bqsr.saved.bai ${bamf}/${1}.bqsr.saved.bam.bai ${bamf}/${1}.bqsr.bai ${bamf}/${1}.bqsr.bam.bai || echo "All bam index file cleared."; } && \
		# 	ln -s ${bamf}/${1}.bqsr.bam ${bamf}/${1}.saved.bam && \
		# 	samtools index ${bamf}/${1}.bqsr.bam; } || echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Failed saving MQ for this bam."$'\n\n'
		# else
		# 	echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": It seems the MQ saved bam ${bamf}/${1}.bqsr.saved.bam is just not valid. Quit replacing it with original bqsr bam and delete it now."
		# 	rm -f ${bamf}/${1}.bqsr.saved.ba*
		# fi
		>&2 echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": Abandon executing saveMQ script."$'\n\n'
	fi

	# After saving the MQ values. 
	# Remove all the temporary files
	rm ${bamf}/${1}.deduped.ba* 2>/dev/null || echo "${bamf}/${1}.deduped.ba* files not existed."
	rm ${bamf}/${1}.fixed.ba* 2>/dev/null || echo "${bamf}/${1}.fixed.ba* files not existed."

	# Ready for HaplotypeCaller #
	# Call raw gvcf
	# run_haplotype_caller ${sample_ID} ${special_interval}
	parallel_run_HC ${bamf}/${sample_ID}.bqsr.bam ${special_interval}
	check_return_code

	ls -lh ${gvcf}/${sample_ID}.*HC.g.vcf.gz || { echo "No samples existed matching pattern ${gvcf}/${sample_ID}.*HC.g.vcf.gz" && exit 1; }
	# Merge the realigned part bam file with the original bam file
	# time java -jar ${PICARD} MergeSamFiles \
	# I=$bamf/"$1".realigned.bam \
	# I=$bamf/"$1".bqsr.bam \
	# O=$bamf/"$1".merged.bam
	migrate_from_staging_per_sample
}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
	declare -a func_names=($(typeset -f | awk '!/^main[ (]/ && /^[^ {}]+ *\(\)/ { gsub(/[()]/, "", $1); printf "%s ", $1;}'))
    declare -a input_func_names=($(return_array_intersection "${func_names[*]}" "$*"))
    if [[ ${#input_func_names[@]} -gt 0 ]]; then
        "$@"
    else
		set -T
		set -E
        trap 'catch_error_report ${?} ${LINENO} ${BASH_COMMAND}' ERR
		trap 'echo "$BASH_COMMAND" > /dev/null' DEBUG
		trap 'catch_exit_status ${?}' EXIT
        main_variant_calling "$@"
    fi
fi