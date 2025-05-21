#!/bin/bash
set -e -o pipefail # Exit on error, and treat pipe errors correctly

# Logging function
# Usage: log "My log message"
# Writes to current stderr. If stderr of the calling function is redirected,
# the log output will go to that redirected location.
log() {
    local msg="$1"
    local call_info line_num func_name script_name timestamp
    call_info=($(caller 0)) # Get info about the call site: LINE FUNCNAME FILENAME
    line_num="${call_info[0]}"
    func_name="${call_info[1]}"
    script_name="${call_info[2]##*/}" # Basename of the script file
    timestamp=$(date +"%Y-%m-%d %H:%M:%S.%3N")
    >&2 echo "[$timestamp] [Script: $script_name] [Func: $func_name] [Line: $line_num] :: $msg"
}


# Worker function to process a single line from the BED file.
# stderr for this entire function's execution will be redirected.
_process_single_bed_line() {
    # --- Initial Argument Reception ---
    local p_chr1_raw="$1"             # BED col 1
    local p_start1_raw="$2"           # BED col 2
    local p_end1_raw="$3"             # BED col 3
    local p_original_strand1_raw="$4" # BED col 7 (Strand for interval 1)
    local p_chr2_raw="$5"             # BED col 4 (Chromosome for interval 2)
    local p_start2_raw="$6"           # BED col 5 (Start for interval 2)
    local p_end2_raw="$7"             # BED col 6 (End for interval 2)
    local p_original_strand2_raw="$8" # BED col 8 (Strand for interval 2)
    local fasta_file="$9"
    local worker_stderr_log_dir="${10}" # Directory to store this worker's stderr log
    local worker_data_parts_dir="${11}" # Directory to store this worker's data output part

    # --- Sanitize Input Parameters ---
    local p_chr1=$(echo "$p_chr1_raw" | tr -d '[:space:]')
    local p_start1=$(echo "$p_start1_raw" | tr -d '[:space:]')
    local p_end1=$(echo "$p_end1_raw" | tr -d '[:space:]')
    local p_original_strand1=$(echo "$p_original_strand1_raw" | tr -d '[:space:]')
    local p_chr2=$(echo "$p_chr2_raw" | tr -d '[:space:]')
    local p_start2=$(echo "$p_start2_raw" | tr -d '[:space:]')
    local p_end2=$(echo "$p_end2_raw" | tr -d '[:space:]')
    local p_original_strand2=$(echo "$p_original_strand2_raw" | tr -d '[:space:]')


    # --- Setup Job-Specific Naming and Logging ---
    local s_p_chr1 s_p_chr2 coord_prefix unique_suffix worker_stderr_file job_identifier
    s_p_chr1=${p_chr1//[^a-zA-Z0-9_.-]/_} 
    s_p_chr2=${p_chr2//[^a-zA-Z0-9_.-]/_}
    coord_prefix="${s_p_chr1}_${p_start1}_${p_end1}_vs_${s_p_chr2}_${p_start2}_${p_end2}"
    unique_suffix=$(date +%s%N)_$$ 
    
    job_identifier="${coord_prefix}_${unique_suffix}" 
    worker_stderr_file="${worker_stderr_log_dir}/err_${job_identifier}.log"

    exec 2> "$worker_stderr_file"

    log "Worker started. Sanitized Inputs -> Seq1: $p_chr1:$p_start1-$p_end1 (InputStrand: $p_original_strand1). Seq2: $p_chr2:$p_start2-$p_end2 (InputStrand: $p_original_strand2)."
    log "Stderr for this job will be in: $worker_stderr_file"
    
    local tmp_seq1_fa tmp_seq2_fa tmp_alignment_paf final_output_part_file
    tmp_seq1_fa=$(mktemp "${worker_data_parts_dir}/minimap_tmp_${job_identifier}_s1_XXXXXX.fa")
    tmp_seq2_fa=$(mktemp "${worker_data_parts_dir}/minimap_tmp_${job_identifier}_s2_XXXXXX.fa")
    tmp_alignment_paf=$(mktemp "${worker_data_parts_dir}/minimap_tmp_${job_identifier}_paf_XXXXXX.paf")
    final_output_part_file="${worker_data_parts_dir}/part_${job_identifier}.txt"


    local cigar_out="N/A"
    local identity_rate_out="N/A"
    local final_s1_strand="${p_original_strand1}" 
    local final_s2_strand="${p_original_strand2}" 

    # --- Sequence 1 Extraction ---
    local s1_sam_start=$((p_start1 + 0)); s1_sam_start=$((s1_sam_start + 1)) 
    local s1_sam_end=$((p_end1 + 0))     
    local s1_region="${p_chr1}:${s1_sam_start}-${s1_sam_end}"
    local seq1_content
    
    log "Extracting seq1 region with command 'samtools faidx \"$fasta_file\" \"$s1_region\"'"
    seq1_content=$(samtools faidx "$fasta_file" "$s1_region" || true) 

    # Robust check for valid sequence content
    # True if: seq1_content is empty OR (it does not contain any newline OR the part after the first newline is empty)
    if [[ -z "$seq1_content" ]] || { [[ "$seq1_content" != *$'\n'* ]] || [[ -z "${seq1_content#*$'\n'}" ]]; }; then
        cigar_out="N/A_SEQ1_SAMTOOLS_FAILED_OR_EMPTY_NO_SEQ_LINES"
        log "SEQ1 FAILED or empty for region: '$s1_region'. Check .fai for '$p_chr1' and bounds."
        log "Raw seq1_content was: [$seq1_content]" # Log raw content for debugging this specific condition
        if ! grep -q -w "^${p_chr1}" "${fasta_file}.fai"; then log "Chromosome '$p_chr1}' NOT FOUND in .fai index."; fi
    else
        echo "$seq1_content" > "${tmp_seq1_fa}"
        log "SEQ1 extracted for '$s1_region' to $tmp_seq1_fa"

        # --- Sequence 2 Extraction ---
        local s2_sam_start=$((p_start2 + 0)); s2_sam_start=$((s2_sam_start + 1)) 
        local s2_sam_end=$((p_end2 + 0))     
        local s2_region="${p_chr2}:${s2_sam_start}-${s2_sam_end}"
        local seq2_content

        log "Extracting seq2 region with command 'samtools faidx \"$fasta_file\" \"$s2_region\"'"
        seq2_content=$(samtools faidx "$fasta_file" "$s2_region" || true) 

        if [[ -z "$seq2_content" ]] || { [[ "$seq2_content" != *$'\n'* ]] || [[ -z "${seq2_content#*$'\n'}" ]]; }; then
            cigar_out="N/A_SEQ2_SAMTOOLS_FAILED_OR_EMPTY_NO_SEQ_LINES"
            log "SEQ2 FAILED or empty for region: '$s2_region'. Check .fai for '$p_chr2' and bounds."
            log "Raw seq2_content was: [$seq2_content]"
            if ! grep -q -w "^${p_chr2}" "${fasta_file}.fai"; then log "Chromosome '$p_chr2}' NOT FOUND in .fai index."; fi
        else
            echo "$seq2_content" > "${tmp_seq2_fa}"
            log "SEQ2 extracted for '$s2_region' to $tmp_seq2_fa"

            log "Aligning '$tmp_seq1_fa' (target) and '$tmp_seq2_fa' (query) with command 'minimap2 -x asm20 -c --eqx ...'"
            minimap2 -x asm20 -c --eqx "${tmp_seq1_fa}" "${tmp_seq2_fa}" > "${tmp_alignment_paf}" || true

            if [[ -s "${tmp_alignment_paf}" ]]; then
                local paf_line_content extracted_cigar paf_strand
                paf_line_content=$(head -n 1 "${tmp_alignment_paf}")
                extracted_cigar=$(echo "${paf_line_content}" | awk '{for(i=1;i<=NF;i++) if($i ~ /^cg:Z:/){sub("cg:Z:","",$i); print $i; exit}}')
                
                if [[ -n "$extracted_cigar" ]]; then
                    cigar_out="$extracted_cigar"
                    paf_strand=$(echo "${paf_line_content}" | awk '{print $5}') 

                    if [[ "$paf_strand" == "+" ]] || [[ "$paf_strand" == "-" ]]; then
                        if [[ "${p_original_strand1}" == "+" ]]; then
                            final_s2_strand="$paf_strand"
                        elif [[ "${p_original_strand1}" == "-" ]]; then
                            if [[ "$paf_strand" == "+" ]]; then final_s2_strand="-"; else final_s2_strand="+"; fi
                        fi
                        log "Strand adjustment: Input s1_strand=$p_original_strand1, Input s2_strand=$p_original_strand2, PAF_relative_strand=$paf_strand -> Final s1_strand=$final_s1_strand, Final s2_strand=$final_s2_strand"
                    else
                        log "Warning: Unexpected PAF strand value '$paf_strand'. Original s2 strand '$p_original_strand2' will be kept for final_s2_strand."
                        final_s2_strand="${p_original_strand2}" 
                    fi

                    identity_rate_out=$(echo "$cigar_out" | awk '
                    {
                        matched_bases = 0; total_aligned_len = 0; current_cigar = $0;
                        while(length(current_cigar) > 0){
                            match(current_cigar, /^[0-9]+/); len_str = substr(current_cigar, RSTART, RLENGTH); current_cigar = substr(current_cigar, RSTART + RLENGTH);
                            match(current_cigar, /^[=XIDNSHP]/); op_char = substr(current_cigar, RSTART, RLENGTH); current_cigar = substr(current_cigar, RSTART + RLENGTH);
                            len_val = len_str + 0;
                            if (op_char == "=") { matched_bases += len_val; total_aligned_len += len_val; }
                            else if (op_char == "X" || op_char == "I" || op_char == "D") { total_aligned_len += len_val; }
                        }
                        if (total_aligned_len > 0) { printf "%.4f", matched_bases / total_aligned_len; } else { print "N/A"; }
                    }')
                    log "CIGAR=$cigar_out, Identity=$identity_rate_out"
                else
                    log "No CIGAR string found in PAF line: [$paf_line_content]" # Log content if no CIGAR
                    cigar_out="N/A_NO_CIGAR_IN_PAF"
                fi
            else
                log "Alignment PAF file $tmp_alignment_paf is empty. No alignment found by minimap2."
                cigar_out="N/A_NO_ALIGNMENT"
            fi 
        fi 
    fi 
    
    echo -e "${p_chr1}\t${p_start1}\t${p_end1}\t${p_chr2}\t${p_start2}\t${p_end2}\t${final_s1_strand}\t${final_s2_strand}\t${cigar_out}\t${identity_rate_out}" > "$final_output_part_file"
    
    echo "$final_output_part_file" 
    
    rm -f "${tmp_seq1_fa}" "${tmp_seq2_fa}" "${tmp_alignment_paf}"
    log "Worker finished. Data part: $final_output_part_file. Stderr log: $worker_stderr_file"
}


# Main function to drive parallel processing
process_bed_homology_parallel() {
    local input_bed_file="$1"
    local fasta_file="$2"
    local final_output_file="$3"
    local num_jobs="${4:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}"

    if [[ -z "$input_bed_file" || -z "$fasta_file" || -z "$final_output_file" ]]; then
        log "Usage: $0 process_bed_homology_parallel <input_bed_file> <reference_fasta> <final_output_file> [num_parallel_jobs]"
        return 1
    fi
    if [[ ! -f "$input_bed_file" ]]; then log "Error: BED file '$input_bed_file' not found."; return 1; fi
    if [[ ! -f "$fasta_file" ]]; then log "Error: FASTA file '$fasta_file' not found."; return 1; fi
    if [[ ! -f "${fasta_file}.fai" ]]; then log "Error: FASTA index '${fasta_file}.fai' not found. Please index with 'samtools faidx $fasta_file'."; return 1; fi
    if ! command -v samtools &> /dev/null; then log "Error: samtools not found."; return 1; fi
    if ! command -v minimap2 &> /dev/null; then log "Error: minimap2 not found."; return 1; fi
    if ! command -v parallel &> /dev/null; then log "Error: GNU parallel not found."; return 1; fi

    local joblog_file="${final_output_file}.joblog"
    local meta_file_path="${final_output_file}.meta_parts_list.txt"
    local worker_stderr_dir="${final_output_file}_worker_stderr_logs"
    local worker_data_parts_dir="${final_output_file}_worker_data_parts"


    log "Preparing worker stderr log directory: $worker_stderr_dir"
    if [[ -d "$worker_stderr_dir" ]]; then rm -rf "$worker_stderr_dir"; fi
    mkdir -p "$worker_stderr_dir" 
    log "Worker stderr log directory created."

    log "Preparing worker data parts directory: $worker_data_parts_dir"
    if [[ -d "$worker_data_parts_dir" ]]; then rm -rf "$worker_data_parts_dir"; fi
    mkdir -p "$worker_data_parts_dir"
    log "Worker data parts directory created."


    touch "$meta_file_path" 
    local header_line=""
    local header_detected=false
    local first_line
    first_line=$(head -n 1 "$input_bed_file" || true) 

    if [[ -n "$first_line" ]] && \
       ( [[ "$first_line" == \#* ]] || \
         [[ "$first_line" =~ ^track[[:space:]] ]] || \
         [[ "$first_line" =~ ^browser[[:space:]] ]] || \
         ( [[ "$first_line" =~ chr ]] && [[ "$first_line" =~ start ]] && [[ "$first_line" =~ end ]] && \
           ! [[ "$first_line" =~ [0-9]+$'\t'[+-]$'\t'.*$'\t'[+-]$ ]] && \
           ! [[ "$first_line" =~ [0-9]+$'\t'[+-]$ ]] ) ); then
        header_line="$first_line"
        header_detected=true
        echo "$header_line" > "$final_output_file"
        log "Detected and wrote header to $final_output_file: $header_line"
    else
        >"$final_output_file" 
        log "No header detected or input file empty. Final output file initialized."
    fi

    log "Starting parallel processing. Num jobs: $num_jobs. Parallel joblog: $joblog_file."
    log "Metafile for data parts: $meta_file_path. Worker stderr logs in: $worker_stderr_dir/. Worker data parts in: $worker_data_parts_dir/"

    >"$meta_file_path" 

    local input_stream_cmd
    if $header_detected; then
        if [[ $(wc -l < "$input_bed_file" | awk '{print $1}') -gt 1 ]]; then
            input_stream_cmd="tail -n +2 '$input_bed_file'"
        else
            input_stream_cmd=":" 
        fi
    elif [[ -s "$input_bed_file" ]]; then
        input_stream_cmd="cat '$input_bed_file'"
    else
        input_stream_cmd=":" 
    fi

    if [[ "$input_stream_cmd" != ":" ]]; then
        eval "$input_stream_cmd" | parallel --no-run-if-empty -j "$num_jobs" --colsep '\t' --joblog "$joblog_file" \
            bash "$0" _process_single_bed_line {1} {2} {3} {7} {4} {5} {6} {8} "$fasta_file" "$worker_stderr_dir" "$worker_data_parts_dir" >> "$meta_file_path"
    else
        log "No data lines to process after header handling or due to empty input file."
    fi
    
    log "Parallel processing finished. Concatenating results into $final_output_file..."

    if [[ -s "$meta_file_path" ]]; then
        while IFS= read -r part_file_path || [[ -n "$part_file_path" ]]; do
            if [[ -f "$part_file_path" ]]; then
                cat "$part_file_path" >> "$final_output_file"
            else
                log "Warning: Data part file '$part_file_path' listed in meta file not found."
            fi
        done < "$meta_file_path"
        log "Concatenation complete."
    else
        log "Meta file '$meta_file_path' is empty. No data output parts were generated or listed."
    fi

    local expected_data_lines=0
    if $header_detected; then
        if [[ $(wc -l < "$input_bed_file" | awk '{print $1}') -gt 1 ]]; then
             expected_data_lines=$(tail -n +2 "$input_bed_file" | wc -l | awk '{print $1}')
        fi
    elif [[ -s "$input_bed_file" ]]; then
        expected_data_lines=$(wc -l < "$input_bed_file" | awk '{print $1}')
    fi
    
    local processed_lines_in_meta=0
    if [[ -f "$meta_file_path" ]]; then
        processed_lines_in_meta=$(wc -l < "$meta_file_path" | awk '{print $1}')
    fi
    
    local lines_in_final_output_data=0
    if [[ -f "$final_output_file" ]]; then
        if $header_detected && [[ $(wc -l < "$final_output_file" | awk '{print $1}') -gt 0 ]]; then
            lines_in_final_output_data=$(tail -n +2 "$final_output_file" | wc -l | awk '{print $1}')
        elif [[ -s "$final_output_file" ]]; then
             lines_in_final_output_data=$(wc -l < "$final_output_file" | awk '{print $1}')
        fi
    fi

    log "Validation: Expected data lines from input: $expected_data_lines"
    log "Validation: Generated data output part files (from meta file): $processed_lines_in_meta"
    log "Validation: Data lines in final output file: $lines_in_final_output_data"

    local non_empty_worker_errors_found=false
    if [[ -d "$worker_stderr_dir" ]] && compgen -G "${worker_stderr_dir}/err_*.log" > /dev/null; then 
        for err_log in "${worker_stderr_dir}"/err_*.log; do
            if [[ -s "$err_log" ]]; then 
                non_empty_worker_errors_found=true
                # No longer log each one here, just set the flag and break
                break 
            fi
        done
    fi

    if $non_empty_worker_errors_found; then
        log "One or more non-empty worker stderr logs found in $worker_stderr_dir. Please review them for details on specific job failures."
    else
        log "No non-empty worker stderr logs found in $worker_stderr_dir."
    fi

    if [[ "$expected_data_lines" -eq "$processed_lines_in_meta" && "$processed_lines_in_meta" -eq "$lines_in_final_output_data" ]]; then
        local input_was_effectively_empty=false
        if [[ "$expected_data_lines" -eq 0 ]]; then
            if ! $header_detected && ! [[ -s "$input_bed_file" ]]; then 
                input_was_effectively_empty=true
            elif $header_detected && [[ $(wc -l < "$input_bed_file" | awk '{print $1}') -le 1 ]]; then 
                input_was_effectively_empty=true
            fi
        fi

        if $input_was_effectively_empty; then
            log "Input file was effectively empty (no data lines). Output is consistent."
        else
            log "Success: All $expected_data_lines data lines appear to be processed and written to $final_output_file."
        fi
        
        log "Cleaning up temporary files and directories..."
        if [[ -f "$meta_file_path" ]]; then
            rm -f "$meta_file_path"
            log "Meta file removed."
        fi
        
        if [[ -d "$worker_data_parts_dir" ]]; then
             log "Removing worker data parts directory: $worker_data_parts_dir"
             rm -rf "$worker_data_parts_dir"
        fi
        
        if ! $non_empty_worker_errors_found && [[ -d "$worker_stderr_dir" ]]; then
            log "No worker errors detected, removing worker stderr log directory: $worker_stderr_dir"
            rm -rf "$worker_stderr_dir"
        else
            log "Worker stderr logs (or directory) kept for review in: $worker_stderr_dir"
        fi
    else
        log "Warning: Line count mismatch or processing issues. Please check joblog and meta file."
        log "Individual data output part files are in '$worker_data_parts_dir'."
        log "Individual worker stderr logs are in '$worker_stderr_dir'."
    fi
    log "Script finished."
}


# Script execution logic
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if declare -f "$1" &> /dev/null; then
        "$@" 
    else
        if [[ "$#" -eq 0 ]]; then
             log "Main usage: $0 process_bed_homology_parallel <input_bed_file> <reference_fasta> <final_output_file> [num_jobs]"
             log " (This script also defines _process_single_bed_line for internal parallel use)."
             exit 1
        elif [[ "$1" == "process_bed_homology_parallel" && "$#" -ge 3 ]]; then
             process_bed_homology_parallel "${@:2}"
        else
            log "Error: Function '$1' not found or invalid arguments for direct execution."
            log "Main usage: $0 process_bed_homology_parallel <input_bed_file> <reference_fasta> <final_output_file> [num_jobs]"
            exit 1
        fi
    fi
fi
