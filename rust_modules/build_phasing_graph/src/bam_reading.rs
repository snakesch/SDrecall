use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::ext::BamRecordExtensions;
use crate::structs::{ReadPair, SortedVecIntervals, ReadPairMap, AlleleDepthMap};

use std::process::{Command, Stdio};
use std::io::{BufRead, BufReader};
use ahash::AHashMap; // Faster HashMap for string keys
use std::path::Path;
use std::thread;
use tempfile::NamedTempFile;
use log::{debug, info, warn, error};
use std::sync::{Arc, Mutex};
use std::collections::VecDeque;


fn should_skip_alignment(read: &Record) -> bool {
    // Skip secondary, supplementary, and duplicate alignments early
    // These are not quality issues but different types of alignments
    read.is_secondary() || read.is_supplementary() || read.is_duplicate()
}

fn is_read_noisy(
    read: &Record,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> bool {
    let qname = String::from_utf8_lossy(read.qname());
    
    // Only evaluate primary alignments for noise - secondary or supplementary alignments are skipped
    if read.is_secondary() || read.is_supplementary() {
        debug!("[is_read_noisy] skip_secondary_supplementary - {} skipped (secondary/supplementary alignment); not considered noisy", qname);
        return false;
    }
    
    // Common fast checks for both paired and unpaired
    if read.is_unmapped() {
        debug!("[is_read_noisy] unmapped_check - {} flagged noisy: unmapped read", qname);
        return true;
    }
    
    if read.is_quality_check_failed() {
        debug!("[is_read_noisy] qc_fail_check - {} flagged noisy: QC fail flag set", qname);
        return true;
    }
    
    if read.mapq() < mapq_filter {
        debug!("[is_read_noisy] mapq_check - {} flagged noisy: MAPQ {} < threshold {}", qname, read.mapq(), mapq_filter);
        return true;
    }
    
    if read.seq_len() == 0 {
        debug!("[is_read_noisy] seq_len_check - {} flagged noisy: missing query_sequence", qname);
        return true;
    }
    
    // Paired-end specific checks
    let aln_len = read.reference_end() - read.reference_start();
    if aln_len < 75 {
        debug!("[is_read_noisy] aln_len_check - {} flagged noisy: alignment span {} < 75", qname, aln_len);
        return true;
    }
    
    // Ensure mate on same reference for proper pairing in this pipeline
    if read.tid() != read.mtid() {
        debug!("[is_read_noisy] mate_tid_check - {} flagged noisy: different reference chromosomes for mate pair", qname);
        return true;
    }
    
    if !read.is_proper_pair() {
        debug!("[is_read_noisy] proper_pair_check - warning:{} is not a proper pair", qname);
    }
    
    // Base quality and soft-clip based checks (controlled by filter_noisy)
    if filter_noisy && !read.qual().is_empty() {
        let mut qual_vec: Vec<u8> = read.qual().to_vec();
        // Use order-stat crate for O(n) median calculation (most efficient)
        let median_qual = if !qual_vec.is_empty() {
            let len = qual_vec.len();
            *order_stat::kth(&mut qual_vec, len / 2)
        } else {
            debug!("[is_read_noisy] median_qual_check - {} flagged noisy: empty qual vector", qname);
            return true;
        };
        
        if median_qual <= basequal_median_filter {
            debug!("[is_read_noisy] median_qual_check - {} flagged noisy: median baseQ {} <= threshold {}", qname, median_qual, basequal_median_filter);
            return true;
        }
        
        let low_qual_count = read.qual().iter().filter(|&&q| q < basequal_median_filter).count();
        if low_qual_count >= 50 {
            debug!("[is_read_noisy] low_qual_count_check - {} flagged noisy: #bases with Q<{} is {} >= 50", qname, basequal_median_filter, low_qual_count);
            return true;
        }
        
        let soft_clip_bases: u32 = read.cigar()
            .iter()
            .filter(|c| c.char() == 'S')
            .map(|c| c.len())
            .sum();
        
        if soft_clip_bases >= 20 {
            debug!("[is_read_noisy] soft_clip_check - {} flagged noisy: total soft-clip length {} >= 20\n", qname, soft_clip_bases);
            return true;
        }
    }
    
    // If none of the noisy conditions triggered
    false
}

/// Optional: Run samtools collate to group reads by qname
/// This produces a temporary file that is automatically deleted when the function returns
fn collate_bam_file(bam_file_path: &str, threads: u8) -> Result<Option<NamedTempFile>, Box<dyn std::error::Error>> {
    // Check if samtools is available
    let samtools_check = Command::new("samtools")
        .arg("--version")
        .output();
    
    if samtools_check.is_err() {
        // samtools not available, return None to use original file
        return Ok(None);
    }
    
    // Create a temporary file for collated output, make sure the temp file end with .bam suffix
    // ? operator: propagates any error from NamedTempFile creation up to the caller
    let temp_file = NamedTempFile::with_suffix_in(".bam", Path::new("."))?;
    let temp_path = temp_file.path().to_str()
        // ok_or: converts Option<&str> to Result<&str, &str> - if None, returns Err with the message
        .ok_or("Failed to convert temp path to string")?; // ? operator: propagates the error if conversion failed
    
    // Run samtools collate with fast mode for efficiency
    let output = Command::new("samtools")
        .args([
            "collate",
            "-f",           // Fast mode - primary alignments only
            "-@", &threads.to_string(),      // Use threads
            bam_file_path,
            "-o", temp_path
        ])
        .output()?; // ? operator: propagates any error from command execution
    
    if !output.status.success() {
        return Err(format!("samtools collate failed: {}", 
            String::from_utf8_lossy(&output.stderr)).into());
    }
    
    Ok(Some(temp_file))
}

/// Process BAM file by iterating through qname-grouped reads
/// This is more efficient than the original single-read iteration
pub fn migrate_bam_to_sorted_intervals_grouped(
    bam_file_path: &str,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
    use_collate: bool,
    threads: u8) -> Result<(ReadPairMap, rust_htslib::bam::HeaderView), Box<dyn std::error::Error>> {

    info!("[migrate_bam_to_sorted_intervals_grouped] Starting BAM processing: file={}, mapq_filter={}, basequal_median_filter={}, filter_noisy={}, use_collate={}, threads={}", 
         bam_file_path, mapq_filter, basequal_median_filter, filter_noisy, use_collate, threads);

    // Optionally collate the BAM file first, make sure the temp file end with .bam suffix
    let (bam_path, _temp_file) = if use_collate {
        // ? operator: propagates any error from collate_bam_file function
        match collate_bam_file(bam_file_path, threads)? {
            Some(temp_file) => {
                let path = temp_file.path().to_str()
                    // ok_or: converts Option<&str> to Result<&str, &str>
                    .ok_or("Failed to convert temp path to string")? // ? operator: propagates error if path conversion fails
                    .to_string();
                (path, Some(temp_file))
            },
            None => (bam_file_path.to_string(), None)
        }
    } else {
        (bam_file_path.to_string(), None)
    };
    
    // ? operator: propagates any error from BAM file opening
    let mut bam = bam::Reader::from_path(&bam_path)?;
    let header = bam.header().clone();
    info!("[migrate_bam_to_sorted_intervals_grouped] BAM file opened successfully, {} chromosomes found", header.target_count());
    
    let mut result = ReadPairMap::new();
    let mut qname_idx_counter = 0usize;
    
    // Pre-allocate chromosome interval trees
    let chrom_count = header.target_count() as usize;
    result.interval_trees.reserve(chrom_count);
    
    for tid in 0..header.target_count() {
        let ref_name = header.tid2name(tid);
        let chrom = String::from_utf8_lossy(ref_name).to_string();
        result.interval_trees.insert(chrom, SortedVecIntervals::new());
    }
    info!("[migrate_bam_to_sorted_intervals_grouped] Interval trees initialized for {} chromosomes", chrom_count);
    
    // Buffer for collecting reads with the same qname
    let mut current_qname: Option<String> = None;
    let mut current_reads: Vec<Record> = Vec::with_capacity(2);
    
    let mut total_reads_processed = 0usize;
    let mut skipped_alignments = 0usize;
    info!("[migrate_bam_to_sorted_intervals_grouped] Starting to process BAM records");
    
    // Process reads grouped by qname
    for read_result in bam.records() {
        let read = read_result?;
        total_reads_processed += 1;
        
        // Skip secondary, supplementary, and duplicate alignments
        if should_skip_alignment(&read) {
            skipped_alignments += 1;
            continue;
        }
        
        let qname = String::from_utf8_lossy(read.qname()).to_string();
        
        // Check if we've moved to a new qname
        // as_ref(): converts &Option<String> to Option<&String> for comparison
        if current_qname.as_ref() != Some(&qname) {
            // Process the previous qname group if any
            if let Some(prev_qname) = current_qname.take() {
                process_qname_group(
                    &mut result,
                    &header,
                    prev_qname,
                    &mut current_reads,
                    &mut qname_idx_counter,
                    mapq_filter,
                    basequal_median_filter,
                    filter_noisy,
                )?; // ? operator: propagates any error from process_qname_group
            }
            
            current_qname = Some(qname.clone());
            current_reads.clear();
        }
        
        // Skip if this qname is already marked as noisy
        if result.noisy_qnames.contains_key(&qname) {
            current_reads.clear();
            continue;
        }
        
        current_reads.push(read);
    }
    
    // Process the last qname group
    if let Some(qname) = current_qname {
        process_qname_group(
            &mut result,
            &header,
            qname,
            &mut current_reads,
            &mut qname_idx_counter,
            mapq_filter,
            basequal_median_filter,
            filter_noisy,
        )?; // ? operator: propagates any error from process_qname_group
    }
    
    // Finalize all interval trees
    for interval_tree in result.interval_trees.values_mut() {
        // ? operator: propagates any error from finalize
        interval_tree.finalize()?;
    }
    
    info!("[migrate_bam_to_sorted_intervals_grouped] BAM processing complete: {} total reads processed, {} alignments skipped, {} read pairs retained, {} noisy qnames filtered",
         total_reads_processed, skipped_alignments, result.readpair_dict.len(), result.noisy_qnames.len());
    
    // The temp file will be automatically deleted when _temp_file goes out of scope
    Ok((result, header))
}

/// Process all reads for a single qname
fn process_qname_group(
    result: &mut ReadPairMap,
    header: &bam::HeaderView,
    qname: String,
    reads: &mut Vec<Record>,
    qname_idx_counter: &mut usize,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool) -> Result<(), Box<dyn std::error::Error>> {
    // Check if any read in the group is noisy
    let is_noisy = reads.iter().any(|read| 
        is_read_noisy(read, mapq_filter, basequal_median_filter, filter_noisy)
    );
    
    if is_noisy {
        debug!("[process_qname_group] This qname {} is noisy. Skip it.\n", qname);
        result.noisy_qnames.insert(qname.clone(), ());
        return Ok(());
    }
    
    // Filter to get exactly 2 primary reads (read1 and read2)
    let mut read1_opt: Option<Record> = None;
    let mut read2_opt: Option<Record> = None;
    
    for read in reads.iter() {
        if read.is_first_in_template() {
            read1_opt = Some(read.clone());
        } else if read.is_last_in_template() {
            read2_opt = Some(read.clone());
        }
    }
    
    // We need both reads for a complete pair
    let (read1, read2) = match (read1_opt, read2_opt) {
        (Some(r1), Some(r2)) => (r1, r2),
        _ => return Ok(()), // Skip incomplete pairs
    };
    
    // Assign qname_idx and create ReadPair
    // * operator: dereferences the mutable reference to get the actual usize value
    let qname_idx = *qname_idx_counter;
    result.qname_to_idx.insert(qname.clone(), qname_idx);
    result.idx_to_qname.insert(qname_idx, qname.clone());
    // * operator: dereferences the mutable reference to modify the actual usize value
    *qname_idx_counter += 1;
    
    // Add intervals for both reads
    for read in [&read1, &read2] {
        let tid = read.tid();
        let ref_name = header.tid2name(tid as u32);
        let chrom = String::from_utf8_lossy(ref_name).to_string();
        // ? operator: propagates any error from add_read_interval
        add_read_interval(&mut result.interval_trees, &chrom, read, qname_idx)?;
    }
    
    // Create and store the complete ReadPair
    let readpair = ReadPair::new_complete(read1, read2, qname, qname_idx);
    result.readpair_dict.insert(qname_idx, readpair);
    
    Ok(())
}


/// Generalized function to add interval for any read record
fn add_read_interval(
    interval_trees: &mut AHashMap<String, SortedVecIntervals>,
    chrom: &str,
    read: &Record,
    qname_idx: usize,
) -> Result<(), String> {
    let start = read.reference_start();
    let end = read.reference_end();
    
    if let Some(interval_tree) = interval_trees.get_mut(chrom) {
        interval_tree.add_interval(start, end, qname_idx)
    } else {
        Err(format!("Chromosome {} not found in interval trees", chrom))
    }
}

/// Equivalent to Python's stat_ad_to_dict function
/// Builds allele depth information using bcftools mpileup and processes with Polars
/// Returns a more efficient Rust data structure instead of nested HashMaps
pub fn build_allele_depth_map(
    bam_file: &str,
    reference_genome: &str,
    mapq_filter: u8,
    base_qual_filter: u8,
) -> Result<AlleleDepthMap, Box<dyn std::error::Error>> {
    
    info!("[build_allele_depth_map] Starting with BAM: {}, MAPQ>={}, BaseQ>={}", bam_file, mapq_filter, base_qual_filter);
    
    // Stream mpileup stdout directly into query stdin, and write query stdout to the .ad file
    let mut mpileup_child = Command::new("bcftools")
        .args([
            "mpileup", "-Ou",
            "--fasta-ref", reference_genome,
            "-a", "FORMAT/AD",
            "--indels-2.0",
            "-A",
            "-q", &mapq_filter.to_string(),
            "-Q", &base_qual_filter.to_string(),
            bam_file,
        ])
        .stdout(Stdio::piped())
		.stderr(Stdio::piped())
        .spawn()?;

    // drain mpileup stderr concurrently into a bounded buffer (no logging here to avoid GIL contention)
    let mpileup_stderr_buf: Arc<Mutex<VecDeque<String>>> = Arc::new(Mutex::new(VecDeque::with_capacity(500)));
    let mpileup_stderr_buf_reader = Arc::clone(&mpileup_stderr_buf);
	let mpileup_stderr_jh = if let Some(stderr) = mpileup_child.stderr.take() {
		Some(thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line_res in reader.lines() {
                if let Ok(line) = line_res {
                    let mut buf = mpileup_stderr_buf_reader.lock().unwrap();
                    if buf.len() == buf.capacity() { buf.pop_front(); }
                    buf.push_back(line);
                }
			}
		}))
	} else { None };

    let mpileup_stdout = mpileup_child.stdout.take().ok_or("Failed to capture mpileup stdout")?;

    let mut query_child = Command::new("bcftools")
        .args(["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\\n", "-"])
        .stdin(Stdio::from(mpileup_stdout))
        .stdout(Stdio::piped())
		.stderr(Stdio::piped())
        .spawn()?;

    // drain query stderr concurrently into a bounded buffer (no logging here)
    let query_stderr_buf: Arc<Mutex<VecDeque<String>>> = Arc::new(Mutex::new(VecDeque::with_capacity(500)));
    let query_stderr_buf_reader = Arc::clone(&query_stderr_buf);
	let query_stderr_jh = if let Some(stderr) = query_child.stderr.take() {
		Some(thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line_res in reader.lines() {
                if let Ok(line) = line_res {
                    let mut buf = query_stderr_buf_reader.lock().unwrap();
                    if buf.len() == buf.capacity() { buf.pop_front(); }
                    buf.push_back(line);
                }
			}
		}))
	} else { None };

    // Stream-parse bcftools query output directly into allele_depth_map
    let query_stdout = query_child.stdout.take().ok_or("Failed to capture bcftools query stdout")?;

    let reader = BufReader::new(query_stdout);
    let mut allele_depth_map = AlleleDepthMap::new();
    let mut line_count = 0usize;
    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() { continue; }
        line_count += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        
        let chrom = fields[0];
        let pos: u32 = match fields[1].parse::<u32>() {
            Ok(p) => p - 1, // Convert to 0-based
            Err(_) => continue,
        };
        let ref_allele = fields[2];
        let alt_alleles = fields[3];
        let ad_str = fields[4];
        
        // Even if ALT is empty, still record REF depth/DP (positions with no ALT)
        
        // Parse AD values
        let ad_values: Vec<u16> = ad_str
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();
            
        if ad_values.is_empty() {
            continue;
        }
        
        let ref_depth = ad_values[0];
        let total_depth: u16 = ad_values.iter().sum();
        
        if total_depth == 0 {
            continue;
        }
        
        // Parse ALT alleles  (filter placeholders and empties)
        let mut alt_list: Vec<&str> = alt_alleles
            .trim_end_matches(",<*>")
            .split(',')
            .collect();

        // Remove placeholders and empty entries
        alt_list.retain(|&alt| !alt.is_empty() && alt != "<*>");

        // Only report positions with at least one ALT allele after filtering
        if alt_list.is_empty() {
            continue;
        }
            
        // Create position entry - just a simple array
        let mut position_data = AlleleDepthMap::new_position_data(total_depth);
        
        // Add reference depth
        AlleleDepthMap::set_allele_depth(
            &mut position_data, 
            base_to_index(ref_allele.chars().next().unwrap_or('N')), 
            ref_depth
        );
        
        // Add alt allele depths
        for (i, &alt_allele) in alt_list.iter().enumerate() {
            if i + 1 < ad_values.len() && !alt_allele.is_empty() {
                let alt_depth = ad_values[i + 1];
                if alt_depth > 0 {
                    // Only handle SNVs for now (same as Python version)
                    if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        let alt_char = alt_allele.chars().next().unwrap();
                        AlleleDepthMap::set_allele_depth(
                            &mut position_data,
                            base_to_index(alt_char), 
                            alt_depth
                        );
                    }
                }
            }
        }
        
        allele_depth_map.insert(chrom, pos, position_data);
        debug!("[build_allele_depth_map] Inserted position data for {} at position {} with ref {} and alt {} where the ADs are {:?}", chrom, pos, ref_allele, alt_alleles, position_data);
    }
    info!("[build_allele_depth_map] Processed {} lines, created map with {} chromosomes", line_count, allele_depth_map.chromosome_count());

    // Ensure both processes exited successfully
    let query_status = query_child.wait()?;
    let mpileup_status = mpileup_child.wait()?;
    if let Some(jh) = mpileup_stderr_jh { let _ = jh.join(); }
	if let Some(jh) = query_stderr_jh { let _ = jh.join(); }

    // Flush buffered stderr lines to logger now (on main thread)
    if let Ok(buf) = mpileup_stderr_buf.lock() {
        for line in buf.iter() {
            let lower = line.to_lowercase();
            if lower.contains("error") { error!("{}", line); }
            else if lower.contains("warn") { warn!("{}", line); }
            else { info!("{}", line); }
        }
    }
    if let Ok(buf) = query_stderr_buf.lock() {
        for line in buf.iter() {
            let lower = line.to_lowercase();
            if lower.contains("error") { error!("{}", line); }
            else if lower.contains("warn") { warn!("{}", line); }
            else { info!("{}", line); }
        }
    }

    if !mpileup_status.success() {
        return Err("bcftools mpileup failed".into());
    }
    if !query_status.success() {
        return Err("bcftools query failed".into());
    }

    Ok(allele_depth_map)
}

/// Convert DNA base to array index (more efficient than HashMap lookups)
#[inline]
fn base_to_index(base: char) -> usize {
    match base.to_ascii_uppercase() {
        'A' => 0,
        'T' => 1, 
        'C' => 2,
        'G' => 3,
        _ => 4,  // N or any other base
    }
}

