use std::collections::HashSet;
use rust_htslib::bam::{self, Read, Record};
use crate::structs::{ReadPair, SortedVecIntervals, Interval, ReadPairMap, AlleleDepthMap, PositionAlleleDepth};
use statrs::statistics::OrderStatistics; // For median calculation
use polars::prelude::*;
use std::process::Command;
use rustc_hash::FxHashMap; // Faster HashMap for integer keys
use ahash::AHashMap; // Faster HashMap for string keys
use std::path::Path;
use tempfile::NamedTempFile;


#[derive(Debug, PartialEq)]
enum ReadPairStatus {
    Incomplete,
    Complete,
}

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
    // Only evaluate primary alignments for noise - paired-end specific checks
    if !read.is_secondary() && !read.is_supplementary() {
        if read.mapq() < mapq_filter ||
           read.is_quality_check_failed() ||
           read.is_unmapped() ||
           read.reference_end() < read.reference_start() + 75 ||
           read.tid() != read.mtid() ||
           !read.is_proper_pair() ||
           read.seq_len() == 0 {
            return true;
        }

        if filter_noisy && !read.qual().is_empty() {
            let qual_vec: Vec<f64> = read.qual().iter().map(|&q| q as f64).collect();
            // Use the mature statrs crate for median calculation
            let median_qual = qual_vec.median() as u8;
            let low_qual_count = read.qual().iter().filter(|&&q| q < basequal_median_filter).count();
            
            let soft_clip_bases: u32 = read.cigar()
                .iter()
                .filter(|c| c.char() == 'S')
                .map(|c| c.len())
                .sum();

            if median_qual <= basequal_median_filter ||
               low_qual_count >= 50 ||
               soft_clip_bases >= 20 {
                return true;
            }
        }
    }

    false
}

/// Optional: Run samtools collate to group reads by qname
/// This produces a temporary file that is automatically deleted when the function returns
fn collate_bam_file(bam_file_path: &str, threads: &u8 = &2) -> Result<Option<NamedTempFile>, Box<dyn std::error::Error>> {
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
    threads: u8 = 4) -> Result<ReadPairMap, Box<dyn std::error::Error>> {

    // Optionally collate the BAM file first, make sure the temp file end with .bam suffix
    let (bam_path, _temp_file) = if use_collate {
        // ? operator: propagates any error from collate_bam_file function
        match collate_bam_file(bam_file_path, &threads)? {
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
    
    let mut result = ReadPairMap::new();
    let mut qname_idx_counter = 0usize;
    
    // Pre-allocate chromosome interval trees
    let chrom_count = header.target_count() as usize;
    result.interval_trees.reserve(chrom_count);
    
    for tid in 0..header.target_count() {
        if let Some(ref_name) = header.tid2name(tid) {
            let chrom = String::from_utf8_lossy(ref_name).to_string();
            result.interval_trees.insert(chrom, SortedVecIntervals::new());
        }
    }
    
    // Buffer for collecting reads with the same qname
    let mut current_qname: Option<String> = None;
    let mut current_reads: Vec<Record> = Vec::with_capacity(2);
    
    // Process reads grouped by qname
    for read in bam.records() {
        // Skip secondary, supplementary, and duplicate alignments
        if should_skip_alignment(&read) {
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
    
    // The temp file will be automatically deleted when _temp_file goes out of scope
    Ok(result)
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
        if let Some(tid) = read.tid() {
            if let Some(ref_name) = header.tid2name(tid as u32) {
                let chrom = String::from_utf8_lossy(ref_name).to_string();
                // ? operator: propagates any error from add_read_interval
                add_read_interval(&mut result.interval_trees, &chrom, read, qname_idx)?;
            }
        }
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
    let ad_file = format!("{}.ad", bam_file);
    
    // Use bcftools mpileup (keeping the efficient C implementation)
    let output = Command::new("bcftools")
        .args([
            "mpileup", "-Ou", 
            "--fasta-ref", reference_genome,
            "-a", "FORMAT/AD",
            "--indels-2.0",
            "-q", &mapq_filter.to_string(),
            "-Q", &base_qual_filter.to_string(),
            bam_file
        ])
        .output()?; // ? operator: propagates any error from command execution
    
    if !output.status.success() {
        return Err(format!("bcftools mpileup failed: {}", String::from_utf8_lossy(&output.stderr)).into());
    }
    
    // Pipe to bcftools query
    let query_output = Command::new("bcftools")
        .args([
            "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n", "-"
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()?; // ? operator: propagates any error from process spawning
    
    // Write the mpileup output to the query command
    let mut child = query_output;
    if let Some(mut stdin) = child.stdin.take() {
        use std::io::Write;
        // ? operator: propagates any error from writing to stdin
        stdin.write_all(&output.stdout)?;
    }
    
    // ? operator: propagates any error from waiting for child process
    let query_output = child.wait_with_output()?;
    
    if !query_output.status.success() {
        return Err(format!("bcftools query failed: {}", String::from_utf8_lossy(&query_output.stderr)).into());
    }
    
    // Write to temporary file for Polars to read
    // ? operator: propagates any error from file writing
    std::fs::write(&ad_file, &query_output.stdout)?;
    
    // Use Polars for efficient data processing (much faster than pandas)
    let mut df = LazyFrame::scan_csv(
        &ad_file,
        ScanArgsCSV::default()
            .with_separator(b'\t')
            .with_has_header(false)
            .with_column_names(Some(vec![
                "chrom".to_string(),
                "pos".to_string(), 
                "ref".to_string(),
                "alt".to_string(),
                "ad".to_string()
            ]))
    )?  // ? operator: propagates any error from CSV scanning
    .filter(col("alt").is_not_null())
    .filter(col("alt").neq(lit("<*>")))
    .collect()?; // ? operator: propagates any error from DataFrame collection
    
    // Convert to efficient Rust data structure
    let mut allele_depth_map = AlleleDepthMap::new();
    
    for row in df.iter() {
        // ? operator: propagates any error from row access and type extraction
        let chrom: &str = row.get(0)?.try_extract()?;
        let pos: u32 = (row.get(1)?.try_extract::<i64>()? - 1) as u32; // Convert to 0-based u32
        let ref_allele: &str = row.get(2)?.try_extract()?;
        let alt_alleles: &str = row.get(3)?.try_extract()?;
        let ad_str: &str = row.get(4)?.try_extract()?;
        
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
        
        // Parse ALT alleles  
        let alt_list: Vec<&str> = alt_alleles
            .trim_end_matches(",<*>")
            .split(',')
            .collect();
            
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
    }
    
    // Clean up temporary file
    std::fs::remove_file(&ad_file).ok();
    
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

