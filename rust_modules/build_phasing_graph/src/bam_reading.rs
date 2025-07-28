use std::collections::{HashMap, HashSet};
use rust_htslib::bam::{self, Read, Record};
use crate::structs::{ReadPair, SortedVecIntervals, Interval, ReadPairMap};
use statrs::statistics::OrderStatistics; // For median calculation


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

pub fn migrate_bam_to_sorted_intervals(
    bam_file_path: &str,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> Result<ReadPairMap, Box<dyn std::error::Error>> {
    let mut bam = bam::Reader::from_path(bam_file_path)?;  // Now mutable
    let header = bam.header().clone();
    
    let mut result = ReadPairMap::new();
    let mut qname_idx_counter = 0usize;
    let mut read_pair_status: HashMap<usize, ReadPairStatus> = HashMap::new();
    let mut temp_readpairs: HashMap<usize, ReadPair> = HashMap::new();
    let mut incomplete_qname_indices: HashSet<usize> = HashSet::new();

    // Pre-allocate HashMap with estimated capacity based on chromosome count
    let chrom_count = header.target_count() as usize;
    result.interval_trees.reserve(chrom_count);

    // Initialize empty SortedVecIntervals for each chromosome
    for tid in 0..header.target_count() {
        if let Some(ref_name) = header.tid2name(tid) {
            let chrom = String::from_utf8_lossy(ref_name).to_string();
            result.interval_trees.insert(chrom, SortedVecIntervals::new());
        }
    }

    // Single pass: process each read and add intervals directly
    for read_result in bam.records() {
        let read = read_result?;
        
        // Skip secondary, supplementary, and duplicate alignments early
        if should_skip_alignment(&read) {
            continue;
        }
        
        let qname = String::from_utf8_lossy(read.qname()).to_string();
        
        if result.noisy_qnames.contains(&qname) {
            continue;
        }

        // Only check for noise in primary alignments
        if is_read_noisy(&read, mapq_filter, basequal_median_filter, filter_noisy) {
            result.noisy_qnames.insert(qname.clone());
            // Efficiently remove all intervals for this qname_idx from all chromosomes
            if let Some(&qname_idx) = result.qname_to_idx.get(&qname) {
                read_pair_status.remove(&qname_idx);
                incomplete_qname_indices.remove(&qname_idx);
                temp_readpairs.remove(&qname_idx);
                
                // Remove from all chromosome interval trees
                for interval_tree in result.interval_trees.values_mut() {
                    interval_tree.remove_qname_idx(qname_idx)?;
                }
            }
            continue;
        }

        // Get or assign qname_idx
        let qname_idx = if let Some(&idx) = result.qname_to_idx.get(&qname) {
            idx
        } else {
            let idx = qname_idx_counter;
            result.qname_to_idx.insert(qname.clone(), idx);
            result.idx_to_qname.insert(idx, qname.clone());
            qname_idx_counter += 1;
            idx
        };

        // Add interval for this individual read - using match
        match read.tid() {
            Some(tid) => {
                match header.tid2name(tid as u32) {
                    Some(ref_name) => {
                        let chrom = String::from_utf8_lossy(ref_name).to_string();
                        add_read_interval(&mut result.interval_trees, &chrom, &read, qname_idx)?;
                    }
                    None => {} // Skip if chromosome name not found
                }
            }
            None => {} // Skip if read has no valid tid
        }

        // Handle read pairing for paired-end data
        match read_pair_status.get(&qname_idx) {
            Some(ReadPairStatus::Incomplete) => {
                // Complete the existing readpair
                if let Some(readpair) = temp_readpairs.get_mut(&qname_idx) {
                    readpair.complete_with_read2(read);
                    read_pair_status.insert(qname_idx, ReadPairStatus::Complete);
                    incomplete_qname_indices.remove(&qname_idx);
                }
            }
            Some(ReadPairStatus::Complete) => continue,
            None => {
                // First read for this qname_idx
                let incomplete_readpair = ReadPair::new_incomplete(read, qname, qname_idx);
                temp_readpairs.insert(qname_idx, incomplete_readpair);
                read_pair_status.insert(qname_idx, ReadPairStatus::Incomplete);
                incomplete_qname_indices.insert(&qname_idx);
            }
        }
    }

    // Handle incomplete pairs by fetching mates
    if !incomplete_qname_indices.is_empty() {
        let mut indexed_bam = bam::IndexedReader::from_path(bam_file_path)?;
        
        let incomplete_indices: Vec<usize> = incomplete_qname_indices.iter().copied().collect();
        for qname_idx in incomplete_indices {
            let Some(readpair) = temp_readpairs.get(&qname_idx) else { continue; };
            let query_read = &readpair.read1;
            
            let Some(mate_read) = fetch_mate_read(&mut indexed_bam, query_read, mapq_filter, basequal_median_filter, filter_noisy)? else { continue; };
            
            // Add mate interval BEFORE completing readpair - handle unmapped mates properly
            let interval_added = mate_read.tid()
                .and_then(|tid| header.tid2name(tid as u32))
                .map(|ref_name| {
                    let chrom = String::from_utf8_lossy(ref_name).to_string();
                    add_read_interval(&mut result.interval_trees, &chrom, &mate_read, qname_idx)
                })
                .transpose()?; // Convert Option<Result<T, E>> to Result<Option<T>, E>
            
            let Some(_) = interval_added else { 
                // Mate is unmapped, remove this qname_idx from incomplete and abandon
                incomplete_qname_indices.remove(&qname_idx);
                temp_readpairs.remove(&qname_idx);
                read_pair_status.remove(&qname_idx);
                continue; 
            };
            
            let Some(readpair_mut) = temp_readpairs.get_mut(&qname_idx) else { continue; };
            
            // Complete the readpair with the mate read
            readpair_mut.complete_with_read2(mate_read);
            read_pair_status.insert(qname_idx, ReadPairStatus::Complete);
            incomplete_qname_indices.remove(&qname_idx); 
        }
    }

    // Remove incomplete pairs efficiently
    for &qname_idx in &incomplete_qname_indices {
        temp_readpairs.remove(&qname_idx);
        for interval_tree in result.interval_trees.values_mut() {
            interval_tree.remove_qname_idx(qname_idx)?;
        }
    }

    // Move completed ReadPairs to final result
    for (qname_idx, status) in read_pair_status {
        if status == ReadPairStatus::Complete {
            if let Some(readpair) = temp_readpairs.remove(&qname_idx) {
                result.readpair_dict.insert(qname_idx, readpair);
            }
        }
    }

    // Finalize all interval trees
    for interval_tree in result.interval_trees.values_mut() {
        interval_tree.finalize()?;
    }

    Ok(result)
}

/// Generalized function to find the mate read for any given read record
fn fetch_mate_read(
    indexed_bam: &mut bam::IndexedReader,
    query_read: &Record,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> Result<Option<Record>, Box<dyn std::error::Error>> {
    let (mate_tid, mate_pos) = match (query_read.mtid(), query_read.mpos()) {
        (Some(tid), pos) if pos >= 0 => (tid, pos),
        _ => return Ok(None),
    };

    let fetch_start = (mate_pos - 5).max(0);
    let fetch_end = mate_pos + 10;
    
    indexed_bam.fetch((mate_tid as u32, fetch_start as u64, fetch_end as u64))?;
    
    let expected_qname = String::from_utf8_lossy(query_read.qname()).to_string();
    
    for mate_result in indexed_bam.records() {
        let mate_read = mate_result?;
        let mate_qname = String::from_utf8_lossy(mate_read.qname()).to_string();
        
        if mate_qname == expected_qname && 
           !should_skip_alignment(&mate_read) &&
           !is_read_noisy(&mate_read, mapq_filter, basequal_median_filter, filter_noisy) {
            return Ok(Some(mate_read));
        }
    }
    
    Ok(None)
}

/// Generalized function to add interval for any read record
fn add_read_interval(
    interval_trees: &mut HashMap<String, SortedVecIntervals>,
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

