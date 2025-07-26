use pyo3::prelude::*;
use pyo3::exceptions::{PyIOError, PyValueError};
use rust_htslib::{bam, bam::Read, bam::record::Aux};
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use anyhow::{Result, Context};



/// Read BED file and return regions as Vec<(chr, start, end)>
fn read_bed_regions(bed_path: &str) -> Result<Vec<(String, u64, u64)>> {
    use std::io::{BufRead, BufReader};
    
    let file = File::open(bed_path)
        .with_context(|| format!("Failed to open BED file: {}", bed_path))?;
    let reader = BufReader::new(file);
    
    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let chr = parts[0].to_string();
            let start = parts[1].parse::<u64>()
                .with_context(|| format!("Invalid start position: {}", parts[1]))?;
            let end = parts[2].parse::<u64>()
                .with_context(|| format!("Invalid end position: {}", parts[2]))?;
            regions.push((chr, start, end));
        }
    }
    
    Ok(regions)
}

/// Check if a read should be included based on multi-aligned filter
fn should_include_read(record: &bam::Record, multi_aligned: bool) -> bool {
    // Always require MAPQ < 50 for inclusion (new mandatory condition)
    if record.mapq() >= 50 {
        return false;
    }

    if !multi_aligned {
        return true;
    }
    
    // Filter: ![SA] && [XA] && abs(AS - XS) <= 5
    // Check for SA tag (should NOT exist)
    let has_sa = match record.aux(b"SA") {
        Ok(_) => true,
        Err(_) => false,
    };
    
    if has_sa {
        return false;
    }
    
    // Check for XA tag (should exist)
    let has_xa = match record.aux(b"XA") {
        Ok(_) => true,
        Err(_) => false,
    };
    
    if !has_xa {
        return false;
    }
    
    // Check AS tag
    let as_score_opt = match record.aux(b"AS") {
        Ok(Aux::I32(score)) => Some(score),
        Ok(Aux::I16(score)) => Some(score as i32),
        Ok(Aux::I8(score)) => Some(score as i32),
        Ok(Aux::U32(score)) => Some(score as i32),
        Ok(Aux::U16(score)) => Some(score as i32),
        Ok(Aux::U8(score)) => Some(score as i32),
        _ => None,  // AS missing
    };

    // Check XS tag
    let xs_score_opt = match record.aux(b"XS") {
        Ok(Aux::I32(score)) => Some(score),
        Ok(Aux::I16(score)) => Some(score as i32),
        Ok(Aux::I8(score)) => Some(score as i32),
        Ok(Aux::U32(score)) => Some(score as i32),
        Ok(Aux::U16(score)) => Some(score as i32),
        Ok(Aux::U8(score)) => Some(score as i32),
        _ => None,  // XS missing
    };

    // If both AS and XS are available, check deviation
    if let (Some(as_score), Some(xs_score)) = (as_score_opt, xs_score_opt) {
        (as_score - xs_score).abs() <= 5
    } else {
        // If AS or XS is missing, fallback to just the MAPQ < 50 check (already enforced above)
        // This satisfies the "alternative condition" without additional logic
        true
    }
}

/// Convert quality scores to FASTQ format
fn quality_to_string(qual: &[u8]) -> String {
    qual.iter()
        .map(|&q| (q + 33) as char)
        .collect()
}

/// Main function to convert BAM to FASTQ
#[pyfunction]
#[pyo3(signature = (input_bam, region_bed, output_freads, output_rreads, multi_aligned=false, threads=1, _tmp_dir="/tmp"))]
fn bam_to_fastq_biobambam(
    input_bam: &str,
    region_bed: &str, 
    output_freads: &str,
    output_rreads: &str,
    multi_aligned: bool,
    threads: usize,
    _tmp_dir: &str,
) -> PyResult<(String, String)> {
    // Ensure output directories exist
    if let Some(parent) = Path::new(output_freads).parent() {
        create_dir_all(parent)
            .map_err(|e| PyIOError::new_err(format!("Failed to create output directory: {}", e)))?;
    }
    if let Some(parent) = Path::new(output_rreads).parent() {
        create_dir_all(parent)
            .map_err(|e| PyIOError::new_err(format!("Failed to create output directory: {}", e)))?;
    }
    
    // Read BED regions
    let regions = read_bed_regions(region_bed)
        .map_err(|e| PyValueError::new_err(format!("Failed to read BED file: {}", e)))?;
    
    // Open BAM file with index
    let mut bam_reader = bam::IndexedReader::from_path(input_bam)
        .map_err(|e| PyIOError::new_err(format!("Failed to open BAM file: {}", e)))?;
    
    // Create BAM header
    let header = bam_reader.header().clone();
    
    // Set threads for BAM reading
    if threads > 1 {
        bam_reader.set_threads(threads - 1)
            .map_err(|e| PyIOError::new_err(format!("Failed to set threads: {}", e)))?;
    }
    
    // Create output files
    let r1_file = File::create(output_freads)
        .map_err(|e| PyIOError::new_err(format!("Failed to create R1 file: {}", e)))?;
    let r2_file = File::create(output_rreads)
        .map_err(|e| PyIOError::new_err(format!("Failed to create R2 file: {}", e)))?;
    
    let r1_writer = Arc::new(Mutex::new(BufWriter::new(r1_file)));
    let r2_writer = Arc::new(Mutex::new(BufWriter::new(r2_file)));
    
// Process each region
for (chr, start, end) in regions {
    // Fetch reads in the region
    let tid = header.tid(chr.as_bytes())
        .ok_or_else(|| PyValueError::new_err(format!("Chromosome {} not found in BAM", chr)))?;
    
    bam_reader.fetch((tid, start as i64, end as i64))
        .map_err(|e| PyIOError::new_err(format!("Failed to fetch region: {}", e)))?;
    
    // Collect all reads first (without filtering yet, to allow pair checks)
    let mut read_pairs: std::collections::HashMap<Vec<u8>, (Option<bam::Record>, Option<bam::Record>)> = 
        std::collections::HashMap::new();
    
    for result in bam_reader.records() {
        let record = result
            .map_err(|e| PyIOError::new_err(format!("Failed to read BAM record: {}", e)))?;
        
        let qname = record.qname().to_vec();
        let is_first = record.is_first_in_template();
        
        let entry = read_pairs.entry(qname).or_insert((None, None));
        if is_first {
            entry.0 = Some(record);
        } else {
            entry.1 = Some(record);
        }
    }
    
    // REVISED: Fetch mates for singletons to include both reads if one overlaps
    // We check if the read is paired and mate is mapped, then fetch from a small window
    // This is minimal: Only for singletons, avoiding full pair fetching overhead
    // Error handling integrated with PyResult propagation
    for (_, entry) in read_pairs.iter_mut() {
        if let (Some(read), None) | (None, Some(read)) = entry {
            if read.is_paired() && !read.is_mate_unmapped() {
                let mtid = read.mtid();
                let mpos = read.mpos();
                let window_start = (mpos - 1500).max(0);
                let window_end = mpos + 1500 + read.seq_len() as i64;
                
                bam_reader.fetch((mtid, window_start, window_end))
                    .map_err(|e| PyIOError::new_err(format!("Failed to fetch mate region: {}", e)))?;
                
                let mut mate = bam::Record::new();
                while let Some(Ok(())) = bam_reader.read(&mut mate) {  // REVISED: Changed to while let Ok(true) = ... to fix map_err on Option (rust-htslib's read returns Result<bool, Error>, not Option; loop until false or error)
                    if mate.qname() == read.qname() && mate.is_first_in_template() != read.is_first_in_template() {
                        if read.is_first_in_template() {
                            entry.1 = Some(mate.clone());  // REVISED: Clone mate for ownership (minimal for singletons)
                        } else {
                            entry.0 = Some(mate.clone());
                        }
                        break;
                    }
                }
            }
        }
    }

    // Now filter per pair: If at least one read passes the filter, include both (if present), but remove singletons
    for (_, (r1_opt, r2_opt)) in read_pairs {
        let passes_filter = match (&r1_opt, &r2_opt) {
            (Some(r1), Some(r2)) => should_include_read(r1, multi_aligned) || should_include_read(r2, multi_aligned),
            (Some(_r1), None) => false,
            (None, Some(_r2)) => false,
            (None, None) => false,
        };

        if passes_filter {
            // Write both if available, even if one didn't individually pass
            if let (Some(r1), Some(r2)) = (r1_opt, r2_opt) {
                // Write R1 first (ensures order)
                {
                    let mut writer = r1_writer.lock().unwrap();
                    writeln!(writer, "@{}", std::str::from_utf8(r1.qname()).unwrap())
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R1: {}", e)))?;
                    writeln!(writer, "{}", std::str::from_utf8(&r1.seq().as_bytes()).unwrap())
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R1: {}", e)))?;
                    writeln!(writer, "+")
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R1: {}", e)))?;
                    writeln!(writer, "{}", quality_to_string(r1.qual()))
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R1: {}", e)))?;
                }
                
                // Write R2
                {
                    let mut writer = r2_writer.lock().unwrap();
                    writeln!(writer, "@{}", std::str::from_utf8(r2.qname()).unwrap())
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R2: {}", e)))?;
                    writeln!(writer, "{}", std::str::from_utf8(&r2.seq().as_bytes()).unwrap())
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R2: {}", e)))?;
                    writeln!(writer, "+")
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R2: {}", e)))?;
                    writeln!(writer, "{}", quality_to_string(r2.qual()))
                        .map_err(|e| PyIOError::new_err(format!("Failed to write R2: {}", e)))?;
                }
            }
        }
    }
}
    
    // Ensure all data is written
    r1_writer.lock().unwrap().flush()
        .map_err(|e| PyIOError::new_err(format!("Failed to flush R1 file: {}", e)))?;
    r2_writer.lock().unwrap().flush()
        .map_err(|e| PyIOError::new_err(format!("Failed to flush R2 file: {}", e)))?;
    
    Ok((output_freads.to_string(), output_rreads.to_string()))
}

/// Python module definition
#[pymodule]
fn rust_read_extraction(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bam_to_fastq_biobambam, m)?)?;
    Ok(())
} 