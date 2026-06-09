/// Standalone test program for BAM to Lapper functionality
/// Tests build_lapper_from_bam and query functions independently
/// Does NOT depend on lib.rs, avoiding compilation of stub modules
///
/// Usage:
///   test_bam_lapper <bam_file> [chrom:start-end,...]
///   test_bam_lapper input.bam chr1:10000-20000

// Include bam_lapper directly — bypasses lib.rs so stub modules are not compiled
#[path = "../bam_lappers.rs"]
mod bam_lappers;

use bam_lappers::{
    build_lapper_from_bam, query_overlapping_qname_indices, query_overlapping_reads,
    query_overlapping_read_pairs,
};
use std::env;
use std::io::Write;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <bam_file> [chrom:start-end,...]", args[0]);
        eprintln!("\nExample:");
        eprintln!("  {} input.bam", args[0]);
        eprintln!("  {} input.bam chr1:1000-2000,chr2:5000-6000", args[0]);
        std::process::exit(1);
    }

    let bam_path = &args[1];

    println!("=== BAM to Lapper Test ===");
    println!("BAM file: {}", bam_path);
    println!();

    // No quality filters — only skip secondary/supplementary/duplicate
    // This matches: samtools view -F 3328 (0x100+0x400+0x800)
    let mapq_filter = 0;
    let basequal_median_filter = 0;
    let paired = false;
    let filter_noisy = false;

    println!("Parameters:");
    println!("  MAPQ filter: {}", mapq_filter);
    println!("  Base quality median filter: {}", basequal_median_filter);
    println!("  Paired-end mode: {}", paired);
    println!("  Filter noisy reads: {}", filter_noisy);
    println!();

    // Build Lapper from BAM
    println!("Building Lapper from BAM...");
    let start_time = std::time::Instant::now();

    let result = match build_lapper_from_bam(
        bam_path,
        mapq_filter,
        basequal_median_filter,
        paired,
        filter_noisy,
    ) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error building Lapper: {}", e);
            std::process::exit(1);
        }
    };

    let elapsed = start_time.elapsed();
    println!("Lapper built in {:.2}s", elapsed.as_secs_f64());
    println!();

    // Print statistics
    println!("=== Processing Statistics ===");
    println!("Total reads processed: {}", result.stats.total_reads_processed);
    println!("Skipped alignments (secondary/supplementary/duplicate): {}", result.stats.skipped_alignments);
    println!("Qnames retained: {}", result.stats.qnames_retained);
    println!("Noisy qnames filtered: {}", result.stats.noisy_qnames_filtered);
    println!();

    // Print chromosome coverage
    println!("=== Chromosome Coverage ===");
    let mut chroms: Vec<_> = result.lapper_dict.keys().collect();
    chroms.sort();

    for chrom in &chroms {
        if let Some(lapper) = result.lapper_dict.get(*chrom) {
            let interval_count = lapper.len();
            if interval_count > 0 {
                println!("  {}: {} intervals", chrom, interval_count);
            }
        }
    }
    println!();

    // Parse test intervals from command line or auto-detect
    let test_intervals = if args.len() > 2 {
        parse_intervals(&args[2])
    } else {
        // Auto-detect: pick first chromosome with data, use first interval region
        let mut default_intervals = Vec::new();
        for chrom in &chroms {
            if let Some(lapper) = result.lapper_dict.get(*chrom) {
                if lapper.len() > 0 {
                    if let Some(first_iv) = lapper.iter().next() {
                        let start = first_iv.start;
                        let end = first_iv.stop;
                        default_intervals.push((
                            chrom.to_string(),
                            start.saturating_sub(1000),
                            end + 1000,
                        ));
                        break;
                    }
                }
            }
        }
        default_intervals
    };

    if test_intervals.is_empty() {
        println!("No test intervals specified and no data found in BAM.");
        return;
    }

    // Test queries
    println!("=== Testing Interval Queries ===");
    for (chrom, start, end) in &test_intervals {
        println!("\nQuery: {}:{}-{}", chrom, start, end);

        // Test 1: Query qname indices from Lapper intervals
        if let Some(lapper) = result.lapper_dict.get(chrom) {
            let qname_indices = query_overlapping_qname_indices(lapper, *start, *end);
            println!("  Lapper interval hits: {} qname indices", qname_indices.len());

            // Deduplicate (Lapper may return duplicates)
            let mut unique_indices: Vec<u32> = qname_indices.clone();
            unique_indices.sort();
            unique_indices.dedup();
            println!("  Unique qname indices: {}", unique_indices.len());

            // Print ALL qnames sorted (for cross-checking with samtools)
            println!("  --- BEGIN QNAME LIST ---");
            let mut sorted_qnames: Vec<String> = unique_indices
                .iter()
                .filter_map(|idx| result.qname_dict.get(idx).cloned())
                .collect();
            sorted_qnames.sort();
            for qname in &sorted_qnames {
                println!("  QNAME: {}", qname);
            }
            println!("  --- END QNAME LIST ---");
        } else {
            println!("  Chromosome not found in Lapper");
            continue;
        }

        // Test 2: Query overlapping reads (with actual overlap verification)
        let reads = query_overlapping_reads(
            &result.lapper_dict,
            &result.read_dict,
            chrom,
            *start,
            *end,
        );
        println!("  Overlapping reads (overlap-verified): {}", reads.len());

        // Count unique qnames from the actual overlapping reads
        let mut read_qnames: Vec<String> = reads
            .iter()
            .map(|r| String::from_utf8_lossy(r.qname()).to_string())
            .collect();
        read_qnames.sort();
        read_qnames.dedup();
        println!("  Unique qnames from overlapping reads: {}", read_qnames.len());

        // Print read details for cross-checking
        println!("  --- BEGIN READ DETAILS ---");
        for read in &reads {
            let qname = String::from_utf8_lossy(read.qname());
            let read_start = read.pos();
            let read_end = read.cigar().end_pos();
            let mapq = read.mapq();
            let flag = read.flags();
            println!(
                "  READ: {}\t{}\t{}\t{}\tMAPQ={}\tFLAG={}",
                qname, chrom, read_start, read_end, mapq, flag
            );
        }
        println!("  --- END READ DETAILS ---");
    }

    // Flush stdout
    std::io::stdout().flush().unwrap();
    println!("\n=== Test Complete ===");
}

/// Parse interval strings like "chr1:1000-2000,chr2:5000-6000"
fn parse_intervals(interval_str: &str) -> Vec<(String, u32, u32)> {
    let mut intervals = Vec::new();
    for part in interval_str.split(',') {
        let part = part.trim();
        if let Some((chrom, range)) = part.split_once(':') {
            if let Some((start_str, end_str)) = range.split_once('-') {
                if let (Ok(start), Ok(end)) = (start_str.parse::<u32>(), end_str.parse::<u32>()) {
                    intervals.push((chrom.to_string(), start, end));
                }
            }
        }
    }
    intervals
}
