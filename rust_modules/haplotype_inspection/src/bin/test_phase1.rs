/// Phase 1 integration test: validates the inspect_haplotypes main loop logic.
///
/// Tests the full Phase 1 data flow using a real BAM file:
/// 1. build_lapper_from_bam → BamLapperResult
/// 2. Option A: qname → &Vec<Record> index construction
/// 3. qname → chrom index construction
/// 4. Simulated haplotype assignment (split by HP tag or by even/odd index)
/// 5. Per-haplotype: read lookup → extract_continuous_regions_dict
/// 6. Per-region: record_hap_err_vectors_per_region → assemble_consensus
///    → judge_misalignment_by_extreme_vardensity → count_var
///
/// Output: structured TSV at /paedyl01/disk1/yangyxt/test_tmp/phase1_rust_test.tsv
///         debug log at /paedyl01/disk1/yangyxt/test_tmp/phase1_rust_test.log

use std::collections::{HashMap, HashSet};
use std::io::Write;
use ndarray::Array1;

use haplotype_inspection::bam_lappers::{build_lapper_from_bam, BamLapperResult};
use haplotype_inspection::pairwise_read_inspection::{
    read_id, count_var, ReadQseqData,
};
use haplotype_inspection::identify_misaligned_haps::{
    extract_continuous_regions_dict,
    record_hap_err_vectors_per_region,
    assemble_consensus,
    judge_misalignment_by_extreme_vardensity,
};
use rust_htslib::bam::Record;

const BAM_PATH: &str = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.clean.bam";
const INTRINSIC_BAM: &str = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/total_intrinsic_alignments.filtered.bam";
const OUTPUT_TSV: &str = "/paedyl01/disk1/yangyxt/test_tmp/phase1_rust_test.tsv";

fn main() {
    env_logger::init();

    eprintln!("=== Phase 1 Rust Integration Test ===");
    eprintln!("BAM: {}", BAM_PATH);

    // ── Step 1: Build Lapper ──────────────────────────────────────────────
    eprintln!("\n--- Step 1: Build Lapper from BAM ---");
    let bam_lapper = build_lapper_from_bam(BAM_PATH, 20, 10, true, true)
        .expect("Failed to build Lapper from BAM");
    eprintln!("Lapper built: {} chroms, {} qnames retained, {} noisy filtered",
              bam_lapper.lapper_dict.len(),
              bam_lapper.stats.qnames_retained,
              bam_lapper.stats.noisy_qnames_filtered);
    eprintln!("read_dict entries: {}", bam_lapper.read_dict.len());
    eprintln!("qname_idx_dict entries: {}", bam_lapper.qname_idx_dict.len());

    // Count total records across all qname_idx entries
    let total_records: usize = bam_lapper.read_dict.values().map(|v| v.len()).sum();
    eprintln!("Total BAM records in read_dict: {}", total_records);

    // ── Step 2: Build Option A index (qname → &Vec<Record>) ──────────────
    eprintln!("\n--- Step 2: Build qname→records index ---");
    let mut qname_to_records: HashMap<&str, &Vec<Record>> = HashMap::new();
    for (qname, &qname_idx) in &bam_lapper.qname_idx_dict {
        if let Some(records) = bam_lapper.read_dict.get(&qname_idx) {
            qname_to_records.insert(qname.as_str(), records);
        }
    }
    eprintln!("qname_to_records index: {} entries", qname_to_records.len());

    // ── Step 3: Build qname → chrom index ─────────────────────────────────
    eprintln!("\n--- Step 3: Build qname→chrom index ---");
    let mut qname_to_chrom: HashMap<&str, String> = HashMap::new();
    for (chrom, lapper) in &bam_lapper.lapper_dict {
        for iv in lapper.iter() {
            if let Some(qname) = bam_lapper.qname_dict.get(&iv.val) {
                qname_to_chrom.entry(qname.as_str()).or_insert_with(|| chrom.clone());
            }
        }
    }
    eprintln!("qname_to_chrom index: {} entries", qname_to_chrom.len());

    // Print first 5 entries
    let mut sorted_entries: Vec<_> = qname_to_chrom.iter().collect();
    sorted_entries.sort_by_key(|(k, _)| *k);
    for (qn, ch) in sorted_entries.iter().take(5) {
        eprintln!("  {} → {}", qn, ch);
    }

    // ── Step 4: Simulate haplotypes ───────────────────────────────────────
    // Split all qnames into 2 "haplotypes" by sorted order (even=hap0, odd=hap1)
    eprintln!("\n--- Step 4: Simulate 2-haplotype split ---");
    let mut all_qnames: Vec<String> = qname_to_records.keys().map(|s| s.to_string()).collect();
    all_qnames.sort();

    let mut hap_qname_info: HashMap<i32, Vec<String>> = HashMap::new();
    for (i, qn) in all_qnames.iter().enumerate() {
        let hid = if i % 2 == 0 { 0 } else { 1 };
        hap_qname_info.entry(hid).or_default().push(qn.clone());
    }

    for (&hid, qnames) in &hap_qname_info {
        eprintln!("  hap {}: {} qnames", hid, qnames.len());
    }

    // ── Step 5-6: Phase 1 main loop ─────────────────────────────────────
    eprintln!("\n--- Step 5-6: Phase 1 main loop ---");

    let total_lowqual_qnames: HashSet<String> = HashSet::new();
    let mut hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
    let mut err_cache: HashMap<String, Array1<f32>> = HashMap::new();

    let mut output_lines: Vec<String> = Vec::new();
    output_lines.push("hap_id\tn_qnames\tn_reads\tn_regions\tregion_chrom\tregion_start\tregion_end\tregion_n_reads\tconsensus_len\tvar_count\textreme_vard\tregion_max_density".to_string());

    let mut total_var_counts: HashMap<i32, i32> = HashMap::new();
    let mut total_extreme_vard: HashMap<i32, bool> = HashMap::new();

    for (&hid, qnames_raw) in hap_qname_info.iter() {
        let qnames: Vec<String> = qnames_raw.iter()
            .filter(|qn| !total_lowqual_qnames.contains(qn.as_str()))
            .cloned()
            .collect();

        let n_qnames = qnames.len();

        // Option A read lookup
        let reads: Vec<&Record> = qnames.iter()
            .flat_map(|qn| {
                qname_to_records.get(qn.as_str())
                    .into_iter()
                    .flat_map(|recs| recs.iter())
            })
            .collect();

        let n_reads = reads.len();

        // Log read IDs for detailed comparison
        let mut read_id_strs: Vec<String> = reads.iter().map(|r| read_id(r)).collect();
        read_id_strs.sort();
        eprintln!("hap {}: {} qnames, {} reads", hid, n_qnames, n_reads);
        eprintln!("  read_ids (first 10): {:?}", &read_id_strs[..std::cmp::min(10, read_id_strs.len())]);

        if reads.is_empty() {
            eprintln!("  SKIPPED: no reads");
            continue;
        }

        // Get chrom for this haplotype
        let hap_chrom = qnames.iter()
            .find_map(|qn| qname_to_chrom.get(qn.as_str()).cloned())
            .unwrap_or_else(|| "unknown".to_string());

        // Extract continuous regions
        let conregion_dict = extract_continuous_regions_dict(&reads);
        let n_regions = conregion_dict.len();
        eprintln!("  {} continuous regions on {}", n_regions, hap_chrom);

        // Per-region inner loop
        for ((span_start, span_end), ref read_indices) in &conregion_dict {
            let region_reads: Vec<&Record> = read_indices.iter().map(|&i| reads[i]).collect();
            if region_reads.is_empty() {
                continue;
            }
            let region_n_reads = region_reads.len();
            let span = (*span_start, *span_end);

            eprintln!("  region {}:{}-{}: {} reads", hap_chrom, span.0, span.1, region_n_reads);

            // a) Record haplotype/error vectors
            let (read_spans, hap_vectors, err_vectors) =
                record_hap_err_vectors_per_region(&region_reads, &mut hap_cache, &mut err_cache);

            eprintln!("    hap_vectors shape: {:?}, err_vectors shape: {:?}, read_spans shape: {:?}",
                      hap_vectors.dim(), err_vectors.dim(), read_spans.dim());

            // b) Assemble consensus
            let consensus_sequence = assemble_consensus(&hap_vectors, &err_vectors, &read_spans);
            let consensus_len = consensus_sequence.len();
            eprintln!("    consensus_len: {}", consensus_len);

            // c) Judge extreme variant density
            let (extreme_vard, region_max_density) =
                judge_misalignment_by_extreme_vardensity(&consensus_sequence);

            // d) Count variants
            let var_count_val = count_var(&consensus_sequence);

            *total_var_counts.entry(hid).or_insert(0) += var_count_val;
            let prev_extreme = *total_extreme_vard.get(&hid).unwrap_or(&false);
            total_extreme_vard.insert(hid, extreme_vard || prev_extreme);

            eprintln!("    var_count={}, extreme_vard={}, max_density={:.6}",
                      var_count_val, extreme_vard, region_max_density);

            output_lines.push(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
                hid, n_qnames, n_reads, n_regions,
                hap_chrom, span.0, span.1, region_n_reads,
                consensus_len, var_count_val, extreme_vard, region_max_density));
        }
    }

    // ── Summary ──────────────────────────────────────────────────────────
    eprintln!("\n=== Summary ===");
    eprintln!("Hap cache size: {}", hap_cache.len());
    eprintln!("Err cache size: {}", err_cache.len());
    for (&hid, &vc) in total_var_counts.iter() {
        eprintln!("  hap {} total_var_count={} extreme_vard={}", hid, vc,
                  total_extreme_vard.get(&hid).unwrap_or(&false));
    }

    // ── Step 7: Build Lapper for intrinsic BAM ───────────────────────────
    eprintln!("\n--- Step 7: Build Lapper from intrinsic BAM ---");
    let intrin_lapper = build_lapper_from_bam(INTRINSIC_BAM, 0, 0, false, false)
        .expect("Failed to build intrinsic Lapper");
    eprintln!("Intrinsic Lapper built: {} chroms, {} qnames",
              intrin_lapper.lapper_dict.len(), intrin_lapper.stats.qnames_retained);

    // ── Write TSV output ─────────────────────────────────────────────────
    let mut f = std::fs::File::create(OUTPUT_TSV).expect("Failed to create output TSV");
    for line in &output_lines {
        writeln!(f, "{}", line).expect("Failed to write");
    }
    eprintln!("\nOutput written to {}", OUTPUT_TSV);
    eprintln!("\n=== Phase 1 Rust Integration Test PASSED ===");
}
