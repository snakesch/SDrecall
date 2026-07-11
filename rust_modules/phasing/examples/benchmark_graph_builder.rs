use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use phasing::bam_reading::{build_allele_depth_map, migrate_bam_to_sorted_intervals_grouped};
use phasing::graph_builder::{build_phasing_graph, build_phasing_graph_legacy};
use phasing::structs::HaplotypeConfig;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Engine {
    Fast,
    Legacy,
}

#[derive(Debug, Parser)]
struct Args {
    #[arg(long)]
    bam: String,
    #[arg(long)]
    intrinsic_bam: String,
    #[arg(long)]
    reference: String,
    #[arg(long, value_enum, default_value_t = Engine::Fast)]
    engine: Engine,
    #[arg(long, default_value_t = 41)]
    mapq_cutoff: u8,
    #[arg(long, default_value_t = 15)]
    basequal_median_cutoff: u8,
    #[arg(long, default_value_t = 148.0)]
    mean_read_length: f32,
    #[arg(long, default_value_t = 25)]
    threads: u8,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    let total_start = Instant::now();

    let pairing_start = Instant::now();
    let (read_pairs, header) = migrate_bam_to_sorted_intervals_grouped(
        &args.bam,
        args.mapq_cutoff,
        args.basequal_median_cutoff,
        true,
        true,
        args.threads,
    )
    .map_err(|error| anyhow::anyhow!(error.to_string()))?;
    let pairing_time = pairing_start.elapsed();

    let primary_start = Instant::now();
    let primary = build_allele_depth_map(
        &args.bam,
        &args.reference,
        args.mapq_cutoff,
        args.basequal_median_cutoff,
    )
    .map_err(|error| anyhow::anyhow!(error.to_string()))?;
    let primary_time = primary_start.elapsed();

    let secondary_start = Instant::now();
    let secondary = build_allele_depth_map(
        &args.intrinsic_bam,
        &args.reference,
        args.mapq_cutoff,
        args.basequal_median_cutoff,
    )
    .map_err(|error| anyhow::anyhow!(error.to_string()))?;
    let secondary_time = secondary_start.elapsed();

    let config = HaplotypeConfig::new(args.mean_read_length);
    let build_start = Instant::now();
    let result = match args.engine {
        Engine::Fast => build_phasing_graph(&read_pairs, &primary, &secondary, &header, &config),
        Engine::Legacy => {
            build_phasing_graph_legacy(&read_pairs, &primary, &secondary, &header, &config)
        }
    }
    .map_err(|error| anyhow::anyhow!(error.to_string()))
    .context("graph construction failed")?;
    let build_time = build_start.elapsed();

    let digest_start = Instant::now();
    let (assignment_digest, edge_digest) = result.content_digests();
    let digest_time = digest_start.elapsed();

    println!("engine={:?}", args.engine);
    println!("read_pairs={}", read_pairs.readpair_dict.len());
    println!("vertices={}", result.vertex_count());
    println!("positive_edges={}", result.edge_count());
    println!("sparse_assignments={}", result.sparse_weight_entry_count());
    println!("excluded_labels={}", result.lowqual_qnames.len());
    println!("assignment_digest={assignment_digest:016x}");
    println!("edge_digest={edge_digest:016x}");
    println!("pairing_seconds={:.6}", pairing_time.as_secs_f64());
    println!("primary_support_seconds={:.6}", primary_time.as_secs_f64());
    println!(
        "secondary_support_seconds={:.6}",
        secondary_time.as_secs_f64()
    );
    println!("build_seconds={:.6}", build_time.as_secs_f64());
    println!("digest_seconds={:.6}", digest_time.as_secs_f64());
    println!("total_seconds={:.6}", total_start.elapsed().as_secs_f64());
    Ok(())
}
