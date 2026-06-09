#!/usr/bin/env python3
"""
Diagnostic: dump the Python total_record_df and LP inputs/outputs for a single chunk,
so we can compare with Rust's intermediate data.
"""
import os, sys, logging, tempfile
import pandas as pd
import numpy as np

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

logging.basicConfig(level=logging.WARNING, format="[%(asctime)s] %(message)s")
diag_logger = logging.getLogger("diag")
diag_logger.setLevel(logging.INFO)

from fp_control.bam_ncls import migrate_bam_to_ncls, calculate_mean_read_length
from fp_control.graph_build import build_phasing_graph
from fp_control.phasing import phasing_realigned_reads
from fp_control.identify_misaligned_haps import inspect_by_haplotypes
from haplotype_inspection import inspect_haplotypes_rust

# ── Config ────────────────────────────────────────────────────
CHUNK_ID = int(sys.argv[1]) if len(sys.argv) > 1 else 37
BAM_DIR = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results"
INTRIN_DIR = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall"
REF_GENOME = "/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta"
OUTDIR = "/paedyl01/disk1/yangyxt/test_tmp/diag"
os.makedirs(OUTDIR, exist_ok=True)

bam = f"{BAM_DIR}/HG002.pooled.raw.deduped.{CHUNK_ID}.bam"
intrin_bam = f"{INTRIN_DIR}/total_intrinsic_alignments.filtered.{CHUNK_ID}.bam"

import numba; numba.set_num_threads(4)

quiet = logging.getLogger("quiet")
quiet.setLevel(logging.WARNING)
quiet.addHandler(logging.NullHandler())

# ── Step 1: shared setup ──────────────────────────────────────
diag_logger.info(f"Chunk {CHUNK_ID}: shared setup...")
mean_rl = calculate_mean_read_length(bam)
(phased_graph, weight_matrix, qname_to_node,
 total_readhap_vector, total_readerr_vector, read_ref_pos_dict,
 total_lowqual_qnames, node_read_ids, read_id_read_dict) = \
    build_phasing_graph(bam, mean_rl, set(), REF_GENOME, mapq_filter=10,
                        basequal_median_filter=15, threads=4, logger=quiet)

qname_hap_info, hap_qname_info = phasing_realigned_reads(
    phased_graph, weight_matrix, 0.301,
    total_readhap_vector=total_readhap_vector,
    total_readerr_vector=total_readerr_vector,
    node_read_ids=node_read_ids, logger=quiet)

diag_logger.info(f"  {len(hap_qname_info)} haplotypes, {sum(len(v) for v in hap_qname_info.values())} total qnames")

# ── Step 2: Python path with TSV dump ────────────────────────
diag_logger.info("Python path...")
py_meta = f"{OUTDIR}/chunk{CHUNK_ID}_python_meta.tsv"
bam_ncls = migrate_bam_to_ncls(bam, mapq_filter=10, basequal_median_filter=15, logger=quiet)
intrin_ncls = migrate_bam_to_ncls(intrin_bam, mapq_filter=0, basequal_median_filter=0,
                                   paired=False, filter_noisy=False, logger=quiet)
intrin_bam_ncls = intrin_ncls[:-1]

py_correct, py_mismap = inspect_by_haplotypes(
    bam, bam_ncls, hap_qname_info, qname_hap_info,
    read_id_read_dict, node_read_ids, intrin_bam_ncls,
    qname_to_node, total_lowqual_qnames,
    total_readhap_vector, total_readerr_vector, {},
    read_ref_pos_dict,
    compare_haplotype_meta_tab=py_meta,
    mean_read_length=mean_rl, logger=quiet)
py_correct, py_mismap = set(py_correct), set(py_mismap)

# ── Step 3: Rust path with TSV dump ──────────────────────────
diag_logger.info("Rust path...")
rust_meta = f"{OUTDIR}/chunk{CHUNK_ID}_rust_meta.tsv"
rust_correct_list, rust_mismap_list = inspect_haplotypes_rust(
    bam_path=bam, intrinsic_bam_path=intrin_bam,
    hap_qname_info={k: list(v) for k, v in hap_qname_info.items()},
    qname_hap_info=dict(qname_hap_info),
    qname_to_node=dict(qname_to_node),
    total_lowqual_qnames=list(total_lowqual_qnames),
    compare_haplotype_meta_tab=rust_meta,
    mean_read_length=float(mean_rl),
    recall_mq_cutoff=10, basequal_median_cutoff=15)
rust_correct, rust_mismap = set(rust_correct_list), set(rust_mismap_list)

# ── Step 4: Compare ──────────────────────────────────────────
diag_logger.info(f"\nPython: {len(py_correct)} correct, {len(py_mismap)} mismap")
diag_logger.info(f"Rust:   {len(rust_correct)} correct, {len(rust_mismap)} mismap")
diag_logger.info(f"correct_only_python: {len(py_correct - rust_correct)}")
diag_logger.info(f"correct_only_rust:   {len(rust_correct - py_correct)}")
diag_logger.info(f"mismap_only_python:  {len(py_mismap - rust_mismap)}")
diag_logger.info(f"mismap_only_rust:    {len(rust_mismap - py_mismap)}")

# ── Step 5: Compare the meta TSVs ────────────────────────────
py_raw = py_meta.replace(".tsv", ".raw.tsv")
if os.path.exists(py_raw):
    py_df = pd.read_csv(py_raw, sep="\t")
    diag_logger.info(f"\nPython raw meta ({py_raw}): {len(py_df)} rows")
    diag_logger.info(f"Columns: {list(py_df.columns)}")
    diag_logger.info(f"\n{py_df.to_string()}")
else:
    diag_logger.warning(f"Python raw meta not found: {py_raw}")

if os.path.exists(rust_meta):
    rust_df = pd.read_csv(rust_meta, sep="\t")
    diag_logger.info(f"\nRust meta ({rust_meta}): {len(rust_df)} rows")
    diag_logger.info(f"Columns: {list(rust_df.columns)}")
    diag_logger.info(f"\n{rust_df.to_string()}")
else:
    diag_logger.warning(f"Rust meta not found: {rust_meta}")

# ── Step 6: Side-by-side per-haplotype comparison ────────────
if os.path.exists(py_raw) and os.path.exists(rust_meta):
    py_df = pd.read_csv(py_raw, sep="\t")
    rust_df = pd.read_csv(rust_meta, sep="\t")

    # Compare hap_ids present
    py_hids = set(py_df["hap_id"].unique()) if "hap_id" in py_df.columns else set()
    rust_hids = set(rust_df["hap_id"].unique()) if "hap_id" in rust_df.columns else set()
    diag_logger.info(f"\nPython hap_ids ({len(py_hids)}): {sorted(py_hids)}")
    diag_logger.info(f"Rust hap_ids ({len(rust_hids)}): {sorted(rust_hids)}")
    diag_logger.info(f"Only in Python: {sorted(py_hids - rust_hids)}")
    diag_logger.info(f"Only in Rust:   {sorted(rust_hids - py_hids)}")

    # For shared haplotypes, compare key columns
    common = sorted(py_hids & rust_hids)
    if common and "coefficient" in py_df.columns and "coefficient" in rust_df.columns:
        diag_logger.info(f"\nPer-haplotype coefficient comparison (shared hids):")
        for hid in common:
            py_rows = py_df[py_df["hap_id"] == hid]
            rs_rows = rust_df[rust_df["hap_id"] == hid]
            py_coeff = py_rows["coefficient"].values
            rs_coeff = rs_rows["coefficient"].values
            py_rank = py_rows["rank"].values if "rank" in py_rows.columns else []
            rs_rank = rs_rows["rank"].values if "rank" in rs_rows.columns else rs_rows["varc_rank"].values if "varc_rank" in rs_rows.columns else []
            diag_logger.info(f"  hid={hid}: py_coeff={py_coeff} rs_coeff={rs_coeff}  py_rank={py_rank} rs_rank={rs_rank}")

diag_logger.info(f"\nDiagnostic files saved to {OUTDIR}/chunk{CHUNK_ID}_*.tsv")
