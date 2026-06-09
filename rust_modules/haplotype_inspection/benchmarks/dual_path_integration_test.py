#!/usr/bin/env python3
"""
Dual-path integration test: Rust inspect_haplotypes_rust vs Python inspect_by_haplotypes.

Runs both implementations on real HG002 chunk BAMs and compares outputs.
Produces a CSV summary and console report with timing and concordance.

Usage:
    cd /paedyl01/disk1/yangyxt/SDrecall-rust-migration
    python rust_modules/haplotype_inspection/benchmarks/dual_path_integration_test.py [--n_chunks 100]
"""

import os
import sys
import time
import csv
import argparse
import logging
import traceback
import tempfile
import numpy as np

# Ensure project root is on path
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

# Configure logging
logging.basicConfig(
    level=logging.WARNING,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
test_logger = logging.getLogger("dual_path_test")
test_logger.setLevel(logging.INFO)

# Silence noisy sub-loggers during test
for name in ("build_phasing_graph", "SubProcess", "pyo3"):
    logging.getLogger(name).setLevel(logging.WARNING)

# ── Imports from SDrecall ──────────────────────────────────────────────
from fp_control.bam_ncls import migrate_bam_to_ncls, calculate_mean_read_length
from fp_control.graph_build import build_phasing_graph
from fp_control.phasing import phasing_realigned_reads
from fp_control.identify_misaligned_haps import inspect_by_haplotypes
from haplotype_inspection import inspect_haplotypes_rust

# ── Paths ──────────────────────────────────────────────────────────────
BAM_DIR = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results"
INTRIN_DIR = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall"
REF_GENOME = "/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta"

# Production parameters (from subprocess logs)
RECALL_MQ_CUTOFF = 10
BASEQUAL_MEDIAN_CUTOFF = 15
EDGE_WEIGHT_CUTOFF = 0.301
THREADS = 4

# ── Output ─────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_CSV = os.path.join(SCRIPT_DIR, "dual_path_results.csv")


def select_chunks(n: int, total: int = 1840) -> list:
    """Select n evenly-spaced chunk IDs from 1..total."""
    if n >= total:
        return list(range(1, total + 1))
    step = total / n
    return [int(1 + i * step) for i in range(n)]


def bam_path(chunk_id: int) -> str:
    return os.path.join(BAM_DIR, f"HG002.pooled.raw.deduped.{chunk_id}.bam")


def intrinsic_bam_path(chunk_id: int) -> str:
    return os.path.join(INTRIN_DIR, f"total_intrinsic_alignments.filtered.{chunk_id}.bam")


def run_chunk(chunk_id: int, quiet_logger) -> dict:
    """Run both Python and Rust paths on one chunk. Returns a result dict."""
    bam = bam_path(chunk_id)
    intrin_bam = intrinsic_bam_path(chunk_id)

    result = {
        "chunk_id": chunk_id,
        "status": "unknown",
        "n_haplotypes": 0,
        "n_qnames": 0,
        "python_time_s": None,
        "rust_time_s": None,
        "speedup": None,
        "correct_match": None,
        "mismap_match": None,
        "correct_only_python": 0,
        "correct_only_rust": 0,
        "mismap_only_python": 0,
        "mismap_only_rust": 0,
    }

    # Verify files exist
    if not os.path.exists(bam):
        result["status"] = "bam_missing"
        return result
    if not os.path.exists(intrin_bam):
        result["status"] = "intrin_bam_missing"
        return result

    # ── Step 1: Shared setup ──────────────────────────────────────
    mean_read_length = calculate_mean_read_length(bam)

    # Build phasing graph (Rust-accelerated, shared by both paths)
    graph_result = build_phasing_graph(
        bam, mean_read_length, set(), REF_GENOME,
        mapq_filter=RECALL_MQ_CUTOFF,
        basequal_median_filter=BASEQUAL_MEDIAN_CUTOFF,
        edge_weight_cutoff=EDGE_WEIGHT_CUTOFF,
        threads=THREADS,
        logger=quiet_logger,
    )

    phased_graph = graph_result[0]
    if phased_graph is None:
        result["status"] = "graph_failed"
        return result

    (phased_graph, weight_matrix, qname_to_node,
     total_readhap_vector, total_readerr_vector, read_ref_pos_dict,
     total_lowqual_qnames, node_read_ids, read_id_read_dict) = graph_result

    # Phase reads into haplotypes
    qname_hap_info, hap_qname_info = phasing_realigned_reads(
        phased_graph, weight_matrix, EDGE_WEIGHT_CUTOFF,
        total_readhap_vector=total_readhap_vector,
        total_readerr_vector=total_readerr_vector,
        node_read_ids=node_read_ids,
        logger=quiet_logger,
    )

    n_haps = len(hap_qname_info)
    result["n_haplotypes"] = n_haps

    if n_haps <= 2:
        result["status"] = "skipped_le2_haps"
        return result

    # Count total qnames across all haplotypes
    all_qnames = set()
    for qns in hap_qname_info.values():
        all_qnames.update(qns)
    result["n_qnames"] = len(all_qnames)

    # ── Step 2: Python NCLS setup (needed only for Python path) ───
    bam_ncls_result = migrate_bam_to_ncls(
        bam, mapq_filter=RECALL_MQ_CUTOFF,
        basequal_median_filter=BASEQUAL_MEDIAN_CUTOFF,
        logger=quiet_logger,
    )
    if bam_ncls_result is None:
        result["status"] = "ncls_failed"
        return result
    bam_ncls = bam_ncls_result

    intrin_ncls_result = migrate_bam_to_ncls(
        intrin_bam, mapq_filter=0, basequal_median_filter=0,
        paired=False, filter_noisy=False,
        logger=quiet_logger,
    )
    if intrin_ncls_result is None:
        result["status"] = "intrin_ncls_failed"
        return result
    # Strip noisy_qnames from intrin_bam_ncls (matches production line 244)
    intrin_bam_ncls = intrin_ncls_result[:-1]

    # ── Step 3: Python path ───────────────────────────────────────
    # Python inspect_by_haplotypes writes CSVs to compare_haplotype_meta_tab;
    # provide a real temp path so it doesn't fail on empty string.
    py_meta_tab = tempfile.mktemp(suffix=f".chunk{chunk_id}.haplotype_meta.tsv")
    t0 = time.perf_counter()
    py_correct, py_mismap = inspect_by_haplotypes(
        bam,
        bam_ncls,
        hap_qname_info,
        qname_hap_info,
        read_id_read_dict,
        node_read_ids,
        intrin_bam_ncls,
        qname_to_node,
        total_lowqual_qnames,
        total_readhap_vector,
        total_readerr_vector,
        {},  # total_genomic_haps
        read_ref_pos_dict,
        compare_haplotype_meta_tab=py_meta_tab,
        mean_read_length=mean_read_length,
        logger=quiet_logger,
    )
    py_time = time.perf_counter() - t0
    py_correct = set(py_correct)
    py_mismap = set(py_mismap)

    # ── Step 4: Rust path ─────────────────────────────────────────
    t0 = time.perf_counter()
    rust_correct_list, rust_mismap_list = inspect_haplotypes_rust(
        bam_path=bam,
        intrinsic_bam_path=intrin_bam,
        hap_qname_info={k: list(v) for k, v in hap_qname_info.items()},
        qname_hap_info=dict(qname_hap_info),
        qname_to_node=dict(qname_to_node),
        total_lowqual_qnames=list(total_lowqual_qnames),
        compare_haplotype_meta_tab="",
        mean_read_length=float(mean_read_length),
        recall_mq_cutoff=RECALL_MQ_CUTOFF,
        basequal_median_cutoff=BASEQUAL_MEDIAN_CUTOFF,
    )
    rust_time = time.perf_counter() - t0
    rust_correct = set(rust_correct_list)
    rust_mismap = set(rust_mismap_list)

    # ── Step 5: Compare ───────────────────────────────────────────
    correct_match = py_correct == rust_correct
    mismap_match = py_mismap == rust_mismap

    result.update({
        "status": "ok",
        "python_time_s": round(py_time, 4),
        "rust_time_s": round(rust_time, 4),
        "speedup": round(py_time / rust_time, 2) if rust_time > 0 else None,
        "correct_match": correct_match,
        "mismap_match": mismap_match,
        "correct_only_python": len(py_correct - rust_correct),
        "correct_only_rust": len(rust_correct - py_correct),
        "mismap_only_python": len(py_mismap - rust_mismap),
        "mismap_only_rust": len(rust_mismap - py_mismap),
    })

    # Clean up temp files from Python path
    for suffix in ("", ".raw"):
        p = py_meta_tab.replace(".tsv", f"{suffix}.tsv")
        if os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                pass

    if not correct_match or not mismap_match:
        result["diff_detail"] = {
            "correct_only_python": sorted(py_correct - rust_correct),
            "correct_only_rust": sorted(rust_correct - py_correct),
            "mismap_only_python": sorted(py_mismap - rust_mismap),
            "mismap_only_rust": sorted(rust_mismap - py_mismap),
        }

    return result


def main():
    parser = argparse.ArgumentParser(description="Dual-path Rust vs Python integration test")
    parser.add_argument("--n_chunks", type=int, default=100, help="Number of chunks to test (default: 100)")
    parser.add_argument("--start", type=int, default=None, help="Start from a specific chunk ID")
    parser.add_argument("--end", type=int, default=None, help="End at a specific chunk ID")
    args = parser.parse_args()

    # Configure numba
    try:
        import numba
        numba.set_num_threads(THREADS)
    except Exception:
        pass

    # Select chunks
    if args.start is not None and args.end is not None:
        chunk_ids = list(range(args.start, args.end + 1))
    else:
        chunk_ids = select_chunks(args.n_chunks)

    test_logger.info(f"Testing {len(chunk_ids)} chunks: {chunk_ids[:5]}...{chunk_ids[-5:]}")
    test_logger.info(f"Output CSV: {OUTPUT_CSV}")

    # Quiet logger for sub-functions
    quiet_logger = logging.getLogger("quiet")
    quiet_logger.setLevel(logging.WARNING)
    if not quiet_logger.handlers:
        quiet_logger.addHandler(logging.NullHandler())

    results = []
    n_ok = 0
    n_match = 0
    n_mismatch = 0
    n_skipped = 0
    n_error = 0
    py_times = []
    rust_times = []

    for i, cid in enumerate(chunk_ids):
        test_logger.info(f"[{i+1}/{len(chunk_ids)}] Chunk {cid}...")
        try:
            r = run_chunk(cid, quiet_logger)
        except Exception as e:
            test_logger.error(f"  EXCEPTION on chunk {cid}: {e}")
            traceback.print_exc()
            r = {"chunk_id": cid, "status": f"error: {e}",
                 "n_haplotypes": 0, "n_qnames": 0,
                 "python_time_s": None, "rust_time_s": None,
                 "speedup": None, "correct_match": None, "mismap_match": None,
                 "correct_only_python": 0, "correct_only_rust": 0,
                 "mismap_only_python": 0, "mismap_only_rust": 0}
            n_error += 1

        results.append(r)

        if r["status"] == "ok":
            n_ok += 1
            py_times.append(r["python_time_s"])
            rust_times.append(r["rust_time_s"])
            if r["correct_match"] and r["mismap_match"]:
                n_match += 1
                test_logger.info(
                    f"  MATCH  haps={r['n_haplotypes']}  qnames={r['n_qnames']}  "
                    f"py={r['python_time_s']:.3f}s  rust={r['rust_time_s']:.3f}s  "
                    f"speedup={r['speedup']:.1f}x"
                )
            else:
                n_mismatch += 1
                test_logger.warning(
                    f"  MISMATCH  chunk={cid}  correct_match={r['correct_match']}  "
                    f"mismap_match={r['mismap_match']}  "
                    f"correct_diff=+{r['correct_only_rust']}/-{r['correct_only_python']}  "
                    f"mismap_diff=+{r['mismap_only_rust']}/-{r['mismap_only_python']}"
                )
                if "diff_detail" in r:
                    for k, v in r["diff_detail"].items():
                        if v:
                            test_logger.warning(f"    {k}: {v[:10]}{'...' if len(v)>10 else ''}")
        elif "skip" in r["status"] or "missing" in r["status"] or "failed" in r["status"]:
            n_skipped += 1
            test_logger.info(f"  SKIPPED: {r['status']}")

    # ── Write CSV ─────────────────────────────────────────────────
    csv_fields = [
        "chunk_id", "status", "n_haplotypes", "n_qnames",
        "python_time_s", "rust_time_s", "speedup",
        "correct_match", "mismap_match",
        "correct_only_python", "correct_only_rust",
        "mismap_only_python", "mismap_only_rust",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    # ── Console summary ───────────────────────────────────────────
    print("\n" + "=" * 70)
    print("DUAL-PATH INTEGRATION TEST SUMMARY")
    print("=" * 70)
    print(f"Chunks tested:    {len(chunk_ids)}")
    print(f"  OK (compared):  {n_ok}")
    print(f"  Skipped:        {n_skipped}")
    print(f"  Errors:         {n_error}")
    print(f"  MATCH:          {n_match}")
    print(f"  MISMATCH:       {n_mismatch}")
    print()

    if py_times:
        py_arr = np.array(py_times)
        rust_arr = np.array(rust_times)
        speedups = py_arr / np.where(rust_arr > 0, rust_arr, 1e-9)

        print(f"Timing (n={len(py_times)} chunks compared):")
        print(f"  Python:  mean={py_arr.mean():.3f}s  median={np.median(py_arr):.3f}s  total={py_arr.sum():.1f}s")
        print(f"  Rust:    mean={rust_arr.mean():.3f}s  median={np.median(rust_arr):.3f}s  total={rust_arr.sum():.1f}s")
        print(f"  Speedup: mean={speedups.mean():.2f}x  median={np.median(speedups):.2f}x  "
              f"min={speedups.min():.2f}x  max={speedups.max():.2f}x")

    print(f"\nResults saved to: {OUTPUT_CSV}")

    if n_mismatch > 0:
        print(f"\nWARNING: {n_mismatch} chunks had mismatches!")
        for r in results:
            if r.get("correct_match") is False or r.get("mismap_match") is False:
                print(f"  Chunk {r['chunk_id']}: correct_match={r['correct_match']} mismap_match={r['mismap_match']}")
        sys.exit(1)
    elif n_match > 0:
        print(f"\nSUCCESS: All {n_match} compared chunks match perfectly.")
    print("=" * 70)


if __name__ == "__main__":
    main()
