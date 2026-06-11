"""
Differential-validation dump hook for the Rust migration (T1 + T3).

This module is inert unless the environment variable ``SDRECALL_DIFF_DUMP_DIR``
points at a writable directory. When set, ``realign_filter_per_cov`` calls into
here once per coverage island to:

  * T3 (phasing): dump the weight matrix + node_read_ids + edge_weight_cutoff and
    the Python ``phasing_realigned_reads`` partition, so the Rust ``phasing`` crate
    can be run on identical inputs and compared up to relabeling.
  * T1 (haplotype inspection): run the *Python* ``inspect_by_haplotypes`` baseline
    on the same inputs the Rust path receives, then compare ``correct``/``mismap``
    set-equality island by island.

Every entry point is wrapped so a validation failure NEVER crashes the pipeline
(it logs and returns). Each island writes to its own sub-directory keyed by the
per-island BAM basename, so the multiprocessing workers never contend on a file.
"""

import os
import glob
import json

import numpy as np


def _dump_dir():
    """Return the dump root if validation is enabled, else None."""
    d = os.environ.get("SDRECALL_DIFF_DUMP_DIR")
    if not d:
        return None
    os.makedirs(d, exist_ok=True)
    return d


def _island_dir(root, bam):
    name = os.path.basename(bam)
    if name.endswith(".bam"):
        name = name[:-4]
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    return d


def _write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh)


def _node_read_ids_to_map(node_read_ids):
    """node_read_ids is a list indexed by vertex; each entry is a (read_id1, read_id2_or_None)
    tuple of string read ids. Normalise to {vertex_str: [read_id, ...]} dropping None."""
    out = {}
    if hasattr(node_read_ids, "items"):  # defensive: dict form
        items = node_read_ids.items()
    else:
        items = enumerate(node_read_ids)
    for vertex, entry in items:
        if entry is None:
            rids = []
        elif isinstance(entry, (list, tuple)):
            rids = [str(r) for r in entry if r is not None]
        else:
            rids = [str(entry)]
        out[str(int(vertex))] = rids
    return out


def dump_phasing(bam,
                 phased_graph,
                 weight_matrix,
                 edge_weight_cutoff,
                 node_read_ids,
                 qname_to_node,
                 qname_hap_info,
                 hap_qname_info,
                 total_readhap_vector,
                 total_readerr_vector,
                 mean_read_length,
                 recall_mq_cutoff,
                 basequal_median_cutoff,
                 intrinsic_bam,
                 logger):
    """T3 dump: phasing inputs + Python partition. No-op unless enabled."""
    root = _dump_dir()
    if root is None:
        return
    try:
        d = _island_dir(root, bam)
        np.save(os.path.join(d, "weight_matrix.npy"),
                np.ascontiguousarray(weight_matrix, dtype=np.float32))

        # node_read_ids: vertex_idx -> [read_id_str, ...]  (string read ids "qname:flag")
        _write_json(os.path.join(d, "node_read_ids.json"),
                    _node_read_ids_to_map(node_read_ids))
        # per-read hap / error vectors keyed by string read id (round-2 variant check needs them)
        _write_json(os.path.join(d, "read_hap.json"),
                    {str(k): [int(x) for x in np.asarray(v).tolist()]
                     for k, v in total_readhap_vector.items()})
        _write_json(os.path.join(d, "read_err.json"),
                    {str(k): [float(x) for x in np.asarray(v).tolist()]
                     for k, v in total_readerr_vector.items()})
        # qname_to_node: qname -> vertex_idx
        _write_json(os.path.join(d, "qname_to_node.json"),
                    {str(k): int(v) for k, v in qname_to_node.items()})
        # vertex_idx -> qname (from the graph vertex property), to compare qname partitions
        vertex_qname = {int(v): phased_graph.vp.qname[v] for v in phased_graph.vertices()}
        _write_json(os.path.join(d, "vertex_qname.json"),
                    {str(k): v for k, v in vertex_qname.items()})

        # Exact graph adjacency (the edge set gt.label_components uses). Components must be
        # reconstructed from THIS, not from the weight matrix: positive-weight overlaps define
        # edges, while 0 (no overlap) and -1 (incompatible) entries are non-edges.
        edges = [[int(e.source()), int(e.target())] for e in phased_graph.edges()]
        _write_json(os.path.join(d, "edges.json"), edges)

        # Python phasing partition outputs
        _write_json(os.path.join(d, "phasing_qname_hap_info.json"),
                    {str(k): int(v) for k, v in qname_hap_info.items()})
        _write_json(os.path.join(d, "phasing_hap_qname_info.json"),
                    {str(k): sorted(v) for k, v in hap_qname_info.items()})

        _write_json(os.path.join(d, "meta.json"), {
            "bam": os.path.abspath(bam),
            "intrinsic_bam": os.path.abspath(intrinsic_bam),
            "edge_weight_cutoff": float(edge_weight_cutoff),
            "mean_read_length": float(mean_read_length),
            "recall_mq_cutoff": int(recall_mq_cutoff),
            "basequal_median_cutoff": int(basequal_median_cutoff),
            "num_vertices": int(phased_graph.num_vertices()),
            "num_edges": int(phased_graph.num_edges()),
            "weight_matrix_shape": list(np.asarray(weight_matrix).shape),
            "n_haplotypes": len(hap_qname_info),
        })
        logger.info(f"[diff_dump] T3 phasing inputs dumped to {d}")
    except Exception as e:
        logger.warning(f"[diff_dump] dump_phasing failed (non-fatal): {e}")


def compare_inspection(bam,
                       intrinsic_bam,
                       hap_qname_info,
                       qname_hap_info,
                       qname_to_node,
                       node_read_ids,
                       read_id_read_dict,
                       total_readhap_vector,
                       total_readerr_vector,
                       read_ref_pos_dict,
                       total_lowqual_qnames,
                       compare_haplotype_meta_tab,
                       mean_read_length,
                       recall_mq_cutoff,
                       basequal_median_cutoff,
                       rust_correct,
                       rust_mismap,
                       logger):
    """T1 differential: run Python inspect_by_haplotypes baseline and compare to
    the Rust result already computed by the caller. No-op unless enabled.

    The Python baseline IS the slow ~50s/island path being replaced, so cap the number of
    comparisons via SDRECALL_DIFF_MAX_INSPECT (default 20; set <0 for unlimited). The cap is
    enforced by counting already-written inspect_diff.json files across all workers (a small
    race may let it overshoot slightly, which is fine for a sampled validation)."""
    root = _dump_dir()
    if root is None:
        return
    try:
        cap = int(os.environ.get("SDRECALL_DIFF_MAX_INSPECT", "20"))
    except ValueError:
        cap = 20
    if cap >= 0 and len(glob.glob(os.path.join(root, "*", "inspect_diff.json"))) >= cap:
        return
    try:
        from fp_control.bam_ncls import migrate_bam_to_ncls
        from fp_control.identify_misaligned_haps import inspect_by_haplotypes

        # Rebuild the NCLS structures exactly as the legacy Python path did
        # (realign_filter_per_cov.py lines 222 + 235-244).
        bam_ncls = migrate_bam_to_ncls(bam,
                                       mapq_filter=recall_mq_cutoff,
                                       basequal_median_filter=basequal_median_cutoff,
                                       logger=logger)
        if bam_ncls is None:
            logger.warning("[diff_dump] bam_ncls is None; skipping inspect comparison")
            return
        intrin_bam_ncls = migrate_bam_to_ncls(intrinsic_bam,
                                              mapq_filter=0,
                                              basequal_median_filter=0,
                                              paired=False,
                                              filter_noisy=False,
                                              logger=logger)
        intrin_bam_ncls = intrin_bam_ncls[:-1]

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
            {},  # total_genomic_haps (empty, as in the production path)
            read_ref_pos_dict,
            compare_haplotype_meta_tab=compare_haplotype_meta_tab,
            mean_read_length=mean_read_length,
            logger=logger,
        )
        py_correct = set(py_correct)
        py_mismap = set(py_mismap)
        rust_correct = set(rust_correct)
        rust_mismap = set(rust_mismap)

        correct_equal = (py_correct == rust_correct)
        mismap_equal = (py_mismap == rust_mismap)

        d = _island_dir(root, bam)
        _write_json(os.path.join(d, "inspect_diff.json"), {
            "n_haplotypes": len(hap_qname_info),
            "py_correct_n": len(py_correct),
            "py_mismap_n": len(py_mismap),
            "rust_correct_n": len(rust_correct),
            "rust_mismap_n": len(rust_mismap),
            "correct_equal": correct_equal,
            "mismap_equal": mismap_equal,
            "correct_only_py": sorted(py_correct - rust_correct),
            "correct_only_rust": sorted(rust_correct - py_correct),
            "mismap_only_py": sorted(py_mismap - rust_mismap),
            "mismap_only_rust": sorted(rust_mismap - py_mismap),
        })
        status = "MATCH" if (correct_equal and mismap_equal) else "MISMATCH"
        logger.info(
            f"[diff_dump] T1 inspect {status} for {os.path.basename(bam)}: "
            f"correct_equal={correct_equal} mismap_equal={mismap_equal} "
            f"(py {len(py_correct)}/{len(py_mismap)} vs rust {len(rust_correct)}/{len(rust_mismap)})"
        )
    except Exception as e:
        import traceback
        logger.warning(f"[diff_dump] compare_inspection failed (non-fatal): {e}\n{traceback.format_exc()}")
