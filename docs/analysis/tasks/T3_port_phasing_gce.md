# T3 — Port phasing + GCE → `phasing` crate

**Crate:** `phasing` (new — lib + bin)
**Status (2026-07-17):** Production-integrated. Graph build, sparse phasing/GCE, and HP writing run in Rust; 58 focused phasing tests passed in the production verification set.
**Depends on:** T0, T2
**Track:** Historical June track B; complete.
**Replaces (Python):** `fp_control/phasing.py` (`phasing_realigned_reads`), `fp_control/gce_algorithm.py` (Greedy Clique Expansion), and the GCE-related `@njit` kernels in `fp_control/numba_operators.py`.

## Goal & scope boundary

Produce the two phasing outputs **purely in Rust**:
- `qname_hap_info`: `HM<vertex_idx i32 → hap_id i32>` (note: keyed by vertex index despite the name — confirmed at `phasing.py:224,260`)
- `hap_qname_info`: `HM<hap_id i32 → {qname}>`

from `(weight_matrix, edge_weight_cutoff, hap/err vectors, node_read_ids)` + the graph adjacency. This removes the **last Python step in Phase 2c** and the **last numba dependency** in the hot path.

Out of scope: the BAM read and graph build (that's `phasing-graph`); the consensus/similarity inspection (that's `haplotype-inspection`). Those join in T4.

## Data flow

```
weight_matrix (NxN f32, -1 = incompatible pair, 0..1 = edge weight)
        │
        ├─ connected components            (replaces gt.label_components)  → union-find / petgraph
        │
        └─ per component:
              GCE clique expansion          (gce_algorithm.py:266-545)
                seed = global-max-weight vertex
                greedy expand by highest-weight in-mask neighbor
                9-member lookback + cutoff backtrack
              → clique = set of original-graph vertex indices
        │
        ▼
   haplotype assignment → qname_hap_info, hap_qname_info
```

GCE operates on a **CSR sparse** view of the submatrix (Python uses `scipy.sparse`; Rust uses `sprs` or a hand-rolled CSR). All ~22 numba kernels are **pure array/CSR math** (one, `numba_isin`, uses a set) — mechanical ports.

## Dependencies

- Crates: `petgraph` 0.8.3 (components / subgraph filtering), `sprs` 0.11.4 (CSR sparse for GCE), `ndarray` 0.15, `rayon` 1.8 (parallelize across components), `sdrecall-utils`.
- Foundation: `sdrecall-utils`, `sdrecall-io`.
- External tools: none.

## Performance bottleneck / rationale

Historical baseline: the Python/Numba GCE path took about 5 s/island and ran its `para_*` kernels with `parallel=False`. Production now uses the Rust CSR/GCE implementation under the shared phase-aware resource model; additional parallelism is performance work and must be benchmarked.

**Risk:** the only graph-tool coupling (`label_components`, `GraphView` vertex-filter, clique-component split, `qname` vertex-property) lives in `phasing.py`, not GCE. All of it is mechanical to replace (union-find + BFS + index bookkeeping), and once fused in T4 the graph never leaves Rust, so no graph-tool object is needed at all.

## Tests

### Unit (tier 1)
- GCE clique output on hand-built small weight matrices with known maximal cliques.
- Each ported numba kernel vs its known output (reuse the `rank_unique_values` / `calculate_coefficient` cross-validation pattern already proven).
- Connected components vs `gt.label_components` on fixture graphs.

### Differential vs Python (tier 2)
- Add a dump hook in `realign_filter_per_cov.py` (just before `phasing_realigned_reads`, ~line 289) to record the `weight_matrix`, `node_read_ids`, and `edge_weight_cutoff` per island, plus the Python `qname_hap_info`/`hap_qname_info`.
- Run the Rust crate on identical matrices; compare.

**Pass criterion:** identical **haplotype partition** of qnames on 100 % of recorded islands. Compare as a **set partition** (hap_id labels are arbitrary, so compare partitions up to relabeling — not raw ids). Any partition mismatch is a hard fail to root-cause.

**Data:** recorded weight matrices from HG002 + HG006 islands.

## Progress
- [x] Scaffold `phasing` crate (lib + bin: matrix in → maps out) — `rust_modules/phasing/`, added to workspace members
- [x] Port numba kernels (unit-tested individually) — `src/kernels.rs` (hand-rolled CSR + all GCE kernels)
- [x] Connected components (union-find) vs Python — `src/phasing.rs::connected_components` (replaces `gt.label_components`)
- [x] GCE clique expansion (CSR) vs Python cliques — `src/gce.rs` (9-member lookback + cutoff backtrack ported 1:1)
- [x] Phasing assembly → two maps — `src/phasing.rs::phase` + `qname_partition`
- [x] Per-island dump harness + partition-equality differential — `src/main.rs` + `fp_control/diff_dump.py::dump_phasing`
- [x] Replace the June "rayon across components" placeholder with the production resource-managed graph/GCE execution model. Further parallelism changes are performance-only and require controlled benchmarks plus semantic parity.

### Differential result (2026-06-11) — **PASS, 270/270 islands**

Ran the full HG002 CMRG example pipeline (`example/example_run.sh`, hg38) with `SDRECALL_DIFF_DUMP_DIR` set, dumping per-island `weight_matrix.npy` + `edges.json` + `node_read_ids.json` + `read_hap.json` / `read_err.json` + the Python `phasing_realigned_reads` partition. The Rust crate (`phasing --path <dump>`) reproduced the **identical qname partition up to relabeling on all 270 islands**.

- Harness log: `/paedyl01/disk1/yangyxt/test_tmp/phasing_differential_20260611.log` (`Result: 270/270 islands match`).
- Unit tests: 15 pass, clippy clean (`/paedyl01/disk1/yangyxt/test_tmp/phasing_unit_tests_20260611.log`).
- Key fidelity points that held on real data: scipy-CSR zero-dropping (keeps `-1`), `>=` argmax tie-break, the round-2 no-variant gate (needs `read_hap`/`read_err`), and clique-splitting via original-edge + weak-band (`0.1 < w ≤ cutoff`) connectivity.
- `node_read_ids` is a **list indexed by vertex** of `(read_id1, read_id2_or_None)` **string** ids; the read hap/err vectors are dicts keyed by those string ids. (Documented here because it tripped the first dump attempt.)

Remaining formal record: HG006 (needs the external BAM fixture). Further component-level parallelism is optional benchmark-driven optimization, not sign-off work.
