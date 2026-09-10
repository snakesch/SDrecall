# Vector Encoding & Consensus — Technical Reference

Covers `pairwise_read_inspection.rs` (vector extraction, counting) and consensus/batch operations in `identify_misaligned_haps.rs`.

## Haplotype Vector Encoding (i16)

| Value | Meaning |
|-------|---------|
| `1` | Match/Reference (CIGAR `=`) |
| `-4` | SNV (mismatch `X`, excluding N bases) |
| `-6` | Deletion (CIGAR `D`) |
| `>1` | Insertion (value = length * 4, placed at first position of NEXT ref-consuming op) |
| `-10` | Padding (unused positions in fixed-size arrays) |

- Insertion at start of alignment (hap vector empty) is ignored
- CIGAR `M` op (op 0) is rejected — requires `=/X` mode (`--eqx`)

> **This table documents the CURRENT (Python-faithful) encoding** — keep it through T1 differential validation. The TARGET encoding for the consolidated pure-Rust function is the golden scheme below; do not switch until after T1 (see caveat).

## Golden encoding — decided 2026-06-11 (DivA + DivB), TARGET for consolidated `sdrecall-utils`

**DivA — reject `M`, with an error not a panic.** Only `=`/`X`/`I`/`D`/`N`/`S` are valid. An `M` op means the read was not aligned in `--eqx` mode, which never happens in this pipeline (all reads are minimap2-realigned). The strict (`haplotype_inspection`) behavior — reject `M` — is correct and matches Python's `assert operation != 0`; the lenient `build_phasing_graph` copy (treat `M` as match) is wrong and must not be the one we keep. **Mechanism:** the current `panic!` becomes a typed `Result` error that early-stops and propagates to the job boundary — fail-fast with a clear message, never a panic or a silent mis-encode.

**DivB — summation insertion encoding (fixes the overwrite information loss).**

| Value | Meaning (golden) |
|-------|------------------|
| `1` | Match (`=`) |
| `-4` | SNV (`X`) |
| `-10` | Deletion (`D`) — per deleted ref position (was `-6`) |
| `> 1` | Insertion present: `base + 10*L` where L = inserted length, **added** to the prepending ref-consuming base (was `4*L`, written by overwriting) |

- **Decode:** `value > 1` ⟹ insertion; `structural = round_to_nearest_10(value)` (= `10*L`); `delta = value − structural ∈ {+1 = match, −4 = mismatch, 0 = pure deletion}`. Examples: match + 2bp ins = `1 + 20 = 21` → (struct 20, delta +1); mismatch + 2bp ins = `−4 + 20 = 16` → (struct 20, delta −4).
- **Why `-10` not `-6`:** deletion positions must round to a clean multiple of 10 (delta 0). `−6` would round to −10 with delta +4, which is not a valid point signal and would decode ambiguously.
- **Improvement over current:** the present scheme overwrites the base under/adjacent to an insertion (the `−4`/`1` is lost); summation makes it recoverable via `delta`. Mathematically unambiguous because the only point signals are `+1` and `−4`, both with |delta| < 5 < (half the indel unit).
- **Required invariants (assert):** no adjacent indels — the aligner emits consecutive mismatches instead, so del-then-ins never co-occupy a base; and no insertion at vector start (nothing to attach to). If either is violated the decode is ambiguous.
- **⚠ Padding-sentinel collision (must fix in the same change):** the current hap-vector **padding** value is `-10` (filtered everywhere by `>= -8` = "is real"; see batch extraction + `assemble_consensus`). Golden deletion is also `-10` → they collide, and deletions would be silently dropped as padding. **Decided fix (2026-06-11): move the hap padding sentinel to `-20`** — `const HAP_PAD: i16 = -20;` — which sits one indel-unit below the most-negative real signal (deletion `-10`), and change the real-value filter from `>= -8` to `!= HAP_PAD` (equiv. `> -20`). `-20` can never be a real value (insertion compounds are positive; the only way to reach `-20` would be two stacked deletions, excluded by the no-adjacent-indel invariant). Sites to update: the `Array2::from_elem(..., -10)` fill in `record_hap_err_vectors_per_region`, the `>= -8` filter in `assemble_consensus`, and the `-10`/`-9` padding tests. The **err** vector's `-10.0` padding does not collide (golden only changes the i16 hap encoding; `assemble_consensus` truncates qual by the seq-derived count, never testing the qual pad) — leave it.

### Golden variant counting — compound signals (decided 2026-06-11)

Under summation, one base can carry a **compound** signal (e.g. SNV *and* a following insertion). The counters must not miss the hidden point signal:

- **`count_snv` (golden):** a mismatch under an insertion is `-4 + 10*L`, which always ends in digit **6** (`6, 16, 26, …`); a match under an insertion is `1 + 10*L`, ending in **1**. So:
  ```rust
  // keep the original vectorized pass (pure SNV) ...
  let pure   = a.iter().filter(|&&v| v == -4).count();
  // ... and ADD a second vectorized pass for compound SNV (mismatch + insertion):
  let compound = a.iter().filter(|&&v| v > 1 && v % 10 == 6).count();
  let snv = (pure + compound) as i32;
  ```
- **`count_continuous_indel_blocks` (golden):** deletion is now `-10` (was `-6`); insertions are still detected by `v > 1` (this already catches both `…1` match-ins and `…6` mismatch-ins). So the predicate becomes `v == -10 || v > 1`.
- **Consequence (intended):** a compound base counts as **both** an SNV and an indel block — which is correct (it really is two variant events). The current overwrite encoding loses the SNV there, so today's `count_var` *undercounts* compounds. Golden fixes that, and it is exactly why golden agrees with the CIGAR oracle below while the current encoding does not.
- **No loophole** other than del-adjacent-to-ins (`-10 + 10*L`, ends in `0`), which the aligner prevents (emits mismatches) and which we assert against.

**⚠ Validation caveat — why we don't switch yet.** The golden encoding is numerically different from Python's, so it cannot be checked by Rust-vs-Python hap-vector parity. Plan: (1) keep the Python-faithful encoding above through **T1** (vector-level + output-level parity); (2) switch to golden in **T4** / once Python is retired, re-validating only at the final `correct`/`mismap` set level. The switch is **all-or-nothing**: `assemble_consensus`, `count_var`/`count_snv`/`count_continuous_indel_blocks`, `numba_shared_variant_positions`, and `ref_genome_similarity` must all decode the new scheme in the same change.

## Golden-encoding validation — CIGAR/pileup oracle (decided 2026-06-11)

Because golden diverges from Python, **Python is retired as the oracle for these encoding-dependent functions** at the switch. The replacement ground truth is the raw alignment itself (CIGAR + pileup against the reference) — which is also **encoding-agnostic**, so the same harness validates the current encoding *now* and golden *after the switch*. Build it during **T1** (over a committed real-island BAM slice + reference); it becomes the **golden-switch safety net in T4** unchanged.

- **Variant counts (`count_snv`/`count_continuous_indel_blocks`/`count_var`):** oracle walks each read's CIGAR — `#SNV = Σ len(X ops)`; `#indel blocks = number of maximal runs of consecutive I/D ops` (robust to any adjacency). Assert against the counters over every read. Note: this oracle expresses the *true* count, so it agrees with **golden** and will (correctly) disagree with the **current/Python** encoding precisely on compound SNV+ins bases — that gap is the documented bug golden fixes, not a harness error.
- **`assemble_consensus`:** oracle is a `rust-htslib` pileup over the member reads + reference FASTA — at each ref position, encode each covering read's base from scratch (`base == ref` ⇒ `1`, else `-4`; indels from pileup flags), take the lowest-error value ≤ 0.2 (default `1`), and assert equality with the consensus over the overlap span. Pair with a brute-force per-position recompute over the same encoded inputs to separately catch vectorization bugs.
- **`numba_shared_variant_positions`:** build each side's position→variant-type map from raw data (homolog refseq from its CIGAR; consensus from the pileup machinery above), intersect by type, compare. Encoding-free on both sides; insertion-length match made explicit rather than relying on integer equality.
- **`ref_genome_similarity`:** validated transitively by the count oracle plus one integration assert that it equals `(var_from_cigar(query), snv_from_cigar(genomic), indel_from_cigar(genomic))`, with a targeted case for the all-`1` short-circuit.

**Through T1 (current encoding):** keep Python parity as the gate; the CIGAR/pileup oracle is built and used for the encoding-agnostic checks (consensus via pileup, base tallies) and its compound-count divergence is recorded as the quantified motivation for the switch. **At T4 (golden):** the CIGAR/pileup oracle becomes the gate; Python parity is dropped for these functions.

## Error Vector Encoding (f32)

| Value | Meaning |
|-------|---------|
| `0.0-1.0` | Error probability (Phred: `10^(-Q/10)`) |
| `0.0` | Placeholder for insertions/deletions |
| `-10.0` | Padding |

- Insertion sets `err_vector[ref_pos - 1] = 0.0` (BEFORE the insertion) — differs from hap_vector which places marker AFTER
- `phred_to_prob` reads from a precomputed 256-entry table (`PHRED_TO_PROB`, built once via `LazyLock`) instead of calling `powf` per base; values are bit-identical (2026-06-10)

## Batch Extraction (`record_hap_err_vectors_per_region`)

Returns 2D padded arrays matching Python's `np.full((n_reads, max_len), -10)` convention:
- `read_spans: Array2<i32>` shape `(N, 2)` — `[ref_start, ref_end]` per read
- `hap_vectors: Array2<i16>` shape `(N, max_len)` — padded with `-10`
- `err_vectors: Array2<f32>` shape `(N, max_len)` — padded with `-10.0`

Caches keyed by read ID (`"qname:flag"`) to avoid recomputation across regions. Cached vectors are copied straight into the padded `Array2` rows — no intermediate per-read `Array1` clone (2026-06-10).

## Consensus Assembly (`assemble_consensus`)

1. Find global start/end from all read spans
2. Initialize consensus arrays (seq: i16, qual: f32)
3. For each read: filter padding (< -8), map to consensus coordinates, update where current read has better quality (lower error prob <= 0.2)

Quality semantics: smaller error probability = better quality. Position updated only when read's error prob <= current AND <= 0.2 threshold.

## Variant Density Filtering (`judge_misalignment_by_extreme_vardensity`)

Three-threshold misalignment detection on consensus:
1. density > 5/85 with padding=42, requires >=1 indel block
2. density > 6/131 with padding=65, requires >=1 indel block
3. density > 10/148 with padding=74, no indel requirement

`count_window_var_density` builds two prefix-sum arrays once (SNV counts + indel-block starts) and answers each sliding window in O(1) — O(n) overall, no per-position allocation — replacing the former O(n·window) re-slice-and-recount. Output is identical, cross-checked against the brute-force reference in `test_count_window_var_density_matches_reference`; ~39× faster on a 2000-base microbench (`bench_count_window_var_density`, `#[ignore]`'d). The helper `count_continuous_blocks` now slices `s![1..arr.len()+1]` (was `s![1..-1]`, valid ndarray but tripped clippy's deny-by-default `reversed_empty_ranges`). (2026-06-10)

## Validation

**Consensus module (HG002 chr1:1633000-1635000):** 140/140 reads pass per-read CIGAR validation; 14 HP tag groups tested; all consensus length checks pass.

**Cross-validation helpers:** `rank_unique_values` and `calculate_coefficient` cross-validated against exact Python/numba output using `/paedyl01/disk1/yangyxt/test_tmp/cross_validate_helpers.py` (12 tests with complex realistic inputs, 6 per function).

**Test logs:**
- `/paedyl01/disk1/yangyxt/test_tmp/test_consensus_debug.log`
- `/paedyl01/disk1/yangyxt/test_tmp/test_helpers_crossval_debug.log`
