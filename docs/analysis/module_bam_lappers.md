# Module: bam_lappers.rs — Technical Reference

## Paired-End BAM Collation Strategy

In coordinate-sorted BAMs, R1 and R2 of the same fragment can be far apart. `samtools collate` groups reads by qname so paired reads are adjacent in the stream. rust-htslib has no native collate — it wraps htslib C bindings, but `samtools collate` is a CLI tool, not an htslib API. Shelling out is the only option.

### 3-Tier Fallback

| Priority | Method | Disk I/O | Memory |
|----------|--------|----------|--------|
| 1st | `spawn_collate_pipe` — pipe stdout via `/dev/fd/N` | None (kernel pipe buffer) | Stream one qname group at a time |
| 2nd | `collate_bam_file` — temp file fallback | Temp BAM written + read | Stream one qname group at a time |
| 3rd | Two-pass in-memory (if samtools unavailable) | None | Entire BAM loaded twice |

### Pipe Mechanics (`spawn_collate_pipe`)

Spawns `samtools collate -f -@ 4 <bam> -o -` with `Stdio::piped()`, takes child stdout fd, opens via `Reader::from_path("/dev/fd/N")`. Uses `std::mem::forget(stdout)` to transfer fd ownership to htslib. Child process reaped after all records consumed.

`/dev/fd/N` is a Linux virtual filesystem entry — opening it gives access to the in-memory pipe, no disk I/O. Data flows: `samtools -> kernel pipe buffer (64KB) -> rust-htslib`.

`Reader::from_stdin()` exists in rust-htslib but only reads fd 0; can't redirect an arbitrary child pipe to it without `dup2`. `/dev/fd/N` is more flexible.

### Per-Read Intervals vs Merged Bounding Boxes

Rust uses per-read intervals instead of Python's merged R1+R2 bounding boxes. Each alignment record gets its own Lapper interval, eliminating false-positive hits in the gap between mate pairs.

### Memory Note

The pipe/collation streaming eliminates the duplicate in-memory `reads_by_qname` HashMap for grouping. However, `read_dict` (all retained Record objects) must remain in memory — the downstream inspection pipeline requires random access to any read at any time.

## Revisions (2026-06-10)

- `fast_median` (used by `is_read_noisy`) now returns `f32` to match Python `np.median`: for even-length quality arrays it no longer integer-truncates the mean of the two middle values (median of `[15,16]` is `15.5`, not `15`). This removes a read-filtering divergence at `median == cutoff + 0.5`, where the old `u8` version dropped reads Python kept. Caller compares `median_q <= basequal_median_filter as f32`.
- Removed dead `cigar_op_to_char` (test-only, no production caller).

## Validation

**Single-end (HG002 chr1:1633000-1635000):** 82/82 qname concordance with samtools, zero missing/extra.

**Paired-end (HG002 full BAM):** `spawn_collate_pipe` via fd 4, 421K reads processed, 196K qname pairs retained, 24 chromosomes, 4.89s (debug build).

**Test logs:**
- `/paedyl01/disk1/yangyxt/test_tmp/test_collate_pipe.log`
