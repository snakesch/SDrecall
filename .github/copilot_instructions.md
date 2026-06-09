# SDrecall Rust Migration Project Instructions

## Project Context

This is the **SDrecall Rust migration worktree** for migrating performance-critical Python code to Rust. The goal is to achieve 2.5-3.7× speedup by replacing Python NCLS and compute-intensive operations with Rust implementations. All migrated code will still need to work under mamba env SDrecall

## Migration Strategy

### Core Principle: One Function at a Time

**CRITICAL:** When working on this migration:

1. **One function per agent** - Each function migration should be handled by a dedicated agent
2. **Keep architecture in mind** - Always consider the total data structure passing architecture
3. **Maintain data flow** - Ensure compatibility between migrated and pending functions
4. **Test incrementally** - Validate each function before moving to the next

### Migration Plan

The detailed migration plan is located at: [docs/analysis/RUST_MIGRATION_PLAN.md](docs/analysis/RUST_MIGRATION_PLAN.md)

**Always refer to this plan when:**
- Starting a new migration task
- Understanding the overall architecture
- Checking what has been completed
- Planning the next steps

### Current Status (Updated 2025-03-05)

✅ **Completed:**
- Worktree created
- `build_phasing_graph` module (already migrated)
- `haplotype_inspection` module structure created (lib.rs, Cargo.toml, build.sh, pyproject.toml)
- Data structures defined (`structs.rs`)
- BAM → Lapper implementation (`bam_lapper.rs` — 499 lines, full filtering logic)
- Test binary for BAM → Lapper (`test_bam_lapper.rs`)
- **BAM Lapper validated:** 100% qname concordance with samtools on HG002 chr1:1633000-1635000 (82/82 qnames match)

🚧 **In Progress (~35% overall):**
- `consensus.rs` — stub only, needs pileup-based consensus assembly
- `similarity.rs` — edit_distance helper done, calculate_similarity and judge_misalignment stubs
- `inspection.rs` — main loop structure outlined, needs implementation
- `python_bindings.rs` — function signature present, data conversion incomplete (7 dicts/lists missing)

⏳ **Pending:**
- Complete consensus assembly implementation
- Complete similarity calculations (`stat_refseq_similarity`)
- Implement main inspection loop
- Complete Python data structure conversion in bindings
- Benchmark BAM → Lapper vs Python NCLS
- Integration testing with Python pipeline
- End-to-end validation on HG006 sample

## Module Structure

### rust_modules/haplotype_inspection/

This is the main module being developed:

```
haplotype_inspection/
├── Cargo.toml              # Dependencies and build config
├── pyproject.toml          # Python packaging config
├── build.sh                # Build script
├── README.md               # Module documentation
└── src/
    ├── lib.rs              # Module entry point
    ├── structs.rs          # Shared data structures
    ├── bam_lapper.rs       # BAM reading + Lapper (replaces NCLS)
    ├── consensus.rs        # Consensus assembly
    ├── similarity.rs       # Reference similarity calculations
    ├── inspection.rs       # Main inspection logic
    └── python_bindings.rs  # Python interface
```

### Key Files to Reference

- **Python source:** `fp_control/identify_misaligned_haps.py`
  - `assemble_consensus` (lines 288-341)
  - `judge_misalignment_by_extreme_vardensity` (lines 369-421)
  - `stat_refseq_similarity` (lines 888-1018)
  - `inspect_by_haplotypes` (lines 1180-1254)

- **Existing Rust module:** `rust_modules/build_phasing_graph/`
  - Reference for code style and patterns
  - Shows how to structure PyO3 bindings

## Development Workflow

### When Migrating a Function

1. **Read the Python source** - Understand the function thoroughly
2. **Identify dependencies** - What data structures and helper functions are needed?
3. **Design Rust equivalent** - Consider performance and safety
4. **Implement incrementally** - Start with core logic, add optimizations later
5. **Add tests** - Unit tests for correctness
6. **Update migration plan** - Check off completed tasks

### Standalone Compilation Principle

**CRITICAL:** When testing ported functions, compile ONLY the completed modules. Do NOT include stub/incomplete modules in the compilation. Use `#[path]` includes to bypass `lib.rs`:

```rust
// In test binary — only pull in the modules you need:
#[path = "../structs.rs"]
mod structs;
#[path = "../bam_lapper.rs"]
mod bam_lapper;
// Do NOT include consensus.rs, similarity.rs, etc. if they are stubs
```

This avoids compilation errors from incomplete modules and ensures test isolation.

### Build Environment (SDrecall conda env)

```bash
# Proper conda activation (NEVER use `source activate`):
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall

# Required environment variables for compilation:
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

# Build a specific binary:
cd rust_modules && cargo build --bin test_bam_lapper

# Build the cdylib for Python:
cd rust_modules/haplotype_inspection && maturin build --release
pip install target/wheels/haplotype_inspection-*.whl
```

**Why these env vars:**
- `LIBCLANG_PATH`: rust-htslib needs libclang for bindgen
- `OPENSSL_NO_VENDOR` + `PKG_CONFIG_PATH`: Prevents curl-sys from building curl from source with incompatible OpenSSL 1.1 API. Uses conda's system curl (built against OpenSSL 3.x) instead.

### Building and Testing

```bash
cd rust_modules/haplotype_inspection
./build.sh

# Or manually:
maturin build --release
pip install target/wheels/haplotype_inspection-*.whl
```

### Integration Pattern

The Rust module will be called from Python with a feature flag:

```python
USE_RUST_HAPLOTYPE_INSPECTION = True

if USE_RUST_HAPLOTYPE_INSPECTION:
    from haplotype_inspection import inspect_haplotypes_rust
    correct_qnames, mismap_qnames = inspect_haplotypes_rust(...)
else:
    # Fallback to Python implementation
    correct_qnames, mismap_qnames = inspect_by_haplotypes(...)
```

## Performance Targets

| Component | Current (Python) | Target (Rust) | Speedup |
|-----------|------------------|---------------|---------|
| BAM → NCLS/Lapper | 19s | 5-8s | 2.4-3.8× |
| inspect_by_haplotypes | 60s | 10-15s | 4-6× |
| **Total per subprocess** | **82s** | **22-30s** | **2.7-3.7×** |

## Important Reminders

### Data Structure Architecture

- **Input from Python:** Phasing results (`qname_hap_info`, `hap_qname_info`)
- **Rust processing:** BAM reading, interval queries, consensus, similarity
- **Output to Python:** Correct and misaligned read names

### Code Quality

- Follow existing Rust patterns in `build_phasing_graph`
- Use `rustc-hash::FxHashMap` for integer keys
- Use `ahash` for string keys
- Add logging with `log` crate (forwarded to Python via `pyo3-log`)
- Write unit tests for all functions

### Safety and Correctness

- Validate against Python implementation
- Keep Python fallback for safety
- Extensive error handling
- Clear error messages

## Next Steps

Refer to [docs/analysis/RUST_MIGRATION_PLAN.md](docs/analysis/RUST_MIGRATION_PLAN.md) for:
- Detailed implementation roadmap
- Function-by-function migration tasks
- Testing strategy
- Success criteria

## Questions?

When in doubt:
1. Check the migration plan
2. Look at `build_phasing_graph` for patterns
3. Read the Python source code
4. Ask for clarification before making assumptions

---

**Remember:** This is a performance-critical migration. Take time to understand the code, maintain correctness, and optimize incrementally.
