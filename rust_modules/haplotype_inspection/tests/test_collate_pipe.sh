#!/bin/bash
# Test: validate pipe-based collation on real HG002 BAM data.
# Runs build_lapper_from_bam in paired mode and checks results.
#
# Usage: bash tests/test_collate_pipe.sh

set -euo pipefail

# Hardcoded test BAM
BAM="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.clean.bam"
LOG="/paedyl01/disk1/yangyxt/test_tmp/test_collate_pipe.log"

if [ ! -f "$BAM" ]; then
    echo "ERROR: Test BAM not found: $BAM"
    exit 1
fi

echo "=== Collate Pipe Test ==="
echo "BAM: $BAM"
echo "Log: $LOG"
echo

# Build and run the test
cd "$(dirname "$0")/.."

# Activate environment
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

echo "Building test binary..."
cargo build --bin test_bam_lapper 2>&1 | tail -5

echo
echo "Running paired-mode test via Rust integration test..."
RUST_LOG=info cargo test --lib test_collate_pipe_real_bam -- --nocapture 2>&1 | tee "$LOG"

echo
echo "Log saved to: $LOG"
