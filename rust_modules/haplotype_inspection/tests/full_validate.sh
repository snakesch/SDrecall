#!/usr/bin/env bash
# =============================================================================
# Full automated comparison: BEDOPS pipeline vs Rust select_regions
#
# Runs both implementations on identical inputs and diffs the outputs.
# Usage: bash tests/full_validate.sh
# Requires: BEDOPS tools (sort-bed, bedops, bedmap) + Rust binary built
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
WORKSPACE_DIR="$(cd "$PROJECT_DIR/../.." && pwd)"
TMPDIR=$(mktemp -d /tmp/validate_full_XXXXXX)
trap "rm -rf $TMPDIR" EXIT

PASS=0
FAIL=0

# ── BEDOPS pipeline (replicates Python's CLI branch) ─────────────────────────
run_bedops() {
    local min_k="$1"; shift
    local output_bed="$1"; shift

    local tagged="$TMPDIR/_tagged.bed"
    local parts="$TMPDIR/_parts.bed"
    local hid=1
    > "$tagged"
    for bed in "$@"; do
        mawk -v id="$hid" 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,id}' "$bed" >> "$tagged"
        hid=$((hid+1))
    done
    sort-bed "$tagged" > "$tagged.sorted" && mv "$tagged.sorted" "$tagged"

    bedops --partition "$tagged" > "$parts"

    bedmap --ec --delim $'\t' --echo --echo-map-id-uniq "$parts" "$tagged" | \
        mawk -F '\t' -v k="$min_k" \
        '{ ids=$NF; n=(ids=="N/A")?0:split(ids,a,";"); if (n>=k) print $1"\t"$2"\t"$3 }' | \
        sort-bed - | bedops --merge - > "$output_bed"
}

# ── Extract Rust output for a given test name ────────────────────────────────
RUST_OUTPUT="$TMPDIR/rust_all.txt"
# Run the Rust binary — we assume it was already built
# target/ lives under rust_modules/ (the workspace root)
RUST_BIN="$PROJECT_DIR/../target/debug/examples/validate_select_regions"
if [[ ! -x "$RUST_BIN" ]]; then
    echo "ERROR: Rust binary not found at $RUST_BIN"
    echo "Build first: cd rust_modules && cargo build --example validate_select_regions"
    exit 1
fi
echo "Using Rust binary: $RUST_BIN"
"$RUST_BIN" > "$RUST_OUTPUT"

extract_rust() {
    local name="$1"
    local outfile="$2"
    # Extract lines between "### <name>" and "---"
    sed -n "/^### $name/,/^---/p" "$RUST_OUTPUT" | grep -v '^###\|^---' > "$outfile" || true
}

compare() {
    local name="$1"
    local bedops_file="$2"
    local rust_file="$3"

    if diff -q "$bedops_file" "$rust_file" > /dev/null 2>&1; then
        local n=$(wc -l < "$bedops_file")
        echo "  PASS: $name ($n regions)"
        PASS=$((PASS+1))
    else
        echo "  FAIL: $name"
        echo "    BEDOPS:"
        cat "$bedops_file" | sed 's/^/      /'
        echo "    Rust:"
        cat "$rust_file" | sed 's/^/      /'
        diff "$bedops_file" "$rust_file" | sed 's/^/      /' || true
        FAIL=$((FAIL+1))
    fi
}

echo "============================================================"
echo " Full BEDOPS vs Rust Comparison"
echo "============================================================"
echo ""

# ── Helper to write test BED files and run both pipelines ────────────────────
run_case() {
    # Args: test_name min_k bed_contents...
    # bed_contents are pairs: filename content
    local test_name="$1"; shift
    local min_k="$1"; shift

    local beds=()
    while [[ $# -ge 2 ]]; do
        local fn="$TMPDIR/$1"
        echo -e "$2" > "$fn"
        beds+=("$fn")
        shift 2
    done

    run_bedops "$min_k" "$TMPDIR/bedops_${test_name}.bed" "${beds[@]}"
    extract_rust "$test_name" "$TMPDIR/rust_${test_name}.bed"
    compare "$test_name" "$TMPDIR/bedops_${test_name}.bed" "$TMPDIR/rust_${test_name}.bed"
}

# Test 1
run_case "Test1_simple_2hap (min=2)" 2 \
    h1.bed "chr1\t100\t300" \
    h2.bed "chr1\t200\t400"

# Test 2
run_case "Test2_staircase_min2 (min=2)" 2 \
    h1.bed "chr1\t100\t400" \
    h2.bed "chr1\t200\t500" \
    h3.bed "chr1\t300\t600"

# Test 3
run_case "Test3_staircase_min3 (min=3)" 3 \
    h1.bed "chr1\t100\t400" \
    h2.bed "chr1\t200\t500" \
    h3.bed "chr1\t300\t600"

# Test 4
run_case "Test4_no_overlap (min=2)" 2 \
    h1.bed "chr1\t100\t200" \
    h2.bed "chr1\t300\t400"

# Test 5
run_case "Test5_identical_3hap (min=2)" 2 \
    h1.bed "chr1\t100\t500" \
    h2.bed "chr1\t100\t500" \
    h3.bed "chr1\t100\t500"

# Test 6
run_case "Test6_contained (min=2)" 2 \
    h1.bed "chr1\t100\t600" \
    h2.bed "chr1\t200\t400"

# Test 7
run_case "Test7_disjoint_zones (min=2)" 2 \
    h1.bed "chr1\t100\t200\nchr1\t500\t600" \
    h2.bed "chr1\t150\t250\nchr1\t550\t650"

# Test 8
run_case "Test8_multi_chrom (min=2)" 2 \
    h1.bed "chr1\t100\t300\nchr2\t500\t700" \
    h2.bed "chr1\t200\t400\nchr2\t600\t800"

# Tests 9-11: Complex 5-haplotype
H1="chr1\t100\t500\nchr2\t1000\t2000\nchr3\t50\t200"
H2="chr1\t200\t600\nchr2\t1500\t2500\nchr3\t100\t300"
H3="chr1\t350\t700\nchr2\t1800\t2200\nchr3\t150\t250"
H4="chr1\t400\t800\nchr2\t100\t500"
H5="chr1\t550\t900\nchr3\t180\t350"

run_case "Test9_complex_min2 (min=2)" 2 \
    h1.bed "$H1" h2.bed "$H2" h3.bed "$H3" h4.bed "$H4" h5.bed "$H5"

run_case "Test10_complex_min3 (min=3)" 3 \
    h1.bed "$H1" h2.bed "$H2" h3.bed "$H3" h4.bed "$H4" h5.bed "$H5"

run_case "Test11_complex_min4 (min=4)" 4 \
    h1.bed "$H1" h2.bed "$H2" h3.bed "$H3" h4.bed "$H4" h5.bed "$H5"

# Test 12
run_case "Test12_adjacent (min=2)" 2 \
    h1.bed "chr1\t100\t200" \
    h2.bed "chr1\t200\t300"

# Test 13
run_case "Test13_single_base (min=2)" 2 \
    h1.bed "chr1\t100\t201" \
    h2.bed "chr1\t200\t300"

# Test 14
run_case "Test14_multi_interval (min=2)" 2 \
    h1.bed "chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600" \
    h2.bed "chr1\t150\t350"

# Test 15
run_case "Test15_large_coords (min=2)" 2 \
    h1.bed "chr1\t100000000\t200000000" \
    h2.bed "chr1\t150000000\t250000000"

echo ""
echo "============================================================"
echo " Summary: $PASS passed, $FAIL failed out of $((PASS+FAIL))"
echo "============================================================"

if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
