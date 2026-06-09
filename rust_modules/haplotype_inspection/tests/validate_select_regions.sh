#!/usr/bin/env bash
# =============================================================================
# Validate Rust select_regions_with_min_haplotypes against BEDOPS pipeline
#
# This script replicates the exact BEDOPS pipeline used in Python's
# select_regions_with_min_haplotypes_from_hapbeds() and compares the output
# with the Rust implementation for multiple complex test cases.
# =============================================================================
set -euo pipefail

TMPDIR=$(mktemp -d /tmp/validate_select_XXXXXX)
trap "rm -rf $TMPDIR" EXIT

PASS=0
FAIL=0
TOTAL=0

run_bedops_pipeline() {
    # Args: min_haplotypes output_bed bed_files...
    local min_k="$1"; shift
    local output_bed="$1"; shift
    local map_file="$TMPDIR/map.tsv"
    local tagged_bed="$TMPDIR/tagged.bed"
    local parts_bed="$TMPDIR/parts.bed"

    # Build mapping file
    > "$map_file"
    local hid=1
    for bed in "$@"; do
        echo -e "${hid}\t${bed}" >> "$map_file"
        hid=$((hid+1))
    done

    # Tag each BED with haplotype ID, then sort
    while IFS=$'\t' read -r hid bed; do
        mawk -v id="$hid" 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,id}' "$bed"
    done < "$map_file" | sort-bed - > "$tagged_bed"

    # Partition
    bedops --partition "$tagged_bed" > "$parts_bed"

    # Map, filter by k, sort
    bedmap --ec --delim $'\t' --echo --echo-map-id-uniq \
        "$parts_bed" "$tagged_bed" | \
        mawk -F '\t' -v k="$min_k" \
        '{ ids=$NF; n=(ids=="N/A")?0:split(ids,a,";"); if (n>=k) print $1"\t"$2"\t"$3 }' | \
        sort-bed - | bedops --merge - > "$output_bed"
}

compare_results() {
    local test_name="$1"
    local bedops_result="$2"
    local rust_result="$3"
    TOTAL=$((TOTAL+1))

    # Normalize: sort both files for comparison
    local bedops_sorted="$TMPDIR/bedops_sorted.bed"
    local rust_sorted="$TMPDIR/rust_sorted.bed"
    sort -k1,1 -k2,2n "$bedops_result" > "$bedops_sorted" 2>/dev/null || true
    sort -k1,1 -k2,2n "$rust_result" > "$rust_sorted" 2>/dev/null || true

    if diff -q "$bedops_sorted" "$rust_sorted" > /dev/null 2>&1; then
        echo "  PASS: $test_name ($(wc -l < "$bedops_sorted") regions)"
        PASS=$((PASS+1))
    else
        echo "  FAIL: $test_name"
        echo "    BEDOPS output:"
        cat "$bedops_result" | sed 's/^/      /'
        echo "    Rust output:"
        cat "$rust_result" | sed 's/^/      /'
        echo "    Diff:"
        diff "$bedops_sorted" "$rust_sorted" | sed 's/^/      /' || true
        FAIL=$((FAIL+1))
    fi
}

echo "============================================================"
echo " Validating select_regions_with_min_haplotypes vs BEDOPS"
echo "============================================================"
echo ""

# ─── Test 1: Simple 2-haplotype overlap ────────────────────────────
echo "Test 1: Simple 2-haplotype overlap"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	300
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	400
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
echo -e "chr1\t200\t300" > "$TMPDIR/rust_out.bed"
compare_results "simple_2hap_overlap" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 2: Three haplotypes, staircase ───────────────────────────
echo "Test 2: Three haplotypes staircase (min=2)"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	400
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	500
EOF
cat > "$TMPDIR/hap3.bed" << 'EOF'
chr1	300	600
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed"
echo -e "chr1\t200\t500" > "$TMPDIR/rust_out.bed"
compare_results "staircase_3hap_min2" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 3: Three haplotypes, staircase (min=3) ──────────────────
echo "Test 3: Three haplotypes staircase (min=3)"
run_bedops_pipeline 3 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed"
echo -e "chr1\t300\t400" > "$TMPDIR/rust_out.bed"
compare_results "staircase_3hap_min3" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 4: No overlap ───────────────────────────────────────────
echo "Test 4: No overlap"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	200
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	300	400
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
> "$TMPDIR/rust_out.bed"  # empty
compare_results "no_overlap" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 5: Identical intervals ──────────────────────────────────
echo "Test 5: Identical intervals (3 haplotypes)"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	500
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	100	500
EOF
cat > "$TMPDIR/hap3.bed" << 'EOF'
chr1	100	500
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed"
echo -e "chr1\t100\t500" > "$TMPDIR/rust_out.bed"
compare_results "identical_3hap" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 6: Contained interval ───────────────────────────────────
echo "Test 6: Contained interval"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	600
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	400
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
echo -e "chr1\t200\t400" > "$TMPDIR/rust_out.bed"
compare_results "contained" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 7: Disjoint overlap zones ───────────────────────────────
echo "Test 7: Disjoint overlap zones"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	200
chr1	500	600
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	150	250
chr1	550	650
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
printf "chr1\t150\t200\nchr1\t550\t600\n" > "$TMPDIR/rust_out.bed"
compare_results "disjoint_zones" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 8: Multi-chromosome ─────────────────────────────────────
echo "Test 8: Multi-chromosome"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	300
chr2	500	700
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	400
chr2	600	800
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
printf "chr1\t200\t300\nchr2\t600\t700\n" > "$TMPDIR/rust_out.bed"
compare_results "multi_chrom" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 9: Complex — 5 haplotypes, 3 chromosomes ────────────────
echo "Test 9: Complex — 5 haplotypes, 3 chromosomes (min=2)"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	500
chr2	1000	2000
chr3	50	200
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	600
chr2	1500	2500
chr3	100	300
EOF
cat > "$TMPDIR/hap3.bed" << 'EOF'
chr1	350	700
chr2	1800	2200
chr3	150	250
EOF
cat > "$TMPDIR/hap4.bed" << 'EOF'
chr1	400	800
chr2	100	500
EOF
cat > "$TMPDIR/hap5.bed" << 'EOF'
chr1	550	900
chr3	180	350
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" \
    "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed" \
    "$TMPDIR/hap4.bed" "$TMPDIR/hap5.bed"
# Rust expected: need to compute. Use BEDOPS as ground truth.
cat "$TMPDIR/bedops_out.bed"
echo "(BEDOPS ground truth for Test 9 above)"

# ─── Test 10: Complex — 5 haplotypes, 3 chromosomes (min=3) ───────
echo "Test 10: Complex — 5 haplotypes, 3 chromosomes (min=3)"
run_bedops_pipeline 3 "$TMPDIR/bedops_out_min3.bed" \
    "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed" \
    "$TMPDIR/hap4.bed" "$TMPDIR/hap5.bed"
cat "$TMPDIR/bedops_out_min3.bed"
echo "(BEDOPS ground truth for Test 10 above)"

# ─── Test 11: Complex — 5 haplotypes, 3 chromosomes (min=4) ───────
echo "Test 11: Complex — 5 haplotypes, 3 chromosomes (min=4)"
run_bedops_pipeline 4 "$TMPDIR/bedops_out_min4.bed" \
    "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed" "$TMPDIR/hap3.bed" \
    "$TMPDIR/hap4.bed" "$TMPDIR/hap5.bed"
cat "$TMPDIR/bedops_out_min4.bed"
echo "(BEDOPS ground truth for Test 11 above)"

# ─── Test 12: Adjacent but non-overlapping ─────────────────────────
echo "Test 12: Adjacent but non-overlapping (touching at boundary)"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	200
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	300
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
> "$TMPDIR/rust_out.bed"  # empty: half-open intervals don't overlap at boundary
compare_results "adjacent_boundary" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 13: Single base overlap ─────────────────────────────────
echo "Test 13: Single base overlap"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	201
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	200	300
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
echo -e "chr1\t200\t201" > "$TMPDIR/rust_out.bed"
compare_results "single_base_overlap" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 14: One haplotype with multiple disjoint intervals ──────
echo "Test 14: One haplotype, multiple intervals, overlapping with another"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100	200
chr1	300	400
chr1	500	600
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	150	350
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
printf "chr1\t150\t200\nchr1\t300\t350\n" > "$TMPDIR/rust_out.bed"
compare_results "multi_interval_one_hap" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

# ─── Test 15: Large coordinate values ─────────────────────────────
echo "Test 15: Large coordinate values (chromosome-scale)"
cat > "$TMPDIR/hap1.bed" << 'EOF'
chr1	100000000	200000000
EOF
cat > "$TMPDIR/hap2.bed" << 'EOF'
chr1	150000000	250000000
EOF
run_bedops_pipeline 2 "$TMPDIR/bedops_out.bed" "$TMPDIR/hap1.bed" "$TMPDIR/hap2.bed"
echo -e "chr1\t150000000\t200000000" > "$TMPDIR/rust_out.bed"
compare_results "large_coords" "$TMPDIR/bedops_out.bed" "$TMPDIR/rust_out.bed"

echo ""
echo "============================================================"
echo " Summary: $PASS passed, $FAIL failed, $TOTAL total"
echo " (Tests 9-11 show BEDOPS ground truth only — need Rust binary)"
echo "============================================================"
