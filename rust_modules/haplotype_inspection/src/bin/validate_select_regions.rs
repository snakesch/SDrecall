/// Validation binary: runs select_regions_with_min_haplotypes on the same inputs
/// as the BEDOPS pipeline in validate_select_regions.sh, printing BED-format output
/// so the two can be compared line-by-line.
use std::collections::HashMap;
use haplotype_inspection::identify_misaligned_haps::select_regions_with_min_haplotypes;

fn make_intervals(data: &[(i32, &[(&str, i64, i64)])]) -> HashMap<i32, Vec<(String, i64, i64)>> {
    let mut m = HashMap::new();
    for &(hid, ivs) in data {
        m.insert(
            hid,
            ivs.iter()
                .map(|(c, s, e)| (c.to_string(), *s, *e))
                .collect(),
        );
    }
    m
}

fn print_bed(result: &Option<Vec<(String, i64, i64)>>) {
    if let Some(regions) = result {
        for (c, s, e) in regions {
            println!("{}\t{}\t{}", c, s, e);
        }
    }
}

fn run_test(name: &str, data: &[(i32, &[(&str, i64, i64)])], min_k: usize) {
    let intervals = make_intervals(data);
    let result = select_regions_with_min_haplotypes(&intervals, min_k);
    println!("### {} (min={})", name, min_k);
    print_bed(&result);
    println!("---");
}

fn main() {
    // Test 1: Simple 2-haplotype overlap
    run_test(
        "Test1_simple_2hap",
        &[
            (1, &[("chr1", 100, 300)]),
            (2, &[("chr1", 200, 400)]),
        ],
        2,
    );

    // Test 2: Staircase, min=2
    run_test(
        "Test2_staircase_min2",
        &[
            (1, &[("chr1", 100, 400)]),
            (2, &[("chr1", 200, 500)]),
            (3, &[("chr1", 300, 600)]),
        ],
        2,
    );

    // Test 3: Staircase, min=3
    run_test(
        "Test3_staircase_min3",
        &[
            (1, &[("chr1", 100, 400)]),
            (2, &[("chr1", 200, 500)]),
            (3, &[("chr1", 300, 600)]),
        ],
        3,
    );

    // Test 4: No overlap
    run_test(
        "Test4_no_overlap",
        &[
            (1, &[("chr1", 100, 200)]),
            (2, &[("chr1", 300, 400)]),
        ],
        2,
    );

    // Test 5: Identical intervals (3 haplotypes)
    run_test(
        "Test5_identical_3hap",
        &[
            (1, &[("chr1", 100, 500)]),
            (2, &[("chr1", 100, 500)]),
            (3, &[("chr1", 100, 500)]),
        ],
        2,
    );

    // Test 6: Contained interval
    run_test(
        "Test6_contained",
        &[
            (1, &[("chr1", 100, 600)]),
            (2, &[("chr1", 200, 400)]),
        ],
        2,
    );

    // Test 7: Disjoint overlap zones
    run_test(
        "Test7_disjoint_zones",
        &[
            (1, &[("chr1", 100, 200), ("chr1", 500, 600)]),
            (2, &[("chr1", 150, 250), ("chr1", 550, 650)]),
        ],
        2,
    );

    // Test 8: Multi-chromosome
    run_test(
        "Test8_multi_chrom",
        &[
            (1, &[("chr1", 100, 300), ("chr2", 500, 700)]),
            (2, &[("chr1", 200, 400), ("chr2", 600, 800)]),
        ],
        2,
    );

    // Test 9: Complex — 5 haplotypes, 3 chromosomes, min=2
    let complex_data: &[(i32, &[(&str, i64, i64)])] = &[
        (1, &[("chr1", 100, 500), ("chr2", 1000, 2000), ("chr3", 50, 200)]),
        (2, &[("chr1", 200, 600), ("chr2", 1500, 2500), ("chr3", 100, 300)]),
        (3, &[("chr1", 350, 700), ("chr2", 1800, 2200), ("chr3", 150, 250)]),
        (4, &[("chr1", 400, 800), ("chr2", 100, 500)]),
        (5, &[("chr1", 550, 900), ("chr3", 180, 350)]),
    ];
    run_test("Test9_complex_min2", complex_data, 2);
    run_test("Test10_complex_min3", complex_data, 3);
    run_test("Test11_complex_min4", complex_data, 4);

    // Test 12: Adjacent boundary
    run_test(
        "Test12_adjacent",
        &[
            (1, &[("chr1", 100, 200)]),
            (2, &[("chr1", 200, 300)]),
        ],
        2,
    );

    // Test 13: Single base overlap
    run_test(
        "Test13_single_base",
        &[
            (1, &[("chr1", 100, 201)]),
            (2, &[("chr1", 200, 300)]),
        ],
        2,
    );

    // Test 14: Multi-interval one haplotype
    run_test(
        "Test14_multi_interval",
        &[
            (1, &[("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 500, 600)]),
            (2, &[("chr1", 150, 350)]),
        ],
        2,
    );

    // Test 15: Large coordinates
    run_test(
        "Test15_large_coords",
        &[
            (1, &[("chr1", 100000000, 200000000)]),
            (2, &[("chr1", 150000000, 250000000)]),
        ],
        2,
    );
}
