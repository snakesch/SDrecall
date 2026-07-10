#!/usr/bin/env python3
"""Classify FN causes from raw per-read fn-hap-scanner records.

This script is intentionally driven by the raw ``*_alt_reads.tsv`` records:
``FN_site, read_qname, mapping_chrom, mapping_pos, MAPQ, source_BAM,
kmer_offset, strand``.  It must not consume stale ``cat1_preliminary`` or
``final_category`` columns from older per-FN summary tables.

Decision tree implemented here:
  * input_bam == 0 and raw_bam > 0 -> Cat1
  * input_bam == 0 and raw_bam == 0 and pacbio > 0 -> Cat1b
  * no evidence in input/raw/PacBio -> Cat1c
  * input_bam > 0 -> Cat2-5 candidate.  If --run-dir is supplied, qnames are
    traced into relevant realigned only_RG BAMs and NFC coverage is checked.

For Cat2-5 candidates, relevant RGs/subgroups are identified from
``realign_groups/RG*/RG*_related_homo_regions.bed`` FC labels.  Matching NFC
beds are then ``realign_groups/RG*/RGX_SUB.nfc.bed``.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


SOURCE_INPUT = "input_bam"
SOURCE_RAW = "raw_bam"
SOURCE_PACBIO = "pacbio"


@dataclass(frozen=True)
class Site:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def bed_pos(self) -> int:
        return self.pos - 1


@dataclass(frozen=True)
class BedInterval:
    chrom: str
    start: int
    end: int

    def contains_1based(self, chrom: str, pos: int) -> bool:
        return self.chrom == chrom and self.start <= pos - 1 < self.end


@dataclass(frozen=True)
class FcGroup:
    rg: str
    subgroup: str
    interval: BedInterval
    nfc_bed: Path


@dataclass
class TraceRecord:
    site: str
    rg: str
    input_qname: str
    bam_qname: str
    chrom: str
    pos: int
    mapq: int
    cigar: str
    overlaps_fn_pos: bool
    supports_alt_at_fn: bool


def read_sites(path: Path) -> list[Site]:
    sites: list[Site] = []
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if len(row) < 5:
                raise ValueError(f"Expected >=5 columns in {path}: {row}")
            sites.append(Site(row[0], row[1], int(row[2]), row[3], row[4]))
    return sites


def read_alt_records(path: Path) -> dict[str, list[dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "FN_site",
            "read_qname",
            "mapping_chrom",
            "mapping_pos",
            "MAPQ",
            "source_BAM",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        by_site: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in reader:
            by_site[row["FN_site"]].append(row)
        return by_site


def source_counts(records: Iterable[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in records:
        counts[row["source_BAM"]] += 1
    return counts


def evidence_category(counts: Counter[str]) -> str:
    if counts[SOURCE_INPUT] > 0:
        return "Cat2-5_candidate_input_ALT"
    if counts[SOURCE_RAW] > 0:
        return "Cat1_raw_only_input_loss"
    if counts[SOURCE_PACBIO] > 0:
        return "Cat1b_PacBio_only"
    return "Cat1c_no_evidence"


def read_bed(path: Path) -> list[BedInterval]:
    intervals: list[BedInterval] = []
    if not path.exists():
        return intervals
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            intervals.append(BedInterval(fields[0], int(fields[1]), int(fields[2])))
    return intervals


def interval_contains(intervals: Iterable[BedInterval], chrom: str, pos: int) -> bool:
    return any(interval.contains_1based(chrom, pos) for interval in intervals)


def find_fc_groups(run_dir: Path, sites: list[Site]) -> dict[str, list[FcGroup]]:
    site_by_key = {site.key: site for site in sites}
    groups: dict[str, list[FcGroup]] = {site.key: [] for site in sites}
    rg_root = run_dir / "realign_groups"
    for related in sorted(rg_root.glob("RG*/RG*_related_homo_regions.bed")):
        rg = related.parent.name
        with related.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or not fields[6].startswith("FC:"):
                    continue
                interval = BedInterval(fields[0], int(fields[1]), int(fields[2]))
                subgroup = fields[6].split(":", 1)[1]
                for key, site in site_by_key.items():
                    if interval.contains_1based(site.chrom, site.pos):
                        groups[key].append(
                            FcGroup(
                                rg=rg,
                                subgroup=subgroup,
                                interval=interval,
                                nfc_bed=related.parent / f"{subgroup}.nfc.bed",
                            )
                        )
    return groups


def qname_base(qname: str) -> str:
    return re.sub(r":RG[0-9]+$", "", qname)


def qname_variants(qname: str, rg: str) -> set[str]:
    base = qname_base(qname)
    return {base, f"{base}:{rg}"}


def cigar_ref_len(cigar: str) -> int:
    if cigar == "*":
        return 0
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(length)
    return total


def record_overlaps_pos(chrom: str, pos: int, cigar: str, site: Site) -> bool:
    if chrom != site.chrom:
        return False
    ref_len = cigar_ref_len(cigar)
    if ref_len == 0:
        return pos == site.pos
    return pos <= site.pos <= pos + ref_len - 1


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    if cigar == "*":
        return []
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def variant_support_from_alignment(fields: list[str], site: Site) -> bool:
    """Return whether a realigned SAM record directly supports the FN ALT allele.

    This covers SNVs/MNPs and simple left-anchored insertions/deletions.  Complex
    alleles that cannot be reconstructed from the local CIGAR return False; the
    raw fn-hap-scanner qname remains available in the trace table for manual
    review.
    """
    chrom = fields[2]
    if chrom != site.chrom:
        return False

    ref_pos = int(fields[3])
    read_pos = 0
    seq = fields[9].upper()
    bases: dict[int, str] = {}
    insertions_after: dict[int, str] = defaultdict(str)

    for length, op in cigar_ops(fields[5]):
        if op in {"M", "=", "X"}:
            for offset in range(length):
                if read_pos + offset < len(seq):
                    bases[ref_pos + offset] = seq[read_pos + offset]
            ref_pos += length
            read_pos += length
        elif op == "I":
            insertions_after[ref_pos - 1] += seq[read_pos : read_pos + length]
            read_pos += length
        elif op in {"D", "N"}:
            for offset in range(length):
                bases[ref_pos + offset] = "-"
            ref_pos += length
        elif op == "S":
            read_pos += length
        elif op in {"H", "P"}:
            continue

    ref = site.ref.upper()
    alt = site.alt.upper()
    pos = site.pos

    if len(ref) == len(alt):
        read_allele = "".join(bases.get(pos + i, "") for i in range(len(ref)))
        return read_allele == alt

    if len(alt) > len(ref) and alt.startswith(ref):
        prefix = "".join(bases.get(pos + i, "") for i in range(len(ref)))
        inserted = insertions_after.get(pos + len(ref) - 1, "")
        return prefix == ref and inserted.startswith(alt[len(ref) :])

    if len(ref) > len(alt) and ref.startswith(alt):
        prefix = "".join(bases.get(pos + i, "") for i in range(len(alt)))
        deleted = "".join(bases.get(pos + i, "") for i in range(len(alt), len(ref)))
        return prefix == alt and deleted == "-" * (len(ref) - len(alt))

    return False


def run_samtools_view_qnames(bam: Path, qnames: set[str]) -> list[list[str]]:
    if not qnames:
        return []
    with tempfile.NamedTemporaryFile("w", delete=False) as handle:
        qname_path = Path(handle.name)
        for qname in sorted(qnames):
            handle.write(f"{qname}\n")
    try:
        proc = subprocess.run(
            ["samtools", "view", "-N", str(qname_path), str(bam)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    finally:
        qname_path.unlink(missing_ok=True)

    rows: list[list[str]] = []
    for line in proc.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) >= 6:
            rows.append(fields)
    return rows


def filtered_mpileup_depth_at(bam: Path, site: Site) -> int:
    region = f"{site.chrom}:{site.pos}-{site.pos}"
    proc = subprocess.run(
        ["samtools", "mpileup", "-A", "-q", "10", "-Q", "15", "-r", region, str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 4 and fields[0] == site.chrom and int(fields[1]) == site.pos:
            return int(fields[3])
    return 0


def parse_mpileup_alt_count(bases: str, site: Site) -> int:
    """Count exact FN ALT observations in one mpileup base string."""
    ref = site.ref.upper()
    alt = site.alt.upper()
    alt_count = 0
    i = 0
    while i < len(bases):
        base = bases[i]
        if base == "^":
            i += 2
            continue
        if base == "$":
            i += 1
            continue

        base_matches_anchor = base in ".,ACGTNacgtn*#"
        if len(ref) == len(alt):
            if len(ref) == 1 and base.upper() == alt:
                alt_count += 1
        i += 1

        if not base_matches_anchor:
            continue
        while i < len(bases) and bases[i] in "+-":
            op = bases[i]
            i += 1
            j = i
            while j < len(bases) and bases[j].isdigit():
                j += 1
            if j == i:
                continue
            length = int(bases[i:j])
            seq = bases[j : j + length].upper()
            i = j + length

            if len(alt) > len(ref) and alt.startswith(ref):
                if op == "+" and seq == alt[len(ref) :]:
                    alt_count += 1
            elif len(ref) > len(alt) and ref.startswith(alt):
                if op == "-" and seq == ref[len(alt) :]:
                    alt_count += 1

    return alt_count


def mpileup_exact_alt_at(
    bam: Path, site: Site, ref_fasta: Path | None
) -> tuple[int, int]:
    region = f"{site.chrom}:{site.pos}-{site.pos}"
    cmd = ["samtools", "mpileup", "-A", "-q", "10", "-Q", "15", "-r", region]
    if ref_fasta is not None:
        cmd.extend(["-f", str(ref_fasta)])
    cmd.append(str(bam))
    proc = subprocess.run(
        cmd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 5 and fields[0] == site.chrom and int(fields[1]) == site.pos:
            depth = int(fields[3])
            alt_count = parse_mpileup_alt_count(fields[4], site)
            return depth, alt_count
    return 0, 0


def trace_qnames(
    run_dir: Path,
    ref_fasta: Path | None,
    sites: list[Site],
    records_by_site: dict[str, list[dict[str, str]]],
    fc_groups_by_site: dict[str, list[FcGroup]],
) -> tuple[dict[str, list[TraceRecord]], dict[tuple[str, str], tuple[int, int]]]:
    site_lookup = {site.key: site for site in sites}
    rg_to_qnames: dict[str, set[str]] = defaultdict(set)
    rg_base_to_sites: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

    for site in sites:
        input_qnames = {
            row["read_qname"]
            for row in records_by_site.get(site.key, [])
            if row["source_BAM"] == SOURCE_INPUT
        }
        if not input_qnames:
            continue
        for group in fc_groups_by_site.get(site.key, []):
            for qname in input_qnames:
                base = qname_base(qname)
                rg_base_to_sites[group.rg][base].add(site.key)
                rg_to_qnames[group.rg].update(qname_variants(qname, group.rg))

    traces_by_site: dict[str, list[TraceRecord]] = defaultdict(list)
    mpileup_by_site_rg: dict[tuple[str, str], tuple[int, int]] = {}
    recall_root = run_dir / "recall_results"

    for rg, qnames in sorted(rg_to_qnames.items(), key=lambda item: item[0]):
        bam = recall_root / f"HG002.sdrecall.only_{rg}.raw.bam"
        if not bam.exists():
            continue
        for fields in run_samtools_view_qnames(bam, qnames):
            bam_qname = fields[0]
            base = qname_base(bam_qname)
            for site_key in rg_base_to_sites[rg].get(base, []):
                site = site_lookup[site_key]
                trace = TraceRecord(
                    site=site_key,
                    rg=rg,
                    input_qname=base,
                    bam_qname=bam_qname,
                    chrom=fields[2],
                    pos=int(fields[3]),
                    mapq=int(fields[4]),
                    cigar=fields[5],
                    overlaps_fn_pos=record_overlaps_pos(
                        fields[2], int(fields[3]), fields[5], site
                    ),
                    supports_alt_at_fn=variant_support_from_alignment(fields, site),
                )
                traces_by_site[site_key].append(trace)

    for site_key, traces in traces_by_site.items():
        site = site_lookup[site_key]
        for rg in {trace.rg for trace in traces if trace.overlaps_fn_pos}:
            bam = recall_root / f"HG002.sdrecall.only_{rg}.raw.bam"
            mpileup_by_site_rg[(site_key, rg)] = mpileup_exact_alt_at(
                bam, site, ref_fasta
            )

    return traces_by_site, mpileup_by_site_rg


def nfc_status(
    site_key: str,
    raw_records: list[dict[str, str]],
    fc_groups: list[FcGroup],
    nfc_cache: dict[Path, list[BedInterval]],
) -> tuple[int, int, int, str]:
    input_rows = [row for row in raw_records if row["source_BAM"] == SOURCE_INPUT]
    outside_fc = 0
    nfc_covered = 0
    nfc_uncovered = 0

    for row in input_rows:
        chrom = row["mapping_chrom"]
        pos = int(row["mapping_pos"])
        in_any_fc = any(group.interval.contains_1based(chrom, pos) for group in fc_groups)
        if in_any_fc:
            continue
        outside_fc += 1
        covered = False
        for group in fc_groups:
            if group.nfc_bed not in nfc_cache:
                nfc_cache[group.nfc_bed] = read_bed(group.nfc_bed)
            if interval_contains(nfc_cache[group.nfc_bed], chrom, pos):
                covered = True
                break
        if covered:
            nfc_covered += 1
        else:
            nfc_uncovered += 1

    if outside_fc == 0:
        status = "no_paralogous_input_locus_outside_FC"
    elif nfc_covered == 0:
        status = "paralogous_input_loci_not_in_NFC"
    elif nfc_uncovered == 0:
        status = "paralogous_input_loci_all_in_NFC"
    else:
        status = "paralogous_input_loci_partly_in_NFC"

    return outside_fc, nfc_covered, nfc_uncovered, status


def summarize_trace(
    site_key: str,
    traces: list[TraceRecord],
    mpileup_by_site_rg: dict[tuple[str, str], tuple[int, int]],
) -> dict[str, str]:
    any_qnames = {trace.input_qname for trace in traces}
    overlap_qnames_by_rg: dict[str, set[str]] = defaultdict(set)
    trace_alt_qnames_by_rg: dict[str, set[str]] = defaultdict(set)
    all_mapqs = [trace.mapq for trace in traces]
    loci = Counter(f"{trace.chrom}:{trace.pos}" for trace in traces)

    for trace in traces:
        if trace.overlaps_fn_pos:
            overlap_qnames_by_rg[trace.rg].add(trace.input_qname)
        if trace.supports_alt_at_fn:
            trace_alt_qnames_by_rg[trace.rg].add(trace.input_qname)

    rg_af_parts: list[str] = []
    min_af = ""
    max_af = ""
    af_values: list[float] = []
    mpileup_alt_total = 0
    for rg in sorted({trace.rg for trace in traces if trace.overlaps_fn_pos}):
        depth, alt_count = mpileup_by_site_rg.get((site_key, rg), (0, 0))
        if alt_count == 0:
            continue
        mpileup_alt_total += alt_count
        af = (alt_count / depth) if depth else 0.0
        af_values.append(af)
        rg_af_parts.append(f"{rg}:{alt_count}/{depth}:{af:.4f}")
    if af_values:
        min_af = f"{min(af_values):.4f}"
        max_af = f"{max(af_values):.4f}"

    return {
        "realigned_qnames_any": str(len(any_qnames)),
        "realigned_records_any": str(len(traces)),
        "realigned_qnames_at_fn": str(
            len(set().union(*overlap_qnames_by_rg.values()))
            if overlap_qnames_by_rg
            else 0
        ),
        "realigned_alt_qnames_at_fn": str(mpileup_alt_total),
        "exact_mpileup_alt_observations_at_fn": str(mpileup_alt_total),
        "trace_cigar_alt_qnames_at_fn": str(
            len(set().union(*trace_alt_qnames_by_rg.values()))
            if trace_alt_qnames_by_rg
            else 0
        ),
        "realigned_rg_af": ";".join(rg_af_parts),
        "min_realigned_rg_af": min_af,
        "max_realigned_rg_af": max_af,
        "realigned_mapq_min": str(min(all_mapqs)) if all_mapqs else "",
        "realigned_mapq_median": (
            f"{statistics.median(all_mapqs):.1f}" if all_mapqs else ""
        ),
        "top_realigned_loci": ";".join(
            f"{locus}:{count}" for locus, count in loci.most_common(5)
        ),
    }


def decide_final_category(
    evidence_cat: str,
    trace_summary: dict[str, str],
    nfc_status_text: str,
    nfc_uncovered: int,
    low_af_threshold: float,
    run_dir_supplied: bool,
) -> tuple[str, str]:
    if evidence_cat != "Cat2-5_candidate_input_ALT":
        return evidence_cat, "final_from_source_BAM_evidence"
    if not run_dir_supplied:
        return "pending_qname_trace", "input_ALT_requires_realigned_RG_trace"

    realigned_any = int(trace_summary["realigned_qnames_any"])
    realigned_at_fn = int(trace_summary["realigned_qnames_at_fn"])
    realigned_alt_at_fn = int(trace_summary["exact_mpileup_alt_observations_at_fn"])

    if realigned_alt_at_fn > 0:
        min_af = float(trace_summary["min_realigned_rg_af"] or "0")
        if min_af <= low_af_threshold:
            return "Cat5_low_per_RG_AF", "input_qnames_realigned_to_FN_low_AF"
        return (
            "pending_nonlow_AF_caller_or_representation",
            "input_qnames_realigned_to_FN_but_AF_not_low",
        )

    if realigned_at_fn > 0:
        return (
            "pending_realigned_at_FN_without_ALT_support",
            "input_qnames_overlap_FN_but_ALT_support_not_confirmed",
        )

    if realigned_any > 0:
        return "Cat4_realignment_ambiguity", "input_qnames_realigned_away_from_FN"

    if nfc_status_text == "paralogous_input_loci_not_in_NFC":
        return "Cat2_missing_NFC_paralog_capture", "input_qnames_absent_from_RG_and_NFC_missing"
    if nfc_uncovered > 0:
        return "Cat2_partial_NFC_paralog_capture", "input_qnames_absent_from_RG_and_NFC_partly_missing"

    return "Cat3_absent_despite_FC_or_NFC_capture", "input_qnames_absent_from_RG_but_source_loci_captured"


def qname_list(rows: Iterable[dict[str, str]], source: str) -> str:
    qnames = sorted({row["read_qname"] for row in rows if row["source_BAM"] == source})
    return ";".join(qnames)


def locus_list(rows: Iterable[dict[str, str]], source: str) -> str:
    loci = sorted(
        {
            f"{row['read_qname']}@{row['mapping_chrom']}:{row['mapping_pos']}:MAPQ{row['MAPQ']}"
            for row in rows
            if row["source_BAM"] == source
        }
    )
    return ";".join(loci)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--fn-sites", required=True, type=Path)
    parser.add_argument("--alt-reads", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--ref-fasta", type=Path)
    parser.add_argument("--trace-output", type=Path)
    parser.add_argument("--low-af-threshold", type=float, default=0.20)
    args = parser.parse_args()

    sites = read_sites(args.fn_sites)
    records_by_site = read_alt_records(args.alt_reads)

    fc_groups_by_site: dict[str, list[FcGroup]] = {site.key: [] for site in sites}
    traces_by_site: dict[str, list[TraceRecord]] = defaultdict(list)
    mpileup_by_site_rg: dict[tuple[str, str], tuple[int, int]] = {}

    if args.run_dir:
        fc_groups_by_site = find_fc_groups(args.run_dir, sites)
        traces_by_site, mpileup_by_site_rg = trace_qnames(
            args.run_dir, args.ref_fasta, sites, records_by_site, fc_groups_by_site
        )

    nfc_cache: dict[Path, list[BedInterval]] = {}
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", newline="") as handle:
        fieldnames = [
            "assembly",
            "FN_site",
            "input_reads",
            "raw_reads",
            "pacbio_reads",
            "input_qnames",
            "input_loci",
            "raw_loci",
            "pacbio_loci",
            "source_evidence_category",
            "fc_groups",
            "paralogous_input_loci_outside_FC",
            "nfc_covered_loci",
            "nfc_uncovered_loci",
            "nfc_status",
            "realigned_qnames_any",
            "realigned_records_any",
            "realigned_qnames_at_fn",
            "realigned_alt_qnames_at_fn",
            "exact_mpileup_alt_observations_at_fn",
            "trace_cigar_alt_qnames_at_fn",
            "realigned_rg_af",
            "min_realigned_rg_af",
            "max_realigned_rg_af",
            "realigned_mapq_min",
            "realigned_mapq_median",
            "top_realigned_loci",
            "final_category",
            "decision_reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for site in sites:
            rows = records_by_site.get(site.key, [])
            counts = source_counts(rows)
            evidence_cat = evidence_category(counts)
            fc_groups = fc_groups_by_site.get(site.key, [])
            outside_fc, nfc_covered, nfc_uncovered, nfc_text = nfc_status(
                site.key, rows, fc_groups, nfc_cache
            )
            trace_summary = summarize_trace(
                site.key, traces_by_site.get(site.key, []), mpileup_by_site_rg
            )
            final_category, reason = decide_final_category(
                evidence_cat,
                trace_summary,
                nfc_text,
                nfc_uncovered,
                args.low_af_threshold,
                args.run_dir is not None,
            )
            writer.writerow(
                {
                    "assembly": args.assembly,
                    "FN_site": site.key,
                    "input_reads": counts[SOURCE_INPUT],
                    "raw_reads": counts[SOURCE_RAW],
                    "pacbio_reads": counts[SOURCE_PACBIO],
                    "input_qnames": qname_list(rows, SOURCE_INPUT),
                    "input_loci": locus_list(rows, SOURCE_INPUT),
                    "raw_loci": locus_list(rows, SOURCE_RAW),
                    "pacbio_loci": locus_list(rows, SOURCE_PACBIO),
                    "source_evidence_category": evidence_cat,
                    "fc_groups": ";".join(
                        f"{group.subgroup}:{group.interval.chrom}:{group.interval.start}-{group.interval.end}"
                        for group in fc_groups
                    ),
                    "paralogous_input_loci_outside_FC": outside_fc,
                    "nfc_covered_loci": nfc_covered,
                    "nfc_uncovered_loci": nfc_uncovered,
                    "nfc_status": nfc_text,
                    **trace_summary,
                    "final_category": final_category,
                    "decision_reason": reason,
                }
            )

    if args.trace_output:
        args.trace_output.parent.mkdir(parents=True, exist_ok=True)
        with args.trace_output.open("w", newline="") as handle:
            fieldnames = [
                "FN_site",
                "rg",
                "input_qname",
                "bam_qname",
                "realigned_chrom",
                "realigned_pos",
                "realigned_MAPQ",
                "CIGAR",
                "overlaps_fn_pos",
                "supports_alt_at_fn",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for site_key in sorted(traces_by_site):
                for trace in traces_by_site[site_key]:
                    writer.writerow(
                        {
                            "FN_site": trace.site,
                            "rg": trace.rg,
                            "input_qname": trace.input_qname,
                            "bam_qname": trace.bam_qname,
                            "realigned_chrom": trace.chrom,
                            "realigned_pos": trace.pos,
                            "realigned_MAPQ": trace.mapq,
                            "CIGAR": trace.cigar,
                            "overlaps_fn_pos": "yes" if trace.overlaps_fn_pos else "no",
                            "supports_alt_at_fn": "yes"
                            if trace.supports_alt_at_fn
                            else "no",
                        }
                    )


if __name__ == "__main__":
    main()
