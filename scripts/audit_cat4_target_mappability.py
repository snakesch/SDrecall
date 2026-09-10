#!/usr/bin/env python3
"""Audit whether hg19 Cat4 reads can align to the isolated FN-containing FC.

This audit deliberately separates two questions:

1. Does the complete paired read align to a reference containing only the FC
   that overlaps the FN coordinate?
2. Is the exact ALT-haplotype marker unique in the complete assembly reference?

An isolated-FC alignment demonstrates sequence compatibility with the target.
It does not establish the biological origin of a read because competing
paralogous loci have intentionally been removed from that alignment.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import openpyxl
import pysam

from audit_cat4_qname_paths import Interval, exact_alt_event
from update_fn_realign_allele_audit import Site, canonical_qname, parse_site


ROOT = Path("/paedyl01/disk1/yangyxt/SDrecall-rust-migration")
AUDIT_ROOT = ROOT / "test_tmp/nfc_fc_conflict_fix_20260723"
RUN_ROOT = AUDIT_ROOT / "runs/hg19/HG002_hg19_nfcfix23_SDrecall"
HAPLOTYPE_DIR = AUDIT_ROOT / "haplotypes/hg19"
REFERENCE = Path("/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta")
DEFAULT_OUTPUT = AUDIT_ROOT / "cat4_target_mappability"

OBSERVATION_COLUMNS = [
    "Assembly",
    "FN allele",
    "Query SD region",
    "Homologous counterpart region",
    "RG index",
    "FC query pair index",
    "NFC/counterpart pair index",
    "Selected input BAM ALT-haplotype qnames",
    "Original raw BAM ALT-haplotype qnames",
    "ALT-haplotype depth / total depth at true ALT site",
    "Realigned exact-allele AD distribution and FN cause",
]
QNAME_COLUMNS = [
    "Assembly",
    "FN allele",
    "ALT-haplotype qname",
    "Evidence source",
    "Selected input BAM alignment records",
    "Original raw BAM alignment records",
    "Realigned BAM alignment records",
    "Should be extracted (FC/NFC pair)",
    "Present in extracted RG FASTQ",
    "Realigned read pair produces exact FN ALT",
]


@dataclass(frozen=True)
class FastqRecord:
    name: str
    sequence: str
    qualities: str


@dataclass(frozen=True)
class FastqPair:
    read1: FastqRecord
    read2: FastqRecord


@dataclass(frozen=True)
class TargetFc:
    rg: str
    label: str
    interval: Interval

    @property
    def contig(self) -> str:
        return self.label.replace(":", "_")

    def display(self) -> str:
        return f"{self.rg} | {self.label} | {self.interval.display()}"


@dataclass(frozen=True)
class ForcedRecord:
    mate: int
    target_fc: str
    genomic_interval: Interval
    mapq: int
    cigar: str
    alignment_score: int | None
    edit_distance: int | None
    identity: float | None
    overlaps_fn: bool
    exact_alt: bool

    def display(self) -> str:
        score = "." if self.alignment_score is None else str(self.alignment_score)
        edits = "." if self.edit_distance is None else str(self.edit_distance)
        identity = "." if self.identity is None else f"{100 * self.identity:.2f}%"
        return (
            f"R{self.mate}:{self.target_fc}:{self.genomic_interval.display()}:"
            f"MQ{self.mapq}:AS{score}:NM{edits}:identity={identity}:"
            f"{self.cigar}:overlaps_FN={'yes' if self.overlaps_fn else 'no'}:"
            f"exact_ALT={'yes' if self.exact_alt else 'no'}"
        )


@dataclass(frozen=True)
class QnameResult:
    assembly: str
    site: Site
    qname: str
    target_fcs: tuple[TargetFc, ...]
    marker_mates: tuple[int, ...]
    marker_occurrences_in_pair: int
    records: tuple[ForcedRecord, ...]
    both_mates_mapped: bool
    proper_pair: bool

    @property
    def marker_records(self) -> tuple[ForcedRecord, ...]:
        return tuple(record for record in self.records if record.mate in self.marker_mates)

    @property
    def marker_reaches_fn(self) -> bool:
        return any(record.overlaps_fn for record in self.marker_records)

    @property
    def exact_alt(self) -> bool:
        return any(record.exact_alt for record in self.marker_records)

    def assertion(self) -> str:
        if self.marker_records:
            records = ";".join(record.display() for record in self.marker_records)
            outcome = "Yes" if self.marker_reaches_fn else "No"
            return (
                f"{outcome} | marker-bearing mate(s) aligned to isolated FN-site FC; "
                f"both mates mapped={'yes' if self.both_mates_mapped else 'no'}; "
                f"proper pair={'yes' if self.proper_pair else 'no'}; {records}"
            )
        return (
            "No | marker-bearing mate has no primary alignment to the isolated "
            f"FN-site FC; both mates mapped={'yes' if self.both_mates_mapped else 'no'}"
        )


@dataclass(frozen=True)
class MarkerHit:
    site_label: str
    marker_length: int
    chrom: str
    start: int
    end: int
    strand: str
    classifications: tuple[str, ...]

    def display(self) -> str:
        labels = ",".join(self.classifications) if self.classifications else "other locus"
        return f"{self.chrom}:{self.start}-{self.end}({self.strand})[{labels}]"


@dataclass(frozen=True)
class CompetitiveResult:
    site_label: str
    qname: str
    primary_roles: tuple[str, ...]
    forced_target_score: int | None
    competing_source_score: int | None

    @property
    def score_delta_source_minus_target(self) -> int | None:
        if self.forced_target_score is None or self.competing_source_score is None:
            return None
        return self.competing_source_score - self.forced_target_score

    @property
    def preference(self) -> str:
        delta = self.score_delta_source_minus_target
        if self.competing_source_score is None:
            return "no non-target FC alignment reported"
        if delta is None:
            return "not comparable"
        if delta > 0:
            return "observed non-target FC has higher alignment score"
        if delta < 0:
            return "FN-site FC has higher alignment score"
        return "equal alignment score"

    def display(self) -> str:
        target = "." if self.forced_target_score is None else str(self.forced_target_score)
        source = "." if self.competing_source_score is None else str(self.competing_source_score)
        delta = self.score_delta_source_minus_target
        delta_text = "." if delta is None else f"{delta:+d}"
        primaries = ",".join(self.primary_roles) if self.primary_roles else "unmapped"
        return (
            f"primary={primaries}; full marker-read AS target/source={target}/{source}; "
            f"source-minus-target={delta_text}; {self.preference}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=ROOT / "FN_CAUSE_DIAGNOSIS.xlsx")
    parser.add_argument("--output-workbook", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-root", type=Path, default=RUN_ROOT)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--seqkit", default="seqkit")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--update-workbook", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def worksheet_headers(ws: openpyxl.worksheet.worksheet.Worksheet) -> list[str]:
    return [str(cell.value) for cell in ws[1]]


def parse_rg_number(rg: str) -> int:
    if not re.fullmatch(r"RG\d+", rg):
        raise ValueError(f"invalid realignment-group label: {rg!r}")
    return int(rg[2:])


def load_target_sites(
    allele_rows: list[dict[str, str]], qname_rows: list[dict[str, str]]
) -> tuple[dict[str, Site], dict[str, set[str]]]:
    sites: dict[str, Site] = {}
    qnames: dict[str, set[str]] = defaultdict(set)
    selected = {
        row["FN allele"]
        for row in allele_rows
        if row["Assembly"] == "hg19"
        and int(row["Qnames without query-linked extraction path"]) > 0
    }
    for label in sorted(selected):
        sites[label] = parse_site("hg19", label)
    for row in qname_rows:
        label = row["FN allele"]
        if row["Assembly"] == "hg19" and label in sites:
            qnames[label].add(canonical_qname(row["ALT-haplotype qname"]))
    if set(qnames) != set(sites):
        raise ValueError(f"qname/site mismatch: sites={sorted(sites)}, qnames={sorted(qnames)}")
    return sites, qnames


def load_fc_catalog(run_root: Path, rgs: set[str]) -> dict[str, TargetFc]:
    catalog: dict[str, TargetFc] = {}
    for rg in sorted(rgs, key=parse_rg_number):
        bed = run_root / "realign_groups" / rg / f"{rg}_related_homo_regions.bed"
        with bed.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or not fields[6].startswith("FC:"):
                    continue
                label = fields[6]
                value = TargetFc(
                    rg=rg,
                    label=label,
                    interval=Interval(fields[0], int(fields[1]), int(fields[2])),
                )
                previous = catalog.get(label)
                if previous is not None and previous != value:
                    raise ValueError(f"FC label has multiple intervals: {label}")
                catalog[label] = value
    return catalog


def load_site_fcs(
    path_rows: list[dict[str, str]], sites: dict[str, Site], run_root: Path
) -> tuple[dict[str, tuple[TargetFc, ...]], dict[str, tuple[TargetFc, ...]]]:
    target_labels: dict[str, set[str]] = defaultdict(set)
    source_labels: dict[str, set[str]] = defaultdict(set)
    rgs: set[str] = set()
    for row in path_rows:
        label = row["FN allele"]
        if row["Assembly"] != "hg19" or label not in sites:
            continue
        rg = row["RG"]
        rgs.add(rg)
        target_labels[label].update(
            value for value in row["All FN-containing FCs in RG"].split(";") if value
        )
        if row["Intended FC contains FN anchor"] == "no":
            source_labels[label].add(row["Intended FC"])
    catalog = load_fc_catalog(run_root, rgs)
    target_fcs: dict[str, tuple[TargetFc, ...]] = {}
    source_fcs: dict[str, tuple[TargetFc, ...]] = {}
    for label, site in sites.items():
        target_fcs[label] = tuple(
            sorted(
                (catalog[value] for value in target_labels[label]),
                key=lambda fc: (parse_rg_number(fc.rg), fc.label),
            )
        )
        source_fcs[label] = tuple(
            sorted(
                (catalog[value] for value in source_labels[label]),
                key=lambda fc: (parse_rg_number(fc.rg), fc.label),
            )
        )
        if not target_fcs[label]:
            raise ValueError(f"no FN-containing FC found for {label}")
        if not any(
            fc.interval.contains(site.chrom, site.pos - 1) for fc in target_fcs[label]
        ):
            raise ValueError(f"reported target FC does not contain {label}")
    return target_fcs, source_fcs


def read_fastq_record(handle) -> FastqRecord | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\n")
    plus = handle.readline()
    qualities = handle.readline().rstrip("\n")
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError("malformed FASTQ record")
    if len(sequence) != len(qualities):
        raise ValueError(f"FASTQ sequence/quality length mismatch for {header.rstrip()}")
    return FastqRecord(canonical_qname(header[1:]), sequence.upper(), qualities)


def extract_fastq_pairs(
    run_root: Path, qnames_by_rg: dict[str, set[str]]
) -> dict[tuple[str, str], FastqPair]:
    pairs: dict[tuple[str, str], FastqPair] = {}
    for rg, wanted in sorted(qnames_by_rg.items(), key=lambda item: parse_rg_number(item[0])):
        result_dir = run_root / "recall_results"
        paths = (
            result_dir / f"HG002.sdrecall.only_{rg}.r1.fastq",
            result_dir / f"HG002.sdrecall.only_{rg}.r2.fastq",
        )
        with paths[0].open() as read1_handle, paths[1].open() as read2_handle:
            while True:
                read1 = read_fastq_record(read1_handle)
                read2 = read_fastq_record(read2_handle)
                if read1 is None or read2 is None:
                    if read1 is not None or read2 is not None:
                        raise ValueError(f"paired FASTQ files have different lengths for {rg}")
                    break
                if read1.name != read2.name:
                    raise ValueError(f"paired FASTQ qname mismatch: {read1.name} != {read2.name}")
                if read1.name in wanted:
                    key = (rg, read1.name)
                    value = FastqPair(read1, read2)
                    previous = pairs.get(key)
                    if previous is not None and previous != value:
                        raise ValueError(f"qname has inconsistent duplicate FASTQ records: {key}")
                    pairs[key] = value
        missing = sorted(wanted - {qname for observed_rg, qname in pairs if observed_rg == rg})
        if missing:
            raise ValueError(f"{len(missing)} requested qnames absent from {rg} FASTQs: {missing[:3]}")
    return pairs


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def count_overlapping(sequence: str, pattern: str) -> int:
    count = 0
    start = 0
    while True:
        position = sequence.find(pattern, start)
        if position < 0:
            return count
        count += 1
        start = position + 1


def marker_mates(pair: FastqPair, marker: str) -> tuple[tuple[int, ...], int]:
    reverse = reverse_complement(marker)
    mates: list[int] = []
    occurrences = 0
    for mate, sequence in ((1, pair.read1.sequence), (2, pair.read2.sequence)):
        mate_occurrences = count_overlapping(sequence, marker)
        if reverse != marker:
            mate_occurrences += count_overlapping(sequence, reverse)
        if mate_occurrences:
            mates.append(mate)
            occurrences += mate_occurrences
    return tuple(mates), occurrences


def write_fastq(path: Path, records: list[FastqRecord]) -> None:
    with path.open("w") as handle:
        for record in records:
            handle.write(f"@{record.name}\n{record.sequence}\n+\n{record.qualities}\n")


def write_target_fasta(
    path: Path, run_root: Path, target_fcs: tuple[TargetFc, ...]
) -> dict[str, tuple[TargetFc, Interval]]:
    rgs = {fc.rg for fc in target_fcs}
    if len(rgs) != 1:
        raise ValueError(f"target FCs span multiple RGs: {sorted(rgs)}")
    rg = next(iter(rgs))
    masked_path = run_root / "realign_groups" / rg / f"{rg}.masked.fasta"
    contigs: dict[str, tuple[TargetFc, Interval]] = {}
    with pysam.FastaFile(masked_path) as masked, path.open("w") as handle:
        for reference_name, length in zip(masked.references, masked.lengths, strict=True):
            chrom, separator, start_text = reference_name.rpartition(":")
            if not separator or not start_text.isdigit():
                raise ValueError(f"unexpected masked-reference contig: {reference_name}")
            interval = Interval(chrom, int(start_text), int(start_text) + length)
            overlapping = [fc for fc in target_fcs if interval.overlaps(fc.interval)]
            if not overlapping:
                continue
            fc = sorted(overlapping, key=lambda value: value.label)[0]
            contig = f"target_{len(contigs) + 1}_{fc.contig}"
            sequence = masked.fetch(reference_name).upper()
            contigs[contig] = (fc, interval)
            handle.write(f">{contig}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    if not contigs:
        raise ValueError(f"no production masked contig overlaps target FCs in {rg}")
    return contigs


def alignment_identity(read: pysam.AlignedSegment) -> float | None:
    matches = mismatches = insertions = deletions = 0
    for operation, length in read.cigartuples or []:
        if operation == 7:
            matches += length
        elif operation == 8:
            mismatches += length
        elif operation == 1:
            insertions += length
        elif operation == 2:
            deletions += length
        elif operation == 0:
            edit_distance = int(read.get_tag("NM")) if read.has_tag("NM") else 0
            denominator = max(read.query_alignment_length, 1)
            return max(0.0, 1.0 - edit_distance / denominator)
    denominator = matches + mismatches + insertions + deletions
    return matches / denominator if denominator else None


def forced_record(
    read: pysam.AlignedSegment,
    site: Site,
    target_by_contig: dict[str, tuple[TargetFc, Interval]],
) -> ForcedRecord:
    fc, reference_interval = target_by_contig[read.reference_name]
    local_site = Site(
        site.assembly,
        read.reference_name,
        site.pos - reference_interval.start,
        site.ref,
        site.alt,
    )
    exact, _ = exact_alt_event(read, local_site)
    genomic = Interval(
        reference_interval.chrom,
        reference_interval.start + read.reference_start,
        reference_interval.start + int(read.reference_end),
    )
    return ForcedRecord(
        mate=1 if read.is_read1 else 2 if read.is_read2 else 0,
        target_fc=fc.label,
        genomic_interval=genomic,
        mapq=read.mapping_quality,
        cigar=read.cigarstring or "*",
        alignment_score=int(read.get_tag("AS")) if read.has_tag("AS") else None,
        edit_distance=int(read.get_tag("NM")) if read.has_tag("NM") else None,
        identity=alignment_identity(read),
        overlaps_fn=genomic.contains(site.chrom, site.pos - 1),
        exact_alt=exact,
    )


def run_forced_alignment(
    output_dir: Path,
    site: Site,
    target_fcs: tuple[TargetFc, ...],
    qnames: set[str],
    pairs: dict[tuple[str, str], FastqPair],
    marker: str,
    minimap2: str,
    threads: int,
    run_root: Path,
) -> list[QnameResult]:
    site_dir = output_dir / re.sub(r"[^A-Za-z0-9_.-]+", "_", site.label)
    site_dir.mkdir(parents=True, exist_ok=True)
    target_fasta = site_dir / "isolated_fn_site_fc.fasta"
    read1_fastq = site_dir / "selected.r1.fastq"
    read2_fastq = site_dir / "selected.r2.fastq"
    alignment_sam = site_dir / "forced_target_fc.sam"
    minimap_log = site_dir / "minimap2.log"
    target_by_contig = write_target_fasta(target_fasta, run_root, target_fcs)

    rg = target_fcs[0].rg
    site_pairs = {qname: pairs[(rg, qname)] for qname in qnames}
    ordered_qnames = sorted(qnames)
    write_fastq(read1_fastq, [site_pairs[qname].read1 for qname in ordered_qnames])
    write_fastq(read2_fastq, [site_pairs[qname].read2 for qname in ordered_qnames])

    command = [
        minimap2,
        "-ax",
        "sr",
        "--eqx",
        "--MD",
        "-F",
        "1000",
        "--end-bonus",
        "10",
        "-t",
        str(threads),
        str(target_fasta),
        str(read1_fastq),
        str(read2_fastq),
    ]
    with alignment_sam.open("w") as output_handle, minimap_log.open("w") as log_handle:
        subprocess.run(command, check=True, stdout=output_handle, stderr=log_handle)

    records_by_qname: dict[str, list[pysam.AlignedSegment]] = defaultdict(list)
    with pysam.AlignmentFile(alignment_sam, "r") as alignment:
        for read in alignment.fetch(until_eof=True):
            qname = canonical_qname(read.query_name or "")
            if qname not in qnames or read.is_secondary or read.is_supplementary:
                continue
            records_by_qname[qname].append(read)
    results: list[QnameResult] = []
    for qname in ordered_qnames:
        pair = site_pairs[qname]
        mates, marker_occurrences = marker_mates(pair, marker)
        if not mates:
            raise ValueError(f"selected marker is absent from forced-remap FASTQ pair: {site.label} {qname}")
        primary = records_by_qname[qname]
        mapped_mates = {
            1 if read.is_read1 else 2 if read.is_read2 else 0
            for read in primary
            if not read.is_unmapped
        }
        mapped_records = tuple(
            forced_record(read, site, target_by_contig)
            for read in primary
            if not read.is_unmapped
        )
        proper_pair = any(read.is_proper_pair for read in primary if not read.is_unmapped)
        results.append(
            QnameResult(
                assembly="hg19",
                site=site,
                qname=qname,
                target_fcs=target_fcs,
                marker_mates=mates,
                marker_occurrences_in_pair=marker_occurrences,
                records=mapped_records,
                both_mates_mapped={1, 2}.issubset(mapped_mates),
                proper_pair=proper_pair,
            )
        )
    return results


def write_competing_fasta(
    path: Path,
    run_root: Path,
    target_fcs: tuple[TargetFc, ...],
    source_fcs: tuple[TargetFc, ...],
) -> dict[str, tuple[str, TargetFc]]:
    contigs: dict[str, tuple[str, TargetFc]] = {}
    rgs = {fc.rg for fc in target_fcs + source_fcs}
    if len(rgs) != 1:
        raise ValueError(f"competing FCs span multiple RGs: {sorted(rgs)}")
    rg = next(iter(rgs))
    masked_path = run_root / "realign_groups" / rg / f"{rg}.masked.fasta"
    selected: list[tuple[str, str, TargetFc]] = []
    with pysam.FastaFile(masked_path) as masked:
        for reference_name, length in zip(masked.references, masked.lengths, strict=True):
            chrom, separator, start_text = reference_name.rpartition(":")
            if not separator or not start_text.isdigit():
                raise ValueError(f"unexpected masked-reference contig: {reference_name}")
            interval = Interval(chrom, int(start_text), int(start_text) + length)
            roles: list[tuple[str, TargetFc]] = []
            for role, fcs in (("target", target_fcs), ("source", source_fcs)):
                for fc in fcs:
                    if interval.overlaps(fc.interval):
                        roles.append((role, fc))
            if not roles:
                continue
            observed_roles = {role for role, _ in roles}
            if len(observed_roles) != 1:
                raise ValueError(
                    f"one production masked contig contains target and source FCs: {reference_name}"
                )
            role = next(iter(observed_roles))
            representative_fc = sorted(
                (fc for observed_role, fc in roles if observed_role == role),
                key=lambda fc: fc.label,
            )[0]
            selected.append((role, reference_name, representative_fc))

        with path.open("w") as handle:
            for index, (role, reference_name, fc) in enumerate(selected, start=1):
                contig = f"{role}_{index}_{fc.contig}"
                sequence = masked.fetch(reference_name).upper()
                contigs[contig] = (role, fc)
                handle.write(f">{contig}\n")
                for start in range(0, len(sequence), 60):
                    handle.write(sequence[start : start + 60] + "\n")
    if {role for role, _ in contigs.values()} != {"target", "source"}:
        raise ValueError(f"failed to recover target/source production contigs for {rg}")
    return contigs


def run_competing_alignment(
    output_dir: Path,
    site: Site,
    target_fcs: tuple[TargetFc, ...],
    source_fcs: tuple[TargetFc, ...],
    forced_results: list[QnameResult],
    minimap2: str,
    threads: int,
    run_root: Path,
) -> dict[str, CompetitiveResult]:
    site_dir = output_dir / re.sub(r"[^A-Za-z0-9_.-]+", "_", site.label)
    competing_fasta = site_dir / "target_and_observed_source_fcs.fasta"
    read1_fastq = site_dir / "selected.r1.fastq"
    read2_fastq = site_dir / "selected.r2.fastq"
    alignment_sam = site_dir / "target_source_competition.sam"
    minimap_log = site_dir / "target_source_competition.minimap2.log"
    contig_roles = write_competing_fasta(
        competing_fasta, run_root, target_fcs, source_fcs
    )
    command = [
        minimap2,
        "-ax",
        "sr",
        "--eqx",
        "--MD",
        "-F",
        "1000",
        "--end-bonus",
        "10",
        "-t",
        str(threads),
        str(competing_fasta),
        str(read1_fastq),
        str(read2_fastq),
    ]
    with alignment_sam.open("w") as output_handle, minimap_log.open("w") as log_handle:
        subprocess.run(command, check=True, stdout=output_handle, stderr=log_handle)

    forced_by_qname = {result.qname: result for result in forced_results}
    source_scores: dict[str, list[int]] = defaultdict(list)
    primary_roles: dict[str, set[str]] = defaultdict(set)
    with pysam.AlignmentFile(alignment_sam, "r") as alignment:
        for read in alignment.fetch(until_eof=True):
            qname = canonical_qname(read.query_name or "")
            forced = forced_by_qname.get(qname)
            if forced is None or read.is_unmapped or read.is_supplementary:
                continue
            mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
            if mate not in forced.marker_mates:
                continue
            role, _ = contig_roles[read.reference_name]
            if not read.is_secondary:
                primary_roles[qname].add(role)
            if role == "source" and read.has_tag("AS"):
                source_scores[qname].append(int(read.get_tag("AS")))

    results: dict[str, CompetitiveResult] = {}
    for qname, forced in forced_by_qname.items():
        target_scores = [
            record.alignment_score
            for record in forced.marker_records
            if record.alignment_score is not None
        ]
        results[qname] = CompetitiveResult(
            site_label=site.label,
            qname=qname,
            primary_roles=tuple(sorted(primary_roles[qname])),
            forced_target_score=max(target_scores) if target_scores else None,
            competing_source_score=max(source_scores[qname]) if source_scores[qname] else None,
        )
    return results


def write_marker_patterns(path: Path, markers: dict[str, str]) -> dict[str, str]:
    pattern_to_site: dict[str, str] = {}
    with path.open("w") as handle:
        for index, (site_label, sequence) in enumerate(sorted(markers.items()), start=1):
            pattern = f"marker_{index}"
            pattern_to_site[pattern] = site_label
            handle.write(f">{pattern}\n{sequence}\n")
    return pattern_to_site


def overlaps(interval: Interval, chrom: str, start: int, end: int) -> bool:
    return interval.chrom == chrom and interval.start < end and start < interval.end


def exact_marker_hits(
    output_dir: Path,
    markers: dict[str, str],
    target_fcs: dict[str, tuple[TargetFc, ...]],
    source_fcs: dict[str, tuple[TargetFc, ...]],
    reference: Path,
    seqkit: str,
    threads: int,
) -> dict[str, list[MarkerHit]]:
    patterns = output_dir / "alt_haplotype_markers.fasta"
    output = output_dir / "seqkit_exact_marker_locations.tsv"
    log = output_dir / "seqkit_exact_marker_locations.log"
    pattern_to_site = write_marker_patterns(patterns, markers)
    command = [
        seqkit,
        "locate",
        "--ignore-case",
        "--threads",
        str(threads),
        "--pattern-file",
        str(patterns),
        str(reference),
    ]
    with output.open("w") as output_handle, log.open("w") as log_handle:
        subprocess.run(command, check=True, stdout=output_handle, stderr=log_handle)

    hits: dict[str, list[MarkerHit]] = defaultdict(list)
    with output.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"seqID", "patternName", "strand", "start", "end"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"unexpected seqkit locate columns: {reader.fieldnames}")
        for row in reader:
            site_label = pattern_to_site[row["patternName"]]
            start = int(row["start"])
            end = int(row["end"])
            start0 = start - 1
            classifications: list[str] = []
            for fc in target_fcs[site_label]:
                if overlaps(fc.interval, row["seqID"], start0, end):
                    classifications.append(f"FN-site {fc.label}")
            for fc in source_fcs[site_label]:
                if overlaps(fc.interval, row["seqID"], start0, end):
                    classifications.append(f"observed-source {fc.label}")
            hits[site_label].append(
                MarkerHit(
                    site_label=site_label,
                    marker_length=len(markers[site_label]),
                    chrom=row["seqID"],
                    start=start,
                    end=end,
                    strand=row["strand"],
                    classifications=tuple(classifications),
                )
            )
    for site_label in markers:
        hits[site_label].sort(key=lambda hit: (hit.chrom, hit.start, hit.strand))
    return hits


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def replace_prefixed_lines(value: object, replacements: list[str]) -> str:
    prefixes = tuple(line.split(":", 1)[0] + ":" for line in replacements)
    retained = [
        line
        for line in str(value or "").splitlines()
        if line and not line.startswith(prefixes)
    ]
    return "\n".join(retained + replacements)


def save_workbook_atomic(workbook: openpyxl.Workbook, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    workbook.save(temporary)
    os.replace(temporary, output)


def update_workbook(
    input_path: Path,
    output_path: Path,
    target_fcs: dict[str, tuple[TargetFc, ...]],
    hits: dict[str, list[MarkerHit]],
    results_by_site: dict[str, list[QnameResult]],
    competition_by_site: dict[str, dict[str, CompetitiveResult]],
) -> None:
    workbook = openpyxl.load_workbook(input_path)
    if workbook.sheetnames != ["FN observations", "Qname mappings"]:
        raise ValueError(f"unexpected workbook sheets: {workbook.sheetnames}")
    observations = workbook["FN observations"]
    mappings = workbook["Qname mappings"]
    if worksheet_headers(observations) != OBSERVATION_COLUMNS:
        raise ValueError("FN observations schema differs from the agreed format")
    if worksheet_headers(mappings) != QNAME_COLUMNS:
        raise ValueError("Qname mappings schema differs from the agreed format")

    result_by_qname = {
        (result.site.label, result.qname): result
        for site_results in results_by_site.values()
        for result in site_results
    }
    for row_number, row in enumerate(
        observations.iter_rows(min_row=2, values_only=True), start=2
    ):
        if row[0] != "hg19" or row[1] not in results_by_site:
            continue
        site_label = str(row[1])
        site_results = results_by_site[site_label]
        competitions = competition_by_site[site_label]
        site_hits = hits[site_label]
        target_hit_count = sum(
            any(label.startswith("FN-site ") for label in hit.classifications)
            for hit in site_hits
        )
        source_hit_count = sum(
            any(label.startswith("observed-source ") for label in hit.classifications)
            for hit in site_hits
        )
        source_primary = sum(
            "source" in value.primary_roles for value in competitions.values()
        )
        source_higher = sum(
            value.score_delta_source_minus_target is not None
            and value.score_delta_source_minus_target > 0
            for value in competitions.values()
        )
        target_primary = sum(
            "target" in value.primary_roles for value in competitions.values()
        )
        if target_primary:
            cause = (
                "Cause: Cat4 | mixed marker-origin evidence: "
                f"{source_primary}/{len(site_results)} qnames are primary at the observed "
                f"non-target FC and {source_higher}/{len(site_results)} have a higher complete "
                f"marker-read alignment score there; {target_primary}/{len(site_results)} are "
                "primary at the FN-site FC. Forced target-only alignment demonstrates "
                "compatibility but does not establish target origin"
            )
        else:
            cause = (
                "Cause: Cat4 | marker-defined target origin is not established: "
                f"{source_primary}/{len(site_results)} qnames are primary at the observed "
                f"non-target FC and {source_higher}/{len(site_results)} have a higher complete "
                "marker-read alignment score there. Forced target-only alignment demonstrates "
                "compatibility but does not establish target origin"
            )
        replacements = [
            cause,
            "FN-site-overlapping FC: "
            + "; ".join(fc.display() for fc in target_fcs[site_label]),
            (
                f"ALT-haplotype marker exact hg19-reference matches: n={len(site_hits)}; "
                f"FN-site FC hits={target_hit_count}; observed non-target FC hits={source_hit_count}; "
                + (";".join(hit.display() for hit in site_hits) if site_hits else "none")
            ),
            (
                "Marker interpretation: the marker is not a repeated exact match; its sole "
                "exact hg19 reference match is at the observed non-target FC, so marker "
                "presence alone cannot prove FN-site origin"
            ),
            (
                "Forced target-FC remap (isolated target FC; production minimap2 arguments): "
                f"tested qnames={len(site_results)}; marker-bearing mate reaches FN="
                f"{sum(result.marker_reaches_fn for result in site_results)}; exact ALT CIGAR="
                f"{sum(result.exact_alt for result in site_results)}; both mates mapped="
                f"{sum(result.both_mates_mapped for result in site_results)}; proper pairs="
                f"{sum(result.proper_pair for result in site_results)}"
            ),
            (
                "Complete marker-read target/source competition: non-target primary="
                f"{source_primary}; non-target higher AS={source_higher}; "
                f"equal AS={sum(value.score_delta_source_minus_target == 0 for value in competitions.values())}; "
                f"target higher AS={sum(value.score_delta_source_minus_target is not None and value.score_delta_source_minus_target < 0 for value in competitions.values())}"
            ),
        ]
        observations.cell(
            row_number,
            11,
            replace_prefixed_lines(observations.cell(row_number, 11).value, replacements),
        )

    for row_number, row in enumerate(mappings.iter_rows(min_row=2, values_only=True), start=2):
        if row[0] != "hg19":
            continue
        key = (str(row[1]), canonical_qname(str(row[2])))
        result = result_by_qname.get(key)
        if result is None:
            continue
        replacements = [
            "FN-site-overlapping FC: "
            + "; ".join(fc.display() for fc in result.target_fcs),
            "Forced target-FC remap: " + result.assertion(),
            "Complete marker-read target/source competition: "
            + competition_by_site[result.site.label][result.qname].display(),
        ]
        mappings.cell(
            row_number,
            10,
            replace_prefixed_lines(mappings.cell(row_number, 10).value, replacements),
        )
    save_workbook_atomic(workbook, output_path)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cat4_dir = AUDIT_ROOT / "cat4_qname_audit"
    allele_rows = read_tsv(cat4_dir / "cat4_allele_summary.tsv")
    qname_rows = read_tsv(cat4_dir / "cat4_qname_summary.tsv")
    path_rows = read_tsv(cat4_dir / "cat4_qname_extraction_paths.tsv")
    marker_rows = read_tsv(HAPLOTYPE_DIR / "audited_alt_truth_kmers.tsv")
    sites, qnames_by_site = load_target_sites(allele_rows, qname_rows)
    target_fcs, source_fcs = load_site_fcs(path_rows, sites, args.run_root)
    markers = {
        row["FN_site"]: row["kmer_seq"].upper()
        for row in marker_rows
        if row["assembly"] == "hg19" and row["FN_site"] in sites
    }
    if set(markers) != set(sites):
        raise ValueError(f"marker/site mismatch: markers={sorted(markers)}, sites={sorted(sites)}")

    qnames_by_rg: dict[str, set[str]] = defaultdict(set)
    for label, qnames in qnames_by_site.items():
        rgs = {fc.rg for fc in target_fcs[label]}
        if len(rgs) != 1:
            raise ValueError(f"target FCs span multiple RGs for {label}: {sorted(rgs)}")
        qnames_by_rg[next(iter(rgs))].update(qnames)
    pairs = extract_fastq_pairs(args.run_root, qnames_by_rg)

    results_by_site: dict[str, list[QnameResult]] = {}
    competition_by_site: dict[str, dict[str, CompetitiveResult]] = {}
    for label, site in sorted(sites.items()):
        results_by_site[label] = run_forced_alignment(
            args.output_dir,
            site,
            target_fcs[label],
            qnames_by_site[label],
            pairs,
            markers[label],
            args.minimap2,
            args.threads,
            args.run_root,
        )
        competition_by_site[label] = run_competing_alignment(
            args.output_dir,
            site,
            target_fcs[label],
            source_fcs[label],
            results_by_site[label],
            args.minimap2,
            args.threads,
            args.run_root,
        )

    hits = exact_marker_hits(
        args.output_dir,
        markers,
        target_fcs,
        source_fcs,
        args.reference,
        args.seqkit,
        args.threads,
    )

    qname_output_rows: list[dict[str, object]] = []
    summary_output_rows: list[dict[str, object]] = []
    marker_output_rows: list[dict[str, object]] = []
    for label in sorted(sites):
        site_results = results_by_site[label]
        site_hits = hits[label]
        target_hit_count = sum(
            any(value.startswith("FN-site ") for value in hit.classifications)
            for hit in site_hits
        )
        source_hit_count = sum(
            any(value.startswith("observed-source ") for value in hit.classifications)
            for hit in site_hits
        )
        for result in site_results:
            competition = competition_by_site[label][result.qname]
            qname_output_rows.append(
                {
                    "Assembly": "hg19",
                    "FN allele": label,
                    "ALT-haplotype qname": result.qname,
                    "FN-site-overlapping FC": "; ".join(
                        fc.display() for fc in result.target_fcs
                    ),
                    "Marker length": len(markers[label]),
                    "Marker-bearing mate": ",".join(f"R{mate}" for mate in result.marker_mates),
                    "Exact marker occurrences in read pair": result.marker_occurrences_in_pair,
                    "Forced marker mate reaches FN": "yes" if result.marker_reaches_fn else "no",
                    "Forced marker mate produces exact ALT": "yes" if result.exact_alt else "no",
                    "Both mates mapped": "yes" if result.both_mates_mapped else "no",
                    "Proper pair": "yes" if result.proper_pair else "no",
                    "Competing primary locus": ",".join(competition.primary_roles),
                    "Forced target alignment score": competition.forced_target_score,
                    "Observed non-target alignment score": competition.competing_source_score,
                    "Non-target minus target alignment score": competition.score_delta_source_minus_target,
                    "Complete marker-read preference": competition.preference,
                    "Forced target-FC alignment": result.assertion(),
                }
            )
        summary_output_rows.append(
            {
                "Assembly": "hg19",
                "FN allele": label,
                "FN-site-overlapping FC": "; ".join(
                    fc.display() for fc in target_fcs[label]
                ),
                "Tested ALT-haplotype qnames": len(site_results),
                "Forced marker mate reaches FN": sum(
                    result.marker_reaches_fn for result in site_results
                ),
                "Forced marker mate produces exact ALT": sum(
                    result.exact_alt for result in site_results
                ),
                "Both mates mapped": sum(result.both_mates_mapped for result in site_results),
                "Proper pairs": sum(result.proper_pair for result in site_results),
                "Competing non-target primary": sum(
                    "source" in value.primary_roles
                    for value in competition_by_site[label].values()
                ),
                "Complete marker-read non-target higher score": sum(
                    value.score_delta_source_minus_target is not None
                    and value.score_delta_source_minus_target > 0
                    for value in competition_by_site[label].values()
                ),
                "Complete marker-read equal score": sum(
                    value.score_delta_source_minus_target == 0
                    for value in competition_by_site[label].values()
                ),
                "Complete marker-read target higher score": sum(
                    value.score_delta_source_minus_target is not None
                    and value.score_delta_source_minus_target < 0
                    for value in competition_by_site[label].values()
                ),
                "Marker length": len(markers[label]),
                "Exact marker matches in hg19 reference": len(site_hits),
                "Exact marker matches overlapping FN-site FC": target_hit_count,
                "Exact marker matches overlapping observed non-target FC": source_hit_count,
                "Exact marker match loci": ";".join(hit.display() for hit in site_hits),
            }
        )
        for hit in site_hits:
            marker_output_rows.append(
                {
                    "Assembly": "hg19",
                    "FN allele": label,
                    "Marker length": hit.marker_length,
                    "Chromosome": hit.chrom,
                    "Start": hit.start,
                    "End": hit.end,
                    "Strand": hit.strand,
                    "Locus classification": ";".join(hit.classifications) or "other locus",
                }
            )

    qname_output = args.output_dir / "cat4_forced_target_fc_qnames.tsv"
    summary_output = args.output_dir / "cat4_forced_target_fc_summary.tsv"
    marker_output = args.output_dir / "cat4_marker_exact_reference_matches.tsv"
    write_tsv(qname_output, qname_output_rows)
    write_tsv(summary_output, summary_output_rows)
    if marker_output_rows:
        write_tsv(marker_output, marker_output_rows)
    else:
        marker_output.write_text(
            "Assembly\tFN allele\tMarker length\tChromosome\tStart\tEnd\tStrand\tLocus classification\n"
        )

    total_qnames = sum(len(values) for values in results_by_site.values())
    total_reaching = sum(
        result.marker_reaches_fn
        for values in results_by_site.values()
        for result in values
    )
    total_exact = sum(
        result.exact_alt for values in results_by_site.values() for result in values
    )
    print(f"sites={len(sites)}")
    print(f"tested_qname_site_pairs={total_qnames}")
    print(f"forced_marker_mate_reaches_fn={total_reaching}")
    print(f"forced_marker_mate_exact_alt={total_exact}")
    print(f"qname_output={qname_output.resolve()}")
    print(f"summary_output={summary_output.resolve()}")
    print(f"marker_output={marker_output.resolve()}")

    if not args.update_workbook:
        print("workbook_updated=false")
        return
    output_workbook = args.output_workbook or args.workbook
    update_workbook(
        args.workbook,
        output_workbook,
        target_fcs,
        hits,
        results_by_site,
        competition_by_site,
    )
    print(f"workbook_output={output_workbook.resolve()}")
    print("workbook_updated=true")


if __name__ == "__main__":
    main()
