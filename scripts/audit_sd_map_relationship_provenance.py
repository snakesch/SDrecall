#!/usr/bin/env python3
"""Trace FN-site FC to observed remote-FC relationships through SD-map stages.

This audit is deliberately separate from ``FN_CAUSE_DIAGNOSIS.xlsx``.  It uses
the qname-placement evidence already recorded by ``audit_cat4_qname_paths.py``
to define the relationships that need explanation, then asks two questions:

1. At which reference-map preparation stages does a binary SD row directly
   connect the FN-containing FC and the remote FC reached by the ALT-haplotype
   read?
2. Are those FCs in the same connected component of the exact pruned multiplex
   graph built from the target-filtered map used by the production run?

The same candidate definition and overlap tests are used for hg19, hg38, and
CHM13.  An assembly with no observed remote-FC relationship still emits a
zero-count assembly summary instead of being silently skipped.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, TextIO


ROOT = Path("/paedyl01/disk1/yangyxt/SDrecall-rust-migration")
AUDIT_ROOT = ROOT / "test_tmp/nfc_fc_conflict_fix_20260723"
PUBLIC_WGAC = Path(
    "/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/hg19/test_BISER"
)
PUBLIC_CHM13 = Path("/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/t2t")
QNAME_PATHS = AUDIT_ROOT / "cat4_qname_audit/cat4_qname_extraction_paths.tsv"


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str = "."

    def overlaps(self, other: "Interval") -> bool:
        return (
            normalize_chrom(self.chrom) == normalize_chrom(other.chrom)
            and self.start < other.end
            and other.start < self.end
        )

    @property
    def size(self) -> int:
        return self.end - self.start

    def display(self) -> str:
        strand = "" if self.strand == "." else f"({self.strand})"
        return f"{normalize_chrom(self.chrom)}:{self.start}-{self.end}{strand}"

    def node_key(self) -> tuple[str, int, int, str]:
        return (normalize_chrom(self.chrom), self.start, self.end, self.strand)


@dataclass(frozen=True)
class PairFormat:
    chrom1: int = 0
    start1: int = 1
    end1: int = 2
    chrom2: int = 3
    start2: int = 4
    end2: int = 5
    strand1: int | None = 6
    strand2: int | None = 7

    @property
    def required_columns(self) -> int:
        indexes = [
            self.chrom1,
            self.start1,
            self.end1,
            self.chrom2,
            self.start2,
            self.end2,
        ]
        if self.strand1 is not None:
            indexes.append(self.strand1)
        if self.strand2 is not None:
            indexes.append(self.strand2)
        return max(indexes) + 1

    def parse(self, fields: list[str]) -> tuple[Interval, Interval] | None:
        if len(fields) < self.required_columns:
            return None
        try:
            a = Interval(
                fields[self.chrom1],
                int(fields[self.start1]),
                int(fields[self.end1]),
                fields[self.strand1] if self.strand1 is not None else ".",
            )
            b = Interval(
                fields[self.chrom2],
                int(fields[self.start2]),
                int(fields[self.end2]),
                fields[self.strand2] if self.strand2 is not None else ".",
            )
        except ValueError:
            return None
        if a.end <= a.start or b.end <= b.start:
            return None
        return a, b


STANDARD_PAIR = PairFormat()
BISER_RAW_PAIR = PairFormat(strand1=8, strand2=9)
TARGET_FILTERED_PAIR = PairFormat(
    chrom1=0,
    start1=1,
    end1=2,
    strand1=3,
    chrom2=4,
    start2=5,
    end2=6,
    strand2=7,
)


@dataclass(frozen=True)
class Stage:
    order: int
    name: str
    path: Path
    pair_format: PairFormat = STANDARD_PAIR
    producer: str = ""


@dataclass(frozen=True)
class AssemblyConfig:
    assembly: str
    run_root: Path
    prune_cutoff: float
    stages: tuple[Stage, ...]


RUN_ROOTS = {
    "hg19": AUDIT_ROOT / "runs/hg19/HG002_hg19_nfcfix23_SDrecall",
    "hg38": AUDIT_ROOT / "runs/hg38/HG002_hg38_nfcfix23_SDrecall",
    "chm13": AUDIT_ROOT / "runs/t2t/HG002_chm13_nfcfix23_SDrecall",
}


def standard_stages(assembly: str) -> tuple[Stage, ...]:
    prefix = f"WGAC.{assembly}"
    return (
        Stage(
            1,
            "caller_raw_pairs",
            PUBLIC_WGAC / f"{prefix}.bed",
            producer="WGAC pair-coordinate release",
        ),
        Stage(
            2,
            "cigar_generated_pairs",
            PUBLIC_WGAC / f"{prefix}.cigar.bed",
            producer=str(PUBLIC_WGAC / "generate_cigar_for_SD.py"),
        ),
        Stage(
            3,
            "cigar_segmented_pairs",
            PUBLIC_WGAC / f"{prefix}.cigar.trim.bed",
            producer="/paedyl01/disk1/yangyxt/ngs_scripts/condense_trim_coord_genbed.py trim",
        ),
        Stage(
            4,
            "divergence_filtered_pairs",
            PUBLIC_WGAC / f"{prefix}.cigar.trim.filtered.bed",
            producer="/paedyl01/disk1/yangyxt/ngs_scripts/refine_sd_coordinates.sh",
        ),
        Stage(
            5,
            "bidirectional_expanded_pairs",
            PUBLIC_WGAC / f"{prefix}.cigar.trim.homo.expanded.bed",
            producer="/paedyl01/disk1/yangyxt/ngs_scripts/refine_sd_coordinates.sh",
        ),
        Stage(
            6,
            "high_similarity_pairs",
            PUBLIC_WGAC / f"{prefix}.cigar.trim.homo.expanded.highsim.bed",
            producer="post-expansion high-similarity selection",
        ),
        Stage(
            7,
            "production_input_pairs",
            (
                ROOT / f"data/{assembly}/ref_SD/{prefix}.cigar.trim.homo.expanded.highsim.bed"
                if assembly == "hg19"
                else ROOT
                / "data/hg38/ref_SD/WGAC.hg38.cigar.trim.homo.expanded.highsim.with_alt.bed.gz"
            ),
            producer="production SDrecall -m input",
        ),
        Stage(
            8,
            "actual_target_filtered_pairs",
            RUN_ROOTS[assembly] / "realign_groups/filtered_SD_binary_map.tsv",
            pair_format=TARGET_FILTERED_PAIR,
            producer="sd-prep target-overlap + umbrella filter/dedup",
        ),
    )


CHM13_STAGES = (
    Stage(
        1,
        "caller_raw_pairs",
        PUBLIC_CHM13 / "BISER_SD.bed",
        pair_format=BISER_RAW_PAIR,
        producer="BISER pair-coordinate release",
    ),
    Stage(
        2,
        "normalized_8col_pairs",
        PUBLIC_CHM13 / "BISER_SD.8cols.bed",
        producer="BISER 8-column coordinate/strand normalization",
    ),
    Stage(
        3,
        "cigar_attempt_pairs",
        PUBLIC_CHM13 / "BISER_SD.8cols.cigar.bed",
        producer=str(ROOT / "benchmarks/process_refSD_map.sh"),
    ),
    Stage(
        4,
        "cigar_generated_pairs",
        PUBLIC_CHM13 / "BISER_SD.CHM13.cigar.bed",
        producer="successful CIGAR rows from benchmarks/process_refSD_map.sh",
    ),
    Stage(
        5,
        "cigar_segmented_pairs",
        PUBLIC_CHM13 / "BISER_SD.CHM13.cigar.trim.bed",
        producer="/paedyl01/disk1/yangyxt/ngs_scripts/condense_trim_coord_genbed.py trim",
    ),
    Stage(
        6,
        "divergence_filtered_pairs",
        PUBLIC_CHM13 / "BISER_SD.CHM13.cigar.trim.filtered.bed",
        producer="/paedyl01/disk1/yangyxt/ngs_scripts/refine_sd_coordinates.sh",
    ),
    Stage(
        7,
        "bidirectional_expanded_pairs",
        PUBLIC_CHM13 / "BISER_SD.CHM13.cigar.trim.expanded.bed",
        producer="/paedyl01/disk1/yangyxt/ngs_scripts/refine_sd_coordinates.sh",
    ),
    Stage(
        8,
        "production_input_pairs",
        ROOT / "data/chm13/ref_SD/BISER.chm13.cigar.trim.homo.expanded.highsim.bed",
        producer="production SDrecall -m input",
    ),
    Stage(
        9,
        "actual_target_filtered_pairs",
        RUN_ROOTS["chm13"] / "realign_groups/filtered_SD_binary_map.tsv",
        pair_format=TARGET_FILTERED_PAIR,
        producer="sd-prep target-overlap + umbrella filter/dedup",
    ),
)


CONFIGS = {
    "hg19": AssemblyConfig("hg19", RUN_ROOTS["hg19"], 423.17, standard_stages("hg19")),
    "hg38": AssemblyConfig("hg38", RUN_ROOTS["hg38"], 412.50, standard_stages("hg38")),
    "chm13": AssemblyConfig("chm13", RUN_ROOTS["chm13"], 406.77, CHM13_STAGES),
}


@dataclass
class Candidate:
    assembly: str
    allele: str
    rg: str
    target_fc: str
    observed_fc: str
    target_intervals: tuple[Interval, ...]
    observed_intervals: tuple[Interval, ...]
    supporting_qnames: set[str] = field(default_factory=set)
    target_reaching_qnames: set[str] = field(default_factory=set)

    @property
    def key(self) -> tuple[str, str, str, str, str]:
        return (self.assembly, self.allele, self.rg, self.target_fc, self.observed_fc)

    @property
    def remote_only_qnames(self) -> set[str]:
        return self.supporting_qnames - self.target_reaching_qnames

    @property
    def causal_role(self) -> str:
        if not self.target_reaching_qnames:
            return "exclusive_remote_placement"
        if self.remote_only_qnames:
            return "mixed_remote_and_target_placement"
        return "incidental_remote_placement_target_also_reached"


@dataclass
class StageResult:
    candidate: Candidate
    stage: Stage
    rows_read: int
    target_locus_rows: int
    observed_locus_rows: int
    direct_rows: int
    direct_examples: tuple[str, ...]


@dataclass(frozen=True)
class ComponentResult:
    connected: bool | None
    target_anchor_mode: str
    target_anchor_nodes: int
    observed_anchor_mode: str
    observed_anchor_nodes: int
    graph_nodes: int
    graph_edges: int
    graph_sd_edges: int
    graph_po_edges_kept: int


@dataclass
class ReconstructedEdge:
    """The two independent edge flags retained by the Rust multiplex graph."""

    is_sd: bool = False
    is_overlap: bool = False
    overlap_bp: int = 0
    overlap_fraction: float = 0.0


class DisjointSet:
    def __init__(self) -> None:
        self.parent: list[int] = []
        self.rank: list[int] = []

    def add(self) -> int:
        idx = len(self.parent)
        self.parent.append(idx)
        self.rank.append(0)
        return idx

    def find(self, value: int) -> int:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: int, right: int) -> bool:
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return False
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1
        return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True, choices=tuple(CONFIGS))
    parser.add_argument("--qname-paths", type=Path, default=QNAME_PATHS)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def split_labels(value: str) -> list[str]:
    return [label for label in value.split(";") if label]


def display_intervals(intervals: Iterable[Interval]) -> str:
    return ";".join(interval.display() for interval in intervals)


def read_fc_catalog(run_root: Path) -> dict[str, tuple[Interval, ...]]:
    catalog: dict[str, list[Interval]] = defaultdict(list)
    for path in sorted((run_root / "realign_groups").glob("RG*/RG*_related_homo_regions.bed")):
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or not fields[6].startswith("FC:"):
                    continue
                catalog[fields[6]].append(
                    Interval(fields[0], int(fields[1]), int(fields[2]), fields[5])
                )
    return {label: tuple(intervals) for label, intervals in catalog.items()}


def intervals_are_remote(left: Iterable[Interval], right: Iterable[Interval]) -> bool:
    return all(not a.overlaps(b) for a in left for b in right)


def load_candidates(
    assembly: str, qname_paths: Path, fc_catalog: dict[str, tuple[Interval, ...]]
) -> list[Candidate]:
    with qname_paths.open() as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["Assembly"] == assembly
        ]

    reaches_target: dict[tuple[str, str], bool] = defaultdict(bool)
    for row in rows:
        key = (row["FN allele"], row["ALT-haplotype qname"])
        if row["ALT-haplotype record overlaps FN anchor"] == "yes":
            reaches_target[key] = True

    candidates: dict[tuple[str, str, str, str, str], Candidate] = {}
    for row in rows:
        targets = split_labels(row["All FN-containing FCs in RG"])
        observed = split_labels(row["FCs overlapped by ALT-haplotype record"])
        for target_fc in targets:
            for observed_fc in observed:
                if target_fc == observed_fc:
                    continue
                if target_fc not in fc_catalog or observed_fc not in fc_catalog:
                    raise ValueError(
                        f"missing FC catalog entry: {target_fc!r} or {observed_fc!r}"
                    )
                target_intervals = fc_catalog[target_fc]
                observed_intervals = fc_catalog[observed_fc]
                if not intervals_are_remote(target_intervals, observed_intervals):
                    continue
                key = (assembly, row["FN allele"], row["RG"], target_fc, observed_fc)
                candidate = candidates.setdefault(
                    key,
                    Candidate(
                        assembly=assembly,
                        allele=row["FN allele"],
                        rg=row["RG"],
                        target_fc=target_fc,
                        observed_fc=observed_fc,
                        target_intervals=target_intervals,
                        observed_intervals=observed_intervals,
                    ),
                )
                qname = row["ALT-haplotype qname"]
                candidate.supporting_qnames.add(qname)
                if reaches_target[(row["FN allele"], qname)]:
                    candidate.target_reaching_qnames.add(qname)
    return [candidates[key] for key in sorted(candidates)]


def any_overlap(interval: Interval, targets: Iterable[Interval]) -> bool:
    return any(interval.overlaps(target) for target in targets)


def pair_display(a: Interval, b: Interval) -> str:
    return f"{a.display()}<->{b.display()}"


def scan_stage(stage: Stage, candidates: list[Candidate]) -> list[StageResult]:
    if not stage.path.is_file():
        raise FileNotFoundError(stage.path)
    if not candidates:
        return []
    target_rows = [0] * len(candidates)
    observed_rows = [0] * len(candidates)
    direct_rows = [0] * len(candidates)
    examples: list[list[str]] = [[] for _ in candidates]
    rows_read = 0
    with open_text(stage.path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parsed = stage.pair_format.parse(line.rstrip("\n").split("\t"))
            if parsed is None:
                continue
            rows_read += 1
            a, b = parsed
            for idx, candidate in enumerate(candidates):
                a_target = any_overlap(a, candidate.target_intervals)
                b_target = any_overlap(b, candidate.target_intervals)
                a_observed = any_overlap(a, candidate.observed_intervals)
                b_observed = any_overlap(b, candidate.observed_intervals)
                if a_target or b_target:
                    target_rows[idx] += 1
                if a_observed or b_observed:
                    observed_rows[idx] += 1
                if (a_target and b_observed) or (b_target and a_observed):
                    direct_rows[idx] += 1
                    if len(examples[idx]) < 3:
                        examples[idx].append(pair_display(a, b))
    return [
        StageResult(
            candidate=candidate,
            stage=stage,
            rows_read=rows_read,
            target_locus_rows=target_rows[idx],
            observed_locus_rows=observed_rows[idx],
            direct_rows=direct_rows[idx],
            direct_examples=tuple(examples[idx]),
        )
        for idx, candidate in enumerate(candidates)
    ]


def component_anchor_nodes(
    intervals: tuple[Interval, ...],
    node_lookup: dict[tuple[str, int, int, str], int],
) -> tuple[str, set[int]]:
    exact = {
        node_lookup[interval.node_key()]
        for interval in intervals
        if interval.node_key() in node_lookup
    }
    if exact:
        return "exact_FC_node", exact
    return "absent_after_pruning", set()


def build_actual_component_results(
    config: AssemblyConfig, candidates: list[Candidate]
) -> dict[tuple[str, str, str, str, str], ComponentResult]:
    if not candidates:
        return {}
    stage = next(s for s in config.stages if s.name == "actual_target_filtered_pairs")
    rows: list[tuple[Interval, Interval]] = []
    nodes_by_key: dict[tuple[str, int, int, str], Interval] = {}
    first_seen: dict[tuple[str, int, int, str], int] = {}
    chroms: list[str] = []
    seen_chroms: set[str] = set()
    with open_text(stage.path) as handle:
        for line in handle:
            parsed = stage.pair_format.parse(line.rstrip("\n").split("\t"))
            if parsed is None:
                continue
            a, b = parsed
            rows.append((a, b))
            for interval in (a, b):
                key = interval.node_key()
                nodes_by_key.setdefault(key, interval)
                first_seen.setdefault(key, len(first_seen))
                chrom = normalize_chrom(interval.chrom)
                if chrom not in seen_chroms:
                    seen_chroms.add(chrom)
                    chroms.append(chrom)

    # Reproduce graph_build::compose_po_per_chr. In particular, equal-sized
    # overlapping nodes can acquire both directed PO edges after repeated row
    # appearances, and that direction controls whether an SD edge is combined.
    edges: dict[
        tuple[tuple[str, int, int, str], tuple[str, int, int, str]],
        ReconstructedEdge,
    ] = {}
    po_edge_seen: set[
        tuple[tuple[str, int, int, str], tuple[str, int, int, str]]
    ] = set()
    for chrom in chroms:
        inserted: list[Interval] = []
        for a, b in rows:
            if normalize_chrom(a.chrom) != chrom and normalize_chrom(b.chrom) != chrom:
                continue
            for current in (a, b):
                for prior in inserted:
                    if current.node_key() == prior.node_key():
                        continue
                    overlap = min(prior.end, current.end) - max(
                        prior.start, current.start
                    )
                    if normalize_chrom(current.chrom) != normalize_chrom(prior.chrom):
                        continue
                    if overlap <= 0:
                        continue
                    if current.size >= prior.size:
                        large, small = current, prior
                    else:
                        large, small = prior, current
                    edge_key = (large.node_key(), small.node_key())
                    if edge_key in po_edge_seen:
                        continue
                    po_edge_seen.add(edge_key)
                    edges[edge_key] = ReconstructedEdge(
                        is_overlap=True,
                        overlap_bp=overlap,
                        overlap_fraction=overlap / min(large.size, small.size),
                    )
                if normalize_chrom(current.chrom) == chrom:
                    inserted.append(current)

    # Overlay the undirected SD pairs in first-node-appearance orientation.
    # If the same directed PO edge exists, Rust changes only the SD flag and
    # retains the PO attributes used by the later pruning pass.
    sd_pair_seen: set[frozenset[tuple[str, int, int, str]]] = set()
    for a, b in rows:
        a_key = a.node_key()
        b_key = b.node_key()
        unordered = frozenset((a_key, b_key))
        if unordered in sd_pair_seen:
            continue
        sd_pair_seen.add(unordered)
        if first_seen[a_key] <= first_seen[b_key]:
            edge_key = (a_key, b_key)
        else:
            edge_key = (b_key, a_key)
        edge = edges.get(edge_key)
        if edge is None:
            edges[edge_key] = ReconstructedEdge(is_sd=True)
        else:
            edge.is_sd = True

    # Reproduce traversal::prune_graph, then collapse surviving directed edges
    # only for connected-component membership (production also labels the pruned
    # graph as undirected at this point).
    dsu = DisjointSet()
    node_lookup: dict[tuple[str, int, int, str], int] = {}
    for key, interval in nodes_by_key.items():
        if interval.size <= config.prune_cutoff:
            continue
        node_lookup[key] = dsu.add()

    graph_edges = 0
    sd_edges = 0
    po_edges_kept = 0
    for (left, right), edge in edges.items():
        if left == right or left not in node_lookup or right not in node_lookup:
            continue
        if (
            edge.is_overlap
            and edge.overlap_bp < config.prune_cutoff
            and edge.overlap_fraction < 0.5
        ):
            continue
        graph_edges += 1
        sd_edges += int(edge.is_sd)
        po_edges_kept += int(edge.is_overlap)
        dsu.union(node_lookup[left], node_lookup[right])

    results: dict[tuple[str, str, str, str, str], ComponentResult] = {}
    for candidate in candidates:
        target_mode, target_nodes = component_anchor_nodes(
            candidate.target_intervals, node_lookup
        )
        observed_mode, observed_nodes = component_anchor_nodes(
            candidate.observed_intervals, node_lookup
        )
        target_roots = {dsu.find(idx) for idx in target_nodes}
        observed_roots = {dsu.find(idx) for idx in observed_nodes}
        results[candidate.key] = ComponentResult(
            connected=(
                bool(target_roots & observed_roots)
                if target_roots and observed_roots
                else None
            ),
            target_anchor_mode=target_mode,
            target_anchor_nodes=len(target_nodes),
            observed_anchor_mode=observed_mode,
            observed_anchor_nodes=len(observed_nodes),
            graph_nodes=len(node_lookup),
            graph_edges=graph_edges,
            graph_sd_edges=sd_edges,
            graph_po_edges_kept=po_edges_kept,
        )
    return results


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def first_direct_loss(results: list[StageResult]) -> str:
    if not results:
        return "no_stages"
    ordered = sorted(results, key=lambda item: item.stage.order)
    if ordered[0].direct_rows == 0:
        return "absent_from_caller_raw_pairs"
    previous = ordered[0].stage.name
    for result in ordered[1:]:
        present = result.direct_rows > 0
        if present:
            previous = result.stage.name
            continue
        return f"after_{previous}__before_{result.stage.name}"
    return "direct_relationship_survives_all_stages"


def main() -> None:
    args = parse_args()
    config = CONFIGS[args.assembly]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fc_catalog = read_fc_catalog(config.run_root)
    candidates = load_candidates(args.assembly, args.qname_paths, fc_catalog)
    stage_results: list[StageResult] = []
    for stage in config.stages:
        stage_results.extend(scan_stage(stage, candidates))
    components = build_actual_component_results(config, candidates)

    write_tsv(
        args.output_dir / "relationship_candidates.tsv",
        [
            "Assembly",
            "FN allele",
            "RG",
            "FN-site FC",
            "FN-site FC interval",
            "Observed remote FC",
            "Observed remote FC interval",
            "Supporting ALT-haplotype qnames",
            "Qnames also reaching FN anchor",
            "Remote-only qnames",
            "Relationship causal role",
        ],
        (
            {
                "Assembly": c.assembly,
                "FN allele": c.allele,
                "RG": c.rg,
                "FN-site FC": c.target_fc,
                "FN-site FC interval": display_intervals(c.target_intervals),
                "Observed remote FC": c.observed_fc,
                "Observed remote FC interval": display_intervals(c.observed_intervals),
                "Supporting ALT-haplotype qnames": len(c.supporting_qnames),
                "Qnames also reaching FN anchor": len(c.target_reaching_qnames),
                "Remote-only qnames": len(c.remote_only_qnames),
                "Relationship causal role": c.causal_role,
            }
            for c in candidates
        ),
    )

    write_tsv(
        args.output_dir / "stage_audit.tsv",
        [
            "Assembly",
            "FN allele",
            "RG",
            "FN-site FC",
            "Observed remote FC",
            "Stage order",
            "Stage",
            "Stage file",
            "Stage rows read",
            "Rows touching FN-site FC",
            "Rows touching observed remote FC",
            "Direct relationship rows",
            "Direct relationship examples",
        ],
        (
            {
                "Assembly": r.candidate.assembly,
                "FN allele": r.candidate.allele,
                "RG": r.candidate.rg,
                "FN-site FC": r.candidate.target_fc,
                "Observed remote FC": r.candidate.observed_fc,
                "Stage order": r.stage.order,
                "Stage": r.stage.name,
                "Stage file": str(r.stage.path),
                "Stage rows read": r.rows_read,
                "Rows touching FN-site FC": r.target_locus_rows,
                "Rows touching observed remote FC": r.observed_locus_rows,
                "Direct relationship rows": r.direct_rows,
                "Direct relationship examples": ";".join(r.direct_examples),
            }
            for r in stage_results
        ),
    )

    results_by_candidate: dict[tuple[str, str, str, str, str], list[StageResult]] = defaultdict(list)
    for result in stage_results:
        results_by_candidate[result.candidate.key].append(result)

    summary_rows: list[dict[str, object]] = []
    for candidate in candidates:
        results = sorted(results_by_candidate[candidate.key], key=lambda item: item.stage.order)
        component = components[candidate.key]
        summary_rows.append(
            {
                "Assembly": candidate.assembly,
                "FN allele": candidate.allele,
                "RG": candidate.rg,
                "FN-site FC": candidate.target_fc,
                "Observed remote FC": candidate.observed_fc,
                "Relationship causal role": candidate.causal_role,
                "Direct relationship stage path": ";".join(
                    f"{r.stage.name}={r.direct_rows}" for r in results
                ),
                "First direct relationship loss": first_direct_loss(results),
                "Same actual pruned multiplex component": (
                    "yes"
                    if component.connected is True
                    else "no"
                    if component.connected is False
                    else "not_evaluable"
                ),
                "FN-site FC graph anchor": f"{component.target_anchor_mode}:n={component.target_anchor_nodes}",
                "Observed FC graph anchor": f"{component.observed_anchor_mode}:n={component.observed_anchor_nodes}",
                "Production prune cutoff bp": f"{config.prune_cutoff:.2f}",
                "Actual graph nodes": component.graph_nodes,
                "Actual graph edges": component.graph_edges,
                "Actual graph SD edges": component.graph_sd_edges,
                "Actual graph kept PO edges": component.graph_po_edges_kept,
            }
        )
    write_tsv(
        args.output_dir / "relationship_summary.tsv",
        [
            "Assembly",
            "FN allele",
            "RG",
            "FN-site FC",
            "Observed remote FC",
            "Relationship causal role",
            "Direct relationship stage path",
            "First direct relationship loss",
            "Same actual pruned multiplex component",
            "FN-site FC graph anchor",
            "Observed FC graph anchor",
            "Production prune cutoff bp",
            "Actual graph nodes",
            "Actual graph edges",
            "Actual graph SD edges",
            "Actual graph kept PO edges",
        ],
        summary_rows,
    )

    write_tsv(
        args.output_dir / "stage_provenance.tsv",
        ["Assembly", "Stage order", "Stage", "Stage file", "Producer"],
        (
            {
                "Assembly": config.assembly,
                "Stage order": stage.order,
                "Stage": stage.name,
                "Stage file": str(stage.path),
                "Producer": stage.producer,
            }
            for stage in config.stages
        ),
    )

    write_tsv(
        args.output_dir / "assembly_summary.tsv",
        [
            "Assembly",
            "Observed remote-FC relationships",
            "FN alleles with remote-FC relationships",
            "Stages audited",
            "Production prune cutoff bp",
            "Result",
        ],
        [
            {
                "Assembly": config.assembly,
                "Observed remote-FC relationships": len(candidates),
                "FN alleles with remote-FC relationships": len(
                    {candidate.allele for candidate in candidates}
                ),
                "Stages audited": len(config.stages),
                "Production prune cutoff bp": f"{config.prune_cutoff:.2f}",
                "Result": (
                    "audited"
                    if candidates
                    else "audited_no_observed_remote_FC_relationships"
                ),
            }
        ],
    )

    print(
        f"{config.assembly}: candidates={len(candidates)} "
        f"alleles={len({candidate.allele for candidate in candidates})} "
        f"stages={len(config.stages)} output={args.output_dir}"
    )


if __name__ == "__main__":
    main()
