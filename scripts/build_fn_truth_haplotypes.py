#!/usr/bin/env python3
"""Build ALT-only local truth haplotypes for false-negative tracing.

The generated markers always contain the complete target ALT allele. Phased
truth genotypes are projected directly; unphased heterozygous variants are
phased with PacBio HiFi reads when needed. Each marker also carries a per-base
reference-coordinate map so downstream tracing can verify SNVs and indels.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class Site:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int | None, ...]
    phased: bool

    @property
    def label(self) -> str:
        separator = "|" if self.phased else "/"
        gt = separator.join("." if allele is None else str(allele) for allele in self.gt)
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}:{gt}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--hifi-bam", required=True, type=Path)
    parser.add_argument("--output-kmers", required=True, type=Path)
    parser.add_argument("--output-definitions", required=True, type=Path)
    parser.add_argument("--kmer-half", type=int, default=35)
    parser.add_argument("--min-base-quality", type=int, default=20)
    parser.add_argument("--min-phase-reads", type=int, default=2)
    parser.add_argument("--phase-ratio", type=float, default=0.70)
    return parser.parse_args()


def read_sites(path: Path) -> list[Site]:
    sites: list[Site] = []
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or row[0] == "assembly":
                continue
            if len(row) < 5:
                raise ValueError(f"Expected five columns in {path}: {row}")
            sites.append(Site(row[0], row[1], int(row[2]), row[3], row[4]))
    return sites


def sample_name(vcf: pysam.VariantFile) -> str:
    samples = list(vcf.header.samples)
    if len(samples) != 1:
        raise ValueError(f"Expected one truth sample, found {samples}")
    return samples[0]


def load_window_variants(
    vcf: pysam.VariantFile,
    sample: str,
    site: Site,
    start0: int,
    end0: int,
) -> list[Variant]:
    variants: list[Variant] = []
    for record in vcf.fetch(site.chrom, start0, end0):
        if len(record.alts or ()) != 1:
            raise ValueError(f"Expected biallelic normalized truth record: {record}")
        sample_data = record.samples[sample]
        gt = tuple(sample_data.get("GT") or ())
        variants.append(
            Variant(
                chrom=record.chrom,
                pos=record.pos,
                ref=record.ref.upper(),
                alt=record.alts[0].upper(),
                gt=gt,
                phased=sample_data.phased,
            )
        )
    return variants


def base_call(
    read: pysam.AlignedSegment,
    ref_pos0: int,
    min_base_quality: int,
) -> str | None:
    sequence = read.query_sequence
    qualities = read.query_qualities
    if sequence is None or qualities is None:
        return None
    for query_pos, aligned_ref_pos in read.get_aligned_pairs(matches_only=False):
        if aligned_ref_pos != ref_pos0:
            continue
        if query_pos is None or qualities[query_pos] < min_base_quality:
            return None
        return sequence[query_pos].upper()
    return None


def allele_call(
    read: pysam.AlignedSegment,
    variant: Variant,
    min_base_quality: int,
) -> str | None:
    """Return ``ref`` or ``alt`` for a simple normalized allele on one read."""
    if len(variant.ref) == len(variant.alt):
        observed = "".join(
            base_call(read, variant.pos - 1 + offset, min_base_quality) or "?"
            for offset in range(len(variant.ref))
        )
        if observed == variant.ref:
            return "ref"
        if observed == variant.alt:
            return "alt"
        return None

    pairs = read.get_aligned_pairs(matches_only=False)
    sequence = read.query_sequence or ""
    qualities = read.query_qualities
    if qualities is None:
        return None

    anchor_ref0 = variant.pos - 1
    anchor_index = next(
        (
            index
            for index, (query_pos, ref_pos) in enumerate(pairs)
            if query_pos is not None and ref_pos == anchor_ref0
        ),
        None,
    )
    if anchor_index is None:
        return None
    anchor_query = pairs[anchor_index][0]
    if anchor_query is None or qualities[anchor_query] < min_base_quality:
        return None
    if sequence[anchor_query].upper() != variant.ref[0]:
        return None

    # Normalized simple insertion: the ALT starts with the REF anchor.
    if len(variant.ref) == 1 and variant.alt.startswith(variant.ref):
        inserted_query_positions: list[int] = []
        for query_pos, ref_pos in pairs[anchor_index + 1 :]:
            if ref_pos is not None:
                break
            if query_pos is not None:
                inserted_query_positions.append(query_pos)
        if any(qualities[pos] < min_base_quality for pos in inserted_query_positions):
            return None
        inserted = "".join(sequence[pos].upper() for pos in inserted_query_positions)
        if inserted == variant.alt[1:]:
            return "alt"
        if not inserted:
            return "ref"
        return None

    # Normalized simple deletion: the ALT retains only the left anchor.
    if len(variant.alt) == 1 and variant.ref.startswith(variant.alt):
        deleted_ref_positions: list[int] = []
        for query_pos, ref_pos in pairs[anchor_index + 1 :]:
            if query_pos is not None:
                break
            if ref_pos is not None:
                deleted_ref_positions.append(ref_pos)
        expected_deleted = list(range(anchor_ref0 + 1, anchor_ref0 + len(variant.ref)))
        if deleted_ref_positions == expected_deleted:
            return "alt"

        observed_ref = "".join(
            base_call(read, anchor_ref0 + offset, min_base_quality) or "?"
            for offset in range(len(variant.ref))
        )
        if observed_ref == variant.ref:
            return "ref"
        return None

    return None


def resolve_bam_contig(bam: pysam.AlignmentFile, chrom: str) -> str:
    if chrom in bam.references:
        return chrom
    aliases = [chrom.removeprefix("chr")]
    if not chrom.startswith("chr"):
        aliases.append(f"chr{chrom}")
    if chrom in {"chrM", "M"}:
        aliases.extend(["MT", "chrMT"])
    for alias in aliases:
        if alias in bam.references:
            return alias
    raise ValueError(f"Contig {chrom!r} is absent from HiFi BAM")


def infer_same_haplotype(
    bam: pysam.AlignmentFile,
    target: Variant,
    adjacent: Variant,
    min_base_quality: int,
    min_phase_reads: int,
    phase_ratio: float,
) -> tuple[bool, str]:
    start0 = min(target.pos, adjacent.pos) - 1
    end0 = max(target.pos, adjacent.pos)
    bam_chrom = resolve_bam_contig(bam, target.chrom)
    patterns: Counter[tuple[str, str]] = Counter()
    seen_qnames: set[str] = set()

    for read in bam.fetch(bam_chrom, start0, end0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_name in seen_qnames:
            continue
        target_call = allele_call(read, target, min_base_quality)
        adjacent_call = allele_call(read, adjacent, min_base_quality)
        if target_call is None or adjacent_call is None:
            continue
        seen_qnames.add(read.query_name)
        patterns[(target_call, adjacent_call)] += 1

    target_alt_adj_alt = patterns[("alt", "alt")]
    target_alt_adj_ref = patterns[("alt", "ref")]
    informative = target_alt_adj_alt + target_alt_adj_ref
    evidence = ",".join(
        f"{first}/{second}={count}"
        for (first, second), count in sorted(patterns.items())
    ) or "none"

    if informative < min_phase_reads:
        raise ValueError(
            f"Insufficient HiFi phase evidence for {target.label} and {adjacent.label}: "
            f"target-ALT informative reads={informative}; patterns={evidence}"
        )

    same_fraction = target_alt_adj_alt / informative
    opposite_fraction = target_alt_adj_ref / informative
    if same_fraction >= phase_ratio:
        return True, f"same_haplotype;target_ALT_reads={informative};patterns={evidence}"
    if opposite_fraction >= phase_ratio:
        return False, f"opposite_haplotype;target_ALT_reads={informative};patterns={evidence}"
    raise ValueError(
        f"Ambiguous HiFi phase for {target.label} and {adjacent.label}: "
        f"same={same_fraction:.3f}, opposite={opposite_fraction:.3f}; patterns={evidence}"
    )


def alt_reference_map(variant: Variant) -> list[int | None]:
    """Map ALT bases to 1-based reference positions; insertions map to None."""
    if len(variant.ref) == len(variant.alt):
        return [variant.pos + offset for offset in range(len(variant.alt))]

    prefix = 0
    while (
        prefix < len(variant.ref)
        and prefix < len(variant.alt)
        and variant.ref[prefix] == variant.alt[prefix]
    ):
        prefix += 1

    suffix = 0
    while (
        suffix < len(variant.ref) - prefix
        and suffix < len(variant.alt) - prefix
        and variant.ref[-1 - suffix] == variant.alt[-1 - suffix]
    ):
        suffix += 1

    mapping: list[int | None] = [variant.pos + offset for offset in range(prefix)]
    ref_core_len = len(variant.ref) - prefix - suffix
    alt_core_len = len(variant.alt) - prefix - suffix
    if ref_core_len == alt_core_len:
        mapping.extend(variant.pos + prefix + offset for offset in range(alt_core_len))
    else:
        mapping.extend([None] * alt_core_len)
    if suffix:
        suffix_start = variant.pos + len(variant.ref) - suffix
        mapping.extend(suffix_start + offset for offset in range(suffix))
    return mapping


def apply_variants(
    reference_sequence: str,
    window_start0: int,
    variants: list[Variant],
) -> tuple[str, list[int | None]]:
    sequence = reference_sequence
    reference_map: list[int | None] = [
        window_start0 + offset + 1 for offset in range(len(reference_sequence))
    ]
    for variant in sorted(variants, key=lambda item: item.pos, reverse=True):
        offset = variant.pos - 1 - window_start0
        observed_ref = sequence[offset : offset + len(variant.ref)]
        if observed_ref != variant.ref:
            raise ValueError(
                f"Reference mismatch at {variant.label}: expected {variant.ref}, "
                f"observed {observed_ref}"
            )
        sequence = sequence[:offset] + variant.alt + sequence[offset + len(variant.ref) :]
        reference_map = (
            reference_map[:offset]
            + alt_reference_map(variant)
            + reference_map[offset + len(variant.ref) :]
        )
    return sequence, reference_map


def contains_alt_allele(
    sequence: str,
    reference_map: list[int | None],
    target: Variant,
) -> bool:
    """Return whether the target ALT and its coordinate pattern occur together."""
    target_map = alt_reference_map(target)
    allele_length = len(target.alt)
    return any(
        sequence[offset : offset + allele_length] == target.alt
        and reference_map[offset : offset + allele_length] == target_map
        for offset in range(len(sequence) - allele_length + 1)
    )


def format_reference_map(reference_map: list[int | None]) -> str:
    return ",".join("." if position is None else str(position) for position in reference_map)


def called_alleles(variant: Variant) -> list[int]:
    return [allele for allele in variant.gt if allele is not None]


def phased_variant_sets(
    variants: list[Variant],
    target: Variant,
) -> list[tuple[int, list[Variant], list[Variant], list[str]]]:
    target_haplotypes = [index for index, allele in enumerate(target.gt) if allele == 1]
    if not target_haplotypes:
        raise ValueError(f"Target truth genotype does not contain ALT: {target.label}")

    results: list[tuple[int, list[Variant], list[Variant], list[str]]] = []
    for haplotype in target_haplotypes:
        applied: list[Variant] = []
        excluded: list[Variant] = []
        for variant in variants:
            if not variant.phased or len(variant.gt) <= haplotype:
                alleles = called_alleles(variant)
                if alleles and all(allele == 1 for allele in alleles):
                    applied.append(variant)
                    continue
                raise ValueError(
                    f"Cannot project unphased variant {variant.label} onto phased target "
                    f"haplotype {haplotype} for {target.label}"
                )
            allele = variant.gt[haplotype]
            if allele == 1:
                applied.append(variant)
            elif allele == 0:
                excluded.append(variant)
            else:
                raise ValueError(
                    f"Unsupported allele {allele} on haplotype {haplotype} at {variant.label}"
                )
        results.append(
            (
                haplotype,
                applied,
                excluded,
                [f"truth_phased_haplotype={haplotype}"],
            )
        )
    return results


def unphased_variant_sets(
    bam: pysam.AlignmentFile,
    variants: list[Variant],
    target: Variant,
    min_base_quality: int,
    min_phase_reads: int,
    phase_ratio: float,
) -> list[tuple[int, list[Variant], list[Variant], list[str]]]:
    target_alleles = called_alleles(target)
    if sorted(target_alleles) == [0, 1]:
        applied = [target]
        excluded: list[Variant] = []
        phase_notes: list[str] = []
        for variant in variants:
            if variant == target:
                continue
            alleles = called_alleles(variant)
            if alleles and all(allele == 1 for allele in alleles):
                applied.append(variant)
                continue
            if alleles and all(allele == 0 for allele in alleles):
                excluded.append(variant)
                continue
            if sorted(alleles) == [0, 1]:
                same_haplotype, evidence = infer_same_haplotype(
                    bam,
                    target,
                    variant,
                    min_base_quality,
                    min_phase_reads,
                    phase_ratio,
                )
                phase_notes.append(f"{variant.label}:{evidence}")
                if same_haplotype:
                    applied.append(variant)
                else:
                    excluded.append(variant)
                continue
            raise ValueError(f"Unsupported adjacent truth genotype at {variant.label}")
        return [(0, applied, excluded, phase_notes or ["not_needed"])]

    if target_alleles and all(allele == 1 for allele in target_alleles):
        applied = []
        excluded = []
        for variant in variants:
            alleles = called_alleles(variant)
            if alleles and all(allele == 1 for allele in alleles):
                applied.append(variant)
            elif alleles and all(allele == 0 for allele in alleles):
                excluded.append(variant)
            else:
                raise ValueError(
                    f"Unphased homozygous-ALT target {target.label} has a heterozygous "
                    f"adjacent variant {variant.label}; two local haplotypes are required"
                )
        return [(0, applied, excluded, ["unphased_homozygous_ALT"])]

    raise ValueError(f"Unsupported target truth genotype at {target.label}")


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)
    fasta = pysam.FastaFile(str(args.reference))
    truth = pysam.VariantFile(str(args.truth_vcf))
    hifi = pysam.AlignmentFile(str(args.hifi_bam), "rb")
    truth_sample = sample_name(truth)

    args.output_kmers.parent.mkdir(parents=True, exist_ok=True)
    args.output_definitions.parent.mkdir(parents=True, exist_ok=True)

    with args.output_kmers.open("w") as kmer_handle, args.output_definitions.open("w") as definition_handle:
        kmer_writer = csv.writer(kmer_handle, delimiter="\t", lineterminator="\n")
        definition_writer = csv.writer(definition_handle, delimiter="\t", lineterminator="\n")
        kmer_writer.writerow(["assembly", "FN_site", "haplotype", "kmer_seq", "kmer_len"])
        definition_writer.writerow(
            [
                "assembly",
                "FN_site",
                "haplotype",
                "window",
                "truth_variants",
                "applied_ALT_variants",
                "excluded_opposite_haplotype_variants",
                "phase_evidence",
                "kmer_seq",
                "kmer_len",
                "kmer_reference_map",
            ]
        )

        for site in sites:
            start0 = site.pos - 1 - args.kmer_half
            end0 = site.pos - 1 + len(site.ref) + args.kmer_half
            if start0 < 0:
                raise ValueError(f"Marker window starts before {site.chrom}:1 for {site.label}")
            reference_sequence = fasta.fetch(site.chrom, start0, end0).upper()
            expected_reference_length = 2 * args.kmer_half + len(site.ref)
            if len(reference_sequence) != expected_reference_length:
                raise ValueError(f"Unexpected reference window length for {site.label}")

            variants = load_window_variants(truth, truth_sample, site, start0, end0)
            target_matches = [
                variant
                for variant in variants
                if variant.pos == site.pos and variant.ref == site.ref and variant.alt == site.alt
            ]
            if len(target_matches) != 1:
                raise ValueError(f"Expected one exact target truth record for {site.label}")
            target = target_matches[0]

            if target.phased:
                variant_sets = phased_variant_sets(variants, target)
            else:
                variant_sets = unphased_variant_sets(
                    hifi,
                    variants,
                    target,
                    args.min_base_quality,
                    args.min_phase_reads,
                    args.phase_ratio,
                )

            for haplotype, applied, excluded, phase_notes in variant_sets:
                kmer, reference_map = apply_variants(reference_sequence, start0, applied)
                if len(kmer) != len(reference_map):
                    raise ValueError(f"Sequence/map length mismatch for {site.label}")
                if not contains_alt_allele(kmer, reference_map, target):
                    raise ValueError(
                        f"Generated haplotype {haplotype} does not contain target ALT for "
                        f"{site.label}"
                    )

                kmer_writer.writerow(
                    [site.assembly, site.label, haplotype, kmer, len(kmer)]
                )
                definition_writer.writerow(
                    [
                        site.assembly,
                        site.label,
                        haplotype,
                        f"{site.chrom}:{start0 + 1}-{end0}",
                        ";".join(variant.label for variant in variants),
                        ";".join(variant.label for variant in applied),
                        ";".join(variant.label for variant in excluded) or ".",
                        ";".join(phase_notes) or "not_needed",
                        kmer,
                        len(kmer),
                        format_reference_map(reference_map),
                    ]
                )

    hifi.close()
    truth.close()
    fasta.close()


if __name__ == "__main__":
    main()
