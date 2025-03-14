import os
import uuid
import pandas as pd
import logging
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta

from pybedtools import BedTool

from src.utils import executeCmd, update_plain_file_on_md5
from src.log import logger

class Genome:
    """
    A class implementation for genome masking.

    Usage:
    rg = Genome(path)
    rg.mask(bedf)  # bedf needs not to be complemented
    """

    def __init__(self, path):
        self.path = path
        self.fai_index = self._create_index("fai", "samtools faidx")
        self.contigsize = self.get_contig_genome()

    def _create_index(self, ext, cmd):
        index_path = f"{self.path}.{ext}"
        if not os.path.exists(index_path) or os.path.getmtime(self.path) > os.path.getmtime(index_path):
            executeCmd(cmd + f" {self.path}")
        return index_path

    def get_contig_genome(self):
        contig_genome_fn = ".".join(self.path.split(".")[:-1]) + ".contigsize.genome"
        pd.read_csv(self.fai_index, sep="\t", header=None).iloc[:, [0, 1]].to_csv(contig_genome_fn, sep="\t", index=False, header=False)
        return contig_genome_fn

    def mask(self, bedf: str, avg_frag_size=400, std_frag_size=130, genome="hg19", logger=None, path=""):
        if not os.path.exists(bedf):
            raise FileNotFoundError(f"Invalid BED file: {bedf}")
        if path == "":
            raise ValueError("No output path for masked genome. ")

        region = self._extract_region_name(bedf)
        target_bed = self._prepare_target_bed(bedf, avg_frag_size, std_frag_size, genome, logger)

        ref_genome_seq = Fasta(self.path)
        masked_genome_contigs = self._mask_intervals(target_bed, ref_genome_seq, logger)

        self._write_masked_genome(masked_genome_contigs, region, path)

        return path

    def _extract_region_name(self, bedf):
        return os.path.basename(bedf).split("_")[0] if "_related_homo_regions" in bedf else os.path.basename(bedf).split(".")[0]

    def _prepare_target_bed(self, bedf, avg_frag_size, std_frag_size, genome, logger):
        target_bed = (BedTool(bedf).sort().merge()
                      .slop(b=int(avg_frag_size + 2 * std_frag_size + 1000), g=genome)
                      .sort().merge())
        logger.debug(f"Target bed {bedf}, after padding: {target_bed.total_coverage()} bp covered.")
        return target_bed

    def _mask_intervals(self, target_bed, ref_genome_seq, logger):
        mask_seq = "N" * 1000
        return [SeqRecord(self._apply_mask(ref_genome_seq[interval.chrom][interval.start:interval.stop].seq, logger, interval, mask_seq),
                          id=f"{interval.chrom}:{interval.start}", description="") for interval in target_bed]

    def _apply_mask(self, interval_seq, logger, interval, mask_seq):
        masked_seq = mask_seq + str(interval_seq)[1000:-1000] + mask_seq
        if str(masked_seq) != str(interval_seq):
            logger.debug(f"Sequence already contains Ns at both ends: {interval.chrom}:{interval.start}-{interval.stop}")
            if masked_seq.count("N") < 2000:
                logger.error("Sequence already contains N at both ends, SDrecall was unable to add N. ")
                sys.exit(1)
        return Seq(masked_seq)

    def _write_masked_genome(self, masked_genome_contigs, region, path):
        tmp_tag = str(uuid.uuid4())
        tmp_genome = path.replace(".fasta", f".{tmp_tag}.fasta")

        SeqIO.write(masked_genome_contigs, tmp_genome, "fasta")
        self._update_masked_genome(path, tmp_genome)

        return 

    def _update_masked_genome(self, masked_genome, tmp_genome):
        # Update the masked genome based on MD5
        updated = update_plain_file_on_md5(masked_genome, tmp_genome)
        Genome(masked_genome)  # Re-initialize to build indices
        return updated
