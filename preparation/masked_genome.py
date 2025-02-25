import os
import logging
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
from pybedtools import BedTool

from .inferred_depths import create_genome_dict

logger = logging.getLogger("SDrecall")

def create_masked_genome(input_genome, bed_file, output_genome, avg_insert_size=400, std_insert_size=130) -> str:
    """Masks regions in a genome FASTA based on a BED file."""

    if not (os.path.exists(input_genome) and os.path.exists(bed_file)):
        raise FileNotFoundError("Input genome or BED file not found.")

    try:
        fai = f"{input_genome}.fai"
        if not os.path.exists(fai) or os.path.getmtime(input_genome) > os.path.getmtime(fai):
            result = subprocess.run(["samtools", "faidx", input_genome], capture_output=True, text=True, check=True)
            if not os.path.exists(fai):
                raise ValueError(f"Failed to create .fai index. samtools output: {result.stderr}")

        pad_size = int(avg_insert_size + 2 * std_insert_size)
        slop_b = int(avg_insert_size + 2 * std_insert_size + 1000)

        ref = Fasta(input_genome)
        genome_dict = create_genome_dict(fai)

        target_bed = BedTool(bed_file).sort().merge(d=2 * pad_size).slop(b=slop_b, genome=genome_dict).sort()
        logger.debug(f"Target bed {bed_file}, after processing: {target_bed.total_coverage()} bp")

        mask_seq = "N" * 1000
        masked_genome_contigs = []

        for interval in target_bed:
          try:
            original_seq = ref[interval.chrom][interval.start:interval.end].seq
            masked_seq = mask_seq + str(original_seq)[1000:-1000] + mask_seq

            if str(masked_seq) != str(original_seq):
              if masked_seq.count('N') < 2000:
                 logger.error(f"Failed to properly mask region {interval.chrom}:{interval.start}-{interval.end}")
                 return None

            masked_record = SeqRecord(Seq(str(ref[interval.chrom][:].seq)), id=interval.chrom, description="")
            masked_seq_list = list(str(masked_record.seq))
            masked_seq_list[interval.start:interval.end] = list(masked_seq)
            masked_record.seq = Seq("".join(masked_seq_list))

            masked_genome_contigs.append(masked_record)
          except KeyError as e:
            logger.error(f"Chromosome {e} not found in the FASTA index.")
            return None
          except Exception as e:
            logger.exception(f"Error during masking of interval {interval}: {e}")
            return None

        with open(output_genome, "w") as fout:
          SeqIO.write(masked_genome_contigs, fout, "fasta")

        if not os.path.exists(output_genome):
            logger.error(f"Failed to write {output_genome}")
            return None

        result = subprocess.run(["samtools", "faidx", output_genome], capture_output=True, text=True, check=True)
        if not os.path.exists(f"{output_genome}.fai"):
            logger.error(f"Failed to index {output_genome}. samtools output: {result.stderr}")
            return None

        return output_genome

    except Exception as e:
        logger.exception(f"An error occurred: {e}")
        return None
    
def get_reference_seq(bedf_path: str, fastq_path: str, ref_genome: str, padding = 0):
    """
    Extracts sequences from a reference genome based on a BED file,
    optionally padding the regions, and writes the output to a FASTQ file.
    Uses pybedtools and pyfaidx instead of subprocess calls.
    """
    from tempfile import NamedTemporaryFile

    ref_fasta = Fasta(ref_genome)  # Load the reference genome with pyfaidx
    bedtool = BedTool(bedf_path)

    if padding > 0:
        # Create a temporary file for the padded BED
        with NamedTemporaryFile(dir="/tmp", suffix=".bed", delete=False) as tmp_bedf:
            tmp_bedf_path = tmp_bedf.name

        genome_dict = create_genome_dict(ref_genome + ".fai")

        # Use bedtools slop to add padding.  .fai file is not necessary with BedTool object
        bedtool = bedtool.slop(b=padding, genome=genome_dict)  # Use genome string directly
        bedtool.saveas(tmp_bedf_path)
        bedtool = BedTool(tmp_bedf_path)  # Reload from the temporary file

    with open(fastq_path, "w") as outfile:
        for interval in bedtool:
            # Extract sequence using pyfaidx.  Handles chromosome names correctly.
            seq_name = interval.chrom
            start = interval.start
            end = interval.end
            strand = interval.strand if interval.strand in ["+", "-"] else "+"

            sequence = ref_fasta[seq_name][start:end]

            # Correctly handling strand
            if strand == "-":
                sequence = sequence.complement[::-1]
                
            sequence_str = str(sequence).upper() # Convert to uppercase string

            # Write to FASTQ format (simulated, since we don't have quality scores)
            outfile.write(f"@{interval.name}\n") # Use interval name. Create name if it does not exist
            outfile.write(f"{sequence_str}\n")
            outfile.write("+\n")  # Placeholder for quality scores
            outfile.write("F" * len(sequence_str) + "\n")  # Placeholder quality

    if padding > 0:
        # Remove temporary padded BED file
        import os
        os.remove(tmp_bedf_path)
        
    return fastq_path