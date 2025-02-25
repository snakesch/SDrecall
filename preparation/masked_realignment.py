import pysam
from pyfaidx import Fasta
import subprocess
import tempfile
import os

def modify_and_sort_genome_coordinates(mid_align_path, ref_contig_sizes_path, output_bam_path, pc_index="PC1"):
    """
    Modifies masked genome coordinates in a SAM file, sorts by coordinates,
    and outputs a sorted and indexed BAM file.

    Args:
        mid_align_path (str): Path to the input SAM file.
        ref_contig_sizes_path (str): Path to the FAI index file.
        output_bam_path (str): Path to the output BAM file.
        pc_index (str, optional): Prefix for read names. Defaults to "PC1".

    Returns:
        bool: True if successful, False otherwise.
    """

    if not ref_contig_sizes_path:
        print("No reference contig sizes provided.")
        return False

    try:
        contig_sizes = {}
        with open(ref_contig_sizes_path, 'r') as fai_file:
            for line in fai_file:
                fields = line.strip().split('\t')
                contig_sizes[fields[0]] = int(fields[1])
    except FileNotFoundError:
        print(f"Error: FAI index file not found at {ref_contig_sizes_path}")
        return False

    try:
        with pysam.AlignmentFile(mid_align_path, "r") as samfile:
            header = samfile.header.to_dict()

            existing_sq = {sq['SN']: sq for sq in header['SQ']}
            for contig, length in contig_sizes.items():
                sq_entry = {'SN': contig, 'LN': length}
                if "_" in contig:
                    sq_entry['AH'] = '*'
                if contig in existing_sq:
                    existing_sq[contig].update(sq_entry)
                else:
                    header['SQ'].append(sq_entry)
            
            #Add the 'SO' (sort order) tag
            header['HD']['SO'] = 'coordinate'


            with pysam.AlignmentFile(output_bam_path, "wb", header=pysam.AlignmentHeader.from_dict(header)) as outbam:
                for read in samfile:
                    qname = f"{read.query_name}:{pc_index}"

                    rname = read.reference_name
                    rnext = read.next_reference_name

                    if rname is not None:
                        pos = read.reference_start + 1
                    else:
                        pos = 0
                        rname = "*"

                    if rnext == "=":
                        rnext = rname
                        pnext = read.next_reference_start + 1 + (pos - 1)
                    elif rnext is not None:
                        pnext = read.next_reference_start + 1
                    else:
                        pnext = 0
                        rnext = "*"

                    a = pysam.AlignedSegment(outbam.header)
                    a.query_name = qname
                    a.query_sequence = read.query_sequence
                    a.flag = read.flag
                    a.reference_id = outbam.header.get_tid(rname)
                    a.reference_start = pos - 1
                    a.mapping_quality = read.mapping_quality
                    a.cigar = read.cigar
                    a.next_reference_id = outbam.header.get_tid(rnext)
                    a.next_reference_start = pnext - 1
                    a.template_length = read.template_length
                    a.query_qualities = read.query_qualities

                    for tag, value in read.tags:
                        a.set_tag(tag, value)

                    outbam.write(a)

        # Sort and index the BAM file
        pysam.sort("-o", output_bam_path, output_bam_path)  # In-place sort
        pysam.index(output_bam_path)
        return True


    except FileNotFoundError:
        print(f"Error: SAM file not found at {mid_align_path}")
        return False
    except pysam.utils.SamtoolsError as e:
        print(f"Error processing SAM/BAM file: {e}")
        return False

def create_minimap2_index(fasta_path, index_path, threads=1):
    """
    Creates a minimap2 index from a FASTA file.

    Args:
        fasta_path (str): Path to the input FASTA file.
        index_path (str): Path to the output minimap2 index file (.mmi).
        threads (int): Number of threads to use.

    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        command = [
            "minimap2",
            "-d", index_path,  # Output index file
            "-t", str(threads),
            fasta_path      # Input FASTA file
        ]
        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating minimap2 index: {e}")
        return False
    except FileNotFoundError:
        print("Error: minimap2 executable not found.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred creating the index: {e}")
        return False

def minimap2_align(forward_reads, reverse_reads, minimap2_index, mode, pc_label, output_bam_path, ref_contig_sizes_path, threads=1):
    """
    Aligns reads using minimap2, processes output.

    Args:
        forward_reads (str): Forward reads FASTQ.
        reverse_reads (str): Reverse reads FASTQ (or None).
        minimap2_index (str): Minimap2 index.
        mode (str): Minimap2 preset mode.
        pc_label (str): Label for read group and prefix.
        output_bam_path (str): Output sorted and indexed BAM.
        ref_contig_sizes_path (str): Path to FAI index
        threads (int): Number of threads. Default 1.
    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        # Create a temporary file for the SAM output
        with tempfile.NamedTemporaryFile(delete=False, suffix=".sam") as temp_sam_file:
            temp_align = temp_sam_file.name

            # Construct minimap2 command
            command = [
                "minimap2",
                "-ax", mode,
                "--eqx",
                "--MD",
                "-F", "1000",
                "-t", str(threads),
                "-R", f"@RG\\tID:{pc_label}\\tLB:SureSelectXT\\tPL:ILLUMINA\\tPU:1064\\tSM:{pc_label}",
                minimap2_index
            ]
            command.append(forward_reads)
            if reverse_reads:
                command.append(reverse_reads)

            # Run minimap2, capturing output to the temporary SAM file
            subprocess.run(command, stdout=temp_sam_file, check=True)

        # Modify, sort, and index
        if not modify_and_sort_genome_coordinates(temp_align, ref_contig_sizes_path, output_bam_path, pc_label):
            print("Error: Failed to modify and sort coordinates.")
            os.remove(temp_align) # Clean up even on failure
            return False

        # Clean up the temporary file
        os.remove(temp_align)
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e}")
        return False
    except FileNotFoundError:
        print("Error: minimap2 executable not found.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False