import os

from src.utils import executeCmd
from src.log import logger

def bam_to_fastq_biobambam(input_bam, 
                           region_bed, 
                           output_freads, 
                           output_rreads, 
                           multi_aligned=False,
                           threads=1, 
                           tmp_dir="/tmp",
                           logger=logger):
    """
    Extract read pairs where at least one read overlaps the region using biobambam2.
    
    Parameters:
        input_bam: Input BAM file
        region_bed: BED file with regions of interest
        output_freads: Output R1 FASTQ file
        output_rreads: Output R2 FASTQ file
        multi_aligned: Whether the input BAM is multi-aligned
        threads: Number of threads to use (biobambam2 doesn't support threading in the same way)
        tmp_dir: Temporary directory for intermediate files
        logger: Logger object
    
    Returns:
        Tuple containing paths to output R1 and R2 FASTQ files
    """
    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_freads), exist_ok=True)
    os.makedirs(os.path.dirname(output_rreads), exist_ok=True)
            
    # Using biobambam2 bamtofastq with ranges
    tmp_prefix = f"tmp_bb_{os.getpid()}"

    if multi_aligned:
        filter_expr = "[SA] == null and ([XA] != null or ([XS] != null and [AS] != null and ([AS] - [XS] <= 5)) or mapping_quality <= 40)"
        cmd = f"""sambamba view -t {threads} -h -f unpack -L {region_bed} -F "{filter_expr}" {input_bam} | \
                  bamtofastq \
                    filename=stdin \
                    F={output_freads} \
                    F2={output_rreads} \
                    collate=1 \
                    exclude=UNPAIRED \
                    T={tmp_dir}/{tmp_prefix} && \
                  rm -rf {tmp_dir}/{tmp_prefix}"""
    else:
        # Convert BED file to ranges format expected by bamtofastq
        # This is necessary as biobambam2 doesn't accept BED files directly
        ranges_str = ""
        with open(region_bed, 'r') as bed:
            for line in bed:
                if line.strip() and not line.startswith('#'):
                    cols = line.strip().split('\t')
                    if len(cols) >= 3:
                        chrom, start, end = cols[0], int(cols[1])+1, int(cols[2])+1  # Convert 0-based BED to 1-based
                        ranges_str += f"{chrom}:{start}-{end} "
        cmd = f"""bamtofastq \
                    filename={input_bam} \
                    ranges="{ranges_str}" \
                    F={output_freads} \
                    F2={output_rreads} \
                    collate=1 \
                    T={tmp_prefix}"""
    
    executeCmd(cmd, logger=logger)
    
    # Verify output files were created with content
    if os.path.exists(output_freads) and os.path.exists(output_rreads) and \
       os.path.getsize(output_freads) > 0 and os.path.getsize(output_rreads) > 0:
        return output_freads, output_rreads
    else:
        logger.error(f"Failed to extract reads for region {region_bed} from BAM {input_bam}")
        return None, None