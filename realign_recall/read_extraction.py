import os

from src.utils import executeCmd
from src.log import logger

# Try to import Rust module
RUST_AVAILABLE = False
rust_bam_to_fastq = None

try:
    # Try to import the compiled extension directly
    import read_extraction
    if hasattr(read_extraction, 'bam_to_fastq_biobambam'):
        rust_bam_to_fastq = read_extraction.bam_to_fastq_biobambam
        RUST_AVAILABLE = True
        logger.info("Rust module loaded successfully")
except ImportError:
    pass
    
if not RUST_AVAILABLE:
    logger.warning("Rust module not available, falling back to shell commands")


def bam_to_fastq_biobambam(input_bam, 
                           region_bed, 
                           output_freads, 
                           output_rreads, 
                           multi_aligned=False,
                           threads=1, 
                           tmp_dir="/tmp",
                           logger=logger,
                           use_rust=True):
    """
    Extract read pairs where at least one read overlaps the region using biobambam2.
    
    Parameters:
        input_bam: Input BAM file
        region_bed: BED file with regions of interest
        output_freads: Output R1 FASTQ file
        output_rreads: Output R2 FASTQ file
        multi_aligned: Whether the input BAM is multi-aligned
        threads: Number of threads to use
        tmp_dir: Temporary directory for intermediate files
        logger: Logger object
        use_rust: Whether to use Rust implementation if available
    
    Returns:
        Tuple containing paths to output R1 and R2 FASTQ files
    """
    # Use Rust implementation if available and requested
    if use_rust and RUST_AVAILABLE:
        try:
            logger.info("Using Rust implementation for BAM to FASTQ conversion")
            return rust_bam_to_fastq( input_bam, 
                                      region_bed, 
                                      output_freads, 
                                      output_rreads,
                                      multi_aligned, 
                                      threads, 
                                      tmp_dir
                                    )
        except Exception as e:
            logger.error(f"Rust implementation failed: {e}")
            logger.info("Falling back to shell command implementation")
    
    # Original shell command implementation
    logger.info("Using shell command implementation for BAM to FASTQ conversion")
    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_freads), exist_ok=True)
    os.makedirs(os.path.dirname(output_rreads), exist_ok=True)
            
    # Using biobambam2 bamtofastq with ranges
    tmp_prefix = f"tmp_bb_{os.getpid()}"

    if multi_aligned:
        # Filter expression for multi-aligned reads (converted to samtools syntax)
        # Note: Rust implementation uses ![SA] && [XA] && abs(AS-XS) <= 5
        filter_expr = "![SA] && ([XA] || mapq < 50)"
        
        # Samtools command with -P for fetching pairs
        cmd = f"""samtools view -@ {threads} -h -P -L {region_bed} -u -e '{filter_expr}' {input_bam} | \
                  bamtofastq \
                    F={output_freads} \
                    F2={output_rreads} \
                    collate=1 \
                    T={tmp_dir}/{tmp_prefix} && \
                  rm -rf {tmp_dir}/{tmp_prefix}"""
    else:
        # For non-multi-aligned BAMs, use the original approach with regions
        cmd = f"""samtools view -h -P -@ {threads} -L {region_bed} -u {input_bam} | \
                  bamtofastq \
                    F={output_freads} \
                    F2={output_rreads} \
                    collate=1 \
                    T={tmp_dir}/{tmp_prefix} && \
                  rm -rf {tmp_dir}/{tmp_prefix}"""
    
    executeCmd(cmd, logger=logger)
    
    # Verify output files were created with content
    cmd = f"wc -l {output_freads} | awk '{{printf \"%d\", $1}}'"
    num_freads = executeCmd(cmd, logger=logger)
    cmd = f"wc -l {output_rreads} | awk '{{printf \"%d\", $1}}'"
    num_rreads = executeCmd(cmd, logger=logger)
    logger.info(f"The number of reads in {output_freads} is {num_freads}, and the number of reads in {output_rreads} is {num_rreads}\n")
    if os.path.exists(output_freads) and \
       os.path.exists(output_rreads) and \
       os.path.getsize(output_freads) > 0 and \
       os.path.getsize(output_rreads) > 0 and \
       num_freads == num_rreads:
        return output_freads, output_rreads
    else:
        logger.error(f"Failed to extract reads for region {region_bed} from BAM {input_bam}")
        return None, None