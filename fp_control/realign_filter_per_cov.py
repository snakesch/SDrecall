import os
import gc
import pysam
import logging
import traceback

from src.utils import executeCmd, prepare_tmp_file
from realign_recall.annotate_HP_tag_to_vars import annotate_vcf as annotate_vcf_HP_tag
from src.const import shell_utils
from src.log import logger

from fp_control.bam_ncls import migrate_bam_to_ncls, calculate_mean_read_length


def imap_filter_out(args):
    """Worker function that processes a single region and returns results with log file path"""
    # Unpack arguments
    args, log_dir, log_level = args
    raw_bam, output_bam, intrinsic_bam, bam_region_bed, recall_mq_cutoff, basequal_median_cutoff, edge_weight_cutoff, numba_threads, tmp_dir, ref_genome, sample_id, job_id = args
    import numba
    numba.set_num_threads(numba_threads)
    numba.config.THREADING_LAYER = 'tbb'

    # Create unique log file for this subprocess
    log_file = os.path.join(log_dir, f"subprocess_{job_id}_{os.path.basename(raw_bam)}.log")
    
    # Configure file logger
    file_handler = logging.FileHandler(log_file)
    # Try to get formatter from existing handlers
    parent_formatter = None
    if logger.handlers:
        for handler in logger.handlers:
            if handler.formatter:
                parent_formatter = handler.formatter
                break
            
    if parent_formatter:
        file_handler.setFormatter(parent_formatter)
    else:
        file_handler.setFormatter(logging.Formatter("[%(asctime)s] [%(levelname)s] [%(pathname)s:%(funcName)s:%(lineno)d] %(message)s"))
    
    # Create logger that writes to file
    subprocess_logger = logging.getLogger(f"SubProcess-{job_id}")
    subprocess_logger.setLevel(log_level)  # Use logger level from src/log.py logger
    subprocess_logger.addHandler(file_handler)
    subprocess_logger.propagate = False  # Don't send logs to parent loggers
    
    try:
        # Log the start of processing with all parameters for debugging
        subprocess_logger.warning(f"Log level: {subprocess_logger.level}")
        subprocess_logger.warning(f"Thread settings: numba={numba.get_num_threads()}, " + 
                                  f"OMP={os.environ.get('OMP_NUM_THREADS')}, " + 
                                  f"MKL={os.environ.get('MKL_NUM_THREADS')}")
        subprocess_logger.info(f"Processing region in {raw_bam}")
        subprocess_logger.info(f"Parameters: output_bam={output_bam}, intrinsic_bam={intrinsic_bam}, "
                               f"bam_region_bed={bam_region_bed}, "
                               f"recall_mq_cutoff={recall_mq_cutoff}, basequal_median_cutoff={basequal_median_cutoff}, "
                               f"edge_weight_cutoff={edge_weight_cutoff}")
        
        # Call the main processing function
        phased_graph, output_bam, output_vcf = realign_filter_per_cov(
            bam=raw_bam,
            output_bam=output_bam,
            intrinsic_bam=intrinsic_bam,
            bam_region_bed=bam_region_bed,
            recall_mq_cutoff=recall_mq_cutoff,
            basequal_median_cutoff=basequal_median_cutoff,
            edge_weight_cutoff=edge_weight_cutoff,
            tmp_dir=tmp_dir,
            ref_genome=ref_genome,
            sample_id=sample_id,
            logger=subprocess_logger
        )
        gc.collect()
        
        # If successful, return result with comma-separated fields and log file path
        if phased_graph is None:
            subprocess_logger.warning(f"Cannot filter {raw_bam} because the phased graph cannot be built")
            return f"{raw_bam},NaN,NaN,{log_file}"
        
        if output_vcf is None and output_bam is not None:
            subprocess_logger.warning(f"Failed to call variants on the filtered BAM file {output_bam}")
            return f"{raw_bam},{output_bam},NaN,{log_file}"
            
        subprocess_logger.info(f"Successfully processed {raw_bam}")
        return f"{raw_bam},{output_bam},{output_vcf},{log_file}"
        
    except Exception as e:
        # Log the full stack trace for debugging
        subprocess_logger.error(f"Error processing {raw_bam}: {str(e)}")
        subprocess_logger.error(traceback.format_exc())
        gc.collect()
        # Return error result with log file path
        return f"{raw_bam},NaN,NaN,{log_file}"



def weight_matrix_to_dataframe(weight_matrix, graph):
    """
    Convert a weight matrix indexed by graph vertices to a pandas DataFrame
    with query names as row and column labels, using the graph's vertex property map.
    
    Parameters:
    -----------
    weight_matrix : numpy.ndarray
        The square matrix of edge weights, indexed by vertex indices.
    
    graph : graph_tool.Graph
        The graph object containing vertex properties with query names.
        
    Returns:
    --------
    pandas.DataFrame
        A DataFrame where both rows and columns are labeled with query names,
        and values represent edge weights between pairs of reads.
    """
    import pandas as pd
    import numpy as np
    
    # Get direct access to the qname property map
    qname_prop = graph.vertex_properties["qname"]
    
    # Get all vertices that are valid indices in the weight matrix
    vertex_indices = []
    vertex_qnames = []
    
    for v in graph.vertices():
        v_idx = int(v)
        if v_idx < weight_matrix.shape[0]:  # Ensure vertex index is within matrix bounds
            vertex_indices.append(v_idx)
            vertex_qnames.append(qname_prop[v])
    
    # Create the DataFrame using qnames from the property map
    # First create a view of the weight matrix with only the relevant indices
    submatrix = weight_matrix[np.ix_(vertex_indices, vertex_indices)]
    
    # Create the DataFrame with query names as labels
    df = pd.DataFrame(
        data=submatrix,
        index=vertex_qnames,
        columns=vertex_qnames
    )
    
    return df


def realign_filter_per_cov(bam,
                           output_bam = None,
                           intrinsic_bam = None,
                           bam_region_bed = None,
                           recall_mq_cutoff = 10,
                           basequal_median_cutoff = 15,
                           edge_weight_cutoff = 0.201,
                           threads = 4,
                           tmp_dir = "/tmp",
                           ref_genome = None,
                           sample_id = None,
                           logger=logger):
    '''
    Main function for phasing and filtering realigned reads per continuous coverage region.

    This function performs the following main tasks:
    1. Migrates BAM files to NCLS format for efficient query of reads or read pairs overlapping a query genomic region
    2. Builds a read-pair graph from the BAM data for phasing, where each vertex represents a read pair and each edge represents the confidence that two read pairs are originated from the same haplotype
    3. Identifies haplotypes in the graph by finding non-overlapping maximal cliques iteratively (using an approximate but not exact algorithm for efficient clique search, basically we iteratively run Greedy-Clique-Expansion algorithm)
    4. After read pair grouped into different haplotypes, we put them into a binary integer linear programming model to solve the haplotype-level misalignment with HiGHs solver.
    5. Generates output BAM files with filtered and annotate haplotype index assigned to each read pairs.

    Parameters:
    - bam (str): Path to the input BAM file
    - output_bam (str, optional): Path for the output filtered BAM file, if not specified, output_bam will be the input bam with .clean.bam suffix
    - intrinsic_bam (str): Path to the intrinsic BAM file, generated earlier in the SDrecall workflow (basically it aligns the reference sequence of one SD to the other SD)
    - raw_vcf (str): Path to the raw VCF file, the variants detected in the raw pooled alignments.
    - bam_region_bed (str, optional): Path to the BED file defining covered regions of the processing bam file, if not specified, bam_region_bed will be generated with a name after the input bam with .coverage.bed suffix
    - recall_mq_cutoff (int): Mapping quality cutoff to be included in the analysis
    - basequal_median_cutoff (int): Base quality median cutoff to be included in the analysis (if the median BQ of a read is lower than this cutoff, we will discard this read because it is too noisy)
    - edge_weight_cutoff (float): Edge weight cutoff separating two rounds of BK clique searches
    - logger (Logger object): Logger for output messages

    Returns:
    - phased_graph (Graph object): The constructed phasing graph

    This function integrates various analysis steps including BAM processing,
    graph construction, haplotype identification, and read filtering to improve
    the quality of genomic alignments and identify potential misalignments.
    '''

    # Import the auto-selection function that chooses between Rust and Python implementations
    from fp_control.graph_build import build_phasing_graph_auto as build_phasing_graph, RUST_AVAILABLE
    from fp_control.identify_misaligned_haps import inspect_by_haplotypes
    from fp_control.phasing import phasing_realigned_reads

    chunk_id = os.path.basename(bam).split(".")[-2]

    # Use original Python implementation with NCLS preprocessing
    bam_ncls = migrate_bam_to_ncls(bam,
                                   mapq_filter = recall_mq_cutoff,
                                   basequal_median_filter = basequal_median_cutoff,
                                   logger=logger)

    logger.info(f"Successfully migrated the BAM file {bam} to NCLS format, this part is necessary for both Python and Rust implementation\n\n")

    # Now migrate the intrinsic BAM file to NCLS format
    intrin_bam_ncls = migrate_bam_to_ncls(intrinsic_bam,
                                          mapq_filter = 0,
                                          basequal_median_filter = 0,
                                          paired = False,
                                          filter_noisy = False,
                                          logger=logger)
    # Since intrinsic BAM reads are reference sequences, therefore there are no low quality reads
    intrin_bam_ncls = intrin_bam_ncls[:-1]

    # Given the top 1% mismatch count per read (one structrual variant count as 1 mismatch)
    tmp_bam = prepare_tmp_file(suffix=".bam", tmp_dir = tmp_dir).name

    noisy_qnames = set([])
    mismap_qnames = set([])
    norm_qnames = {"":set([])}
    noisy_num = 0
    total_num = 0

    # Check if Rust implementation is available to decide the processing path
    if RUST_AVAILABLE:
        logger.info("Using Rust-accelerated processing - bypassing Python BAM preprocessing")
        # Calculate mean read length directly for Rust implementation
        mean_read_length = calculate_mean_read_length(bam)
        logger.info(f"Average read length: {mean_read_length}")
        
        # Call Rust implementation directly - it handles BAM processing internally
        phased_graph, weight_matrix, qname_to_node, total_readhap_vector, total_readerr_vector, read_ref_pos_dict, total_lowqual_qnames = build_phasing_graph(bam,
                                                                                                                                                              intrinsic_bam,
                                                                                                                                                              None,  # ncls_dict - not used by Rust
                                                                                                                                                              None,  # read_dict - not used by Rust
                                                                                                                                                              None,  # qname_dict - not used by Rust
                                                                                                                                                              mean_read_length,
                                                                                                                                                              set(),  # total_lowqual_qnames - empty initially
                                                                                                                                                              ref_genome,
                                                                                                                                                              logger = logger)
        if phased_graph is None:
            return None, None, None
            
        logger.info(f"Rust implementation: Graph built with {phased_graph.num_vertices()} vertices and {phased_graph.num_edges()} edges")
        
    else:
        logger.info("Using Python implementation - performing BAM preprocessing")

        # parse the results from the tuple returned by migrate_bam_to_ncls
        ncls_dict, read_dict, qname_dict, qname_idx_dict, total_lowqual_qnames = bam_ncls
        
        logger.info(f"Successfully migrated the intrinsic BAM file {intrinsic_bam} to NCLS format\n")
        logger.info(f"Intrinsic BAM contains {len(intrin_bam_ncls[1])} reads in total.\n\n")
        
        # Calculate the mean read length of the input bam file
        mean_read_length = calculate_mean_read_length(bam)
        logger.info(f"Average read length: {mean_read_length}")

        # Create the read-pair graph using Python implementation
        phased_graph, weight_matrix, qname_to_node, total_readhap_vector, total_readerr_vector, read_ref_pos_dict, total_lowqual_qnames = build_phasing_graph(bam,
                                                                                                                                                              intrinsic_bam,
                                                                                                                                                              ncls_dict,
                                                                                                                                                              read_dict,
                                                                                                                                                              qname_dict,
                                                                                                                                                              mean_read_length,
                                                                                                                                                              total_lowqual_qnames,
                                                                                                                                                              ref_genome,
                                                                                                                                                              logger = logger)
        if phased_graph is None:
            return None, None, None

    # Common code for both implementations
    bam_graph = bam.replace(".bam", ".phased.graphml")
    bam_weight_matrix = bam.replace(".bam", ".weight_matrix.tsv")

    logger.info(f"Now succesfully built the phasing graph with {phased_graph.num_vertices()} vertices and {phased_graph.num_edges()} edges. Save it to {bam_graph}, the weight matrix is saved to {bam_weight_matrix}\n\n")
    # Now we need to extract the components in the phased graph
    phased_graph.save(bam_graph)
    # weight_matrix_to_dataframe(weight_matrix, phased_graph).to_csv(bam_weight_matrix, sep = "\t", index = True)

    # Now we need to do local phasing for each component in the graph. (Finding non-overlapping high edge weight cliques inside each component iteratively)
    qname_hap_info, hap_qname_info = phasing_realigned_reads(phased_graph,
                                                             weight_matrix,
                                                             edge_weight_cutoff,
                                                             logger = logger)

    # Inspect the raw BAM corresponding variants to get the high density regions
    # It's like active region identification for GATK HC
    if not bam_region_bed:
        bam_region_bed = bam.replace(".bam", ".coverage.bed")
        cmd = f"bash {shell_utils} samtools_bam_coverage \
                -i {bam} \
                -d 0 \
                -o {bam_region_bed}"
        executeCmd(cmd, logger = logger)

    total_genomic_haps = {}
    compare_haplotype_meta_tab = bam.replace(".bam", ".haplotype_meta.tsv")

    if len(hap_qname_info) <= 2:
        logger.warning(f"Only {len(hap_qname_info)} haplotype clusters are found for bam {bam}. Do not need to choose 2 haplotypes, Skip this region.\n")
        correct_qnames, mismap_qnames = set([qn for hs in list(hap_qname_info.values()) for qn in hs]), set([])
    else:
        correct_qnames, mismap_qnames = inspect_by_haplotypes(bam,
                                                              hap_qname_info,
                                                              qname_hap_info,
                                                              bam_ncls,
                                                              intrin_bam_ncls,
                                                              qname_to_node,
                                                              total_lowqual_qnames,
                                                              total_readhap_vector,
                                                              total_readerr_vector,
                                                              total_genomic_haps,
                                                              read_ref_pos_dict,
                                                              compare_haplotype_meta_tab = compare_haplotype_meta_tab,
                                                              mean_read_length = mean_read_length,
                                                              logger = logger )

    assert len(correct_qnames & mismap_qnames) == 0, f"The correct_qnames and mismap_qnames have overlap: {correct_qnames & mismap_qnames}"
    if len(correct_qnames) == 0:
        logger.warning(f"No read pairs are found to be correctly aligned in the target regions. In total only {len(hap_qname_info)} cliques (haplotypes) are found in the target inspected regions.")
        return None, None, None
    
    logger.info(f"In total has found {len(hap_qname_info)} clique separated components (haplotypes) in the target inspected regions.")
    logger.info(f"We found {len(mismap_qnames)} read pairs that are likely to be misaligned in the target regions.\n {mismap_qnames}\n And {len(correct_qnames)} read pairs that are likely to be correctly aligned in the target regions.\n {correct_qnames}\n")

    with pysam.AlignmentFile(bam, "rb") as bam_handle:
        # Extract the sample name (SM) from the existing read groups in the header
        with pysam.AlignmentFile(tmp_bam, "wb", header = bam_handle.header) as tmp_handle:
            with pysam.AlignmentFile(output_bam, "wb", header = bam_handle.header) as output_handle:
                # viewed_regions = {}
                for read in bam_handle:
                    if read.is_supplementary or \
                        read.is_secondary:
                        continue
                    total_num += 1
                    qname = read.query_name
                    hap_id = qname_hap_info.get(qname_to_node.get(qname, -1), "NA")

                    if qname in mismap_qnames:
                        hap_id = f"{hap_id}_HIGHVD"

                    if qname in total_lowqual_qnames:
                        hap_id = f"LOWQUAL"

                    # Filter out reads with oddly high editing distance that breakthrough the cutoff
                    if read.is_mapped and read.mapping_quality > 10:
                        if qname in total_lowqual_qnames:
                            logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is considered low base quality read, skip this pair of reads")
                            pass
                        elif qname in noisy_qnames:
                            # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is too high to be likely correctly mapped")
                            pass
                        elif qname in mismap_qnames:
                            # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is significantly higher than other reads so it is considered as noisy and misaligned")
                            pass
                        elif qname not in correct_qnames:
                            noisy_num += 1
                            logger.info(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is not in correct_qnames nor mismap_qnames, so it is considered as noisy and misaligned")
                            noisy_qnames.add(qname)
                        elif not read.is_secondary and \
                            not read.is_supplementary and \
                            not read.is_duplicate and \
                            not read.is_qcfail:
                            logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is considered appropriately aligned")
                            read.set_tag('HP', f'chunk{chunk_id}_{hap_id}')
                            qname_pair = norm_qnames.get(qname, set([]))
                            qname_pair.add(read)
                            norm_qnames[qname] = qname_pair
                    else:
                        logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is not mapped or has low mapping quality {read.mapping_quality}, skip this read")

                    if qname in noisy_qnames:
                        hap_id = f"{hap_id}_HIGHVD"
                        
                    read.set_tag('HP', f'chunk{chunk_id}_{hap_id}')
                    tmp_handle.write(read)

                # logger.warning(f"Check {check_odd_num} reads to find oddly high editing distance reads, {zero_odd_num} reads found no odd editing distance read pairs")
                for qname, pair in norm_qnames.items():
                    if qname in noisy_qnames:
                        logger.debug(f"Read pair {qname} is in noisy_qnames, skip this pair of reads")
                        continue
                    if qname in mismap_qnames:
                        logger.debug(f"Read pair {qname} is in mismap_qnames, skip this pair of reads")
                        continue
                    logger.debug(f"Write read pair {qname} to {output_bam}")
                    for read in pair:
                        output_handle.write(read)

    logger.warning(f"Filtered out {len(noisy_qnames)} noisy read-pairs and {len(mismap_qnames - noisy_qnames)} read-pairs with ODD high editing distance, remaining {len(set(norm_qnames.keys()) - noisy_qnames - mismap_qnames)} read-pairs from {bam} (with total {total_num} reads) and output to {output_bam}\n\n")

    # Replace the input BAM file with the tmp BAM file with modified RG tags for visualization of haplotype clusters
    executeCmd(f"samtools sort -O bam -@ {threads} -T {tmp_dir} -o {bam} {tmp_bam} && \
                 samtools index {bam} && \
                 rm {tmp_bam} && \
                 [[ $(samtools view {bam} | wc -l) -ge 1 ]] && \
                 ls -lh {bam}", logger=logger)    

    # Sort and index the output bam file
    try:
        executeCmd(f"samtools sort -O bam -@ {threads} -T {tmp_dir} -o {tmp_bam} {output_bam} && \
                     mv {tmp_bam} {output_bam} && \
                     samtools index {output_bam} && \
                     [[ $(samtools view {output_bam} | wc -l) -ge 1 ]] && \
                     ls -lh {output_bam}", logger=logger)
    except Exception as e:
        logger.warning(f"Failed to sort and index the output bam file {output_bam}, because its empty")
        return None, None, None
    
    # Call variants on the filtered BAM file
    output_vcf = output_bam.replace(".bam", ".vcf.gz")
    logger.info(f"Calling variants on the filtered BAM file {output_bam}")
    cmd = f"bash {shell_utils} bcftools_call_per_RG \
            -m {ref_genome} \
            -b {output_bam} \
            -o {output_vcf} \
            -c {threads} \
            -p {sample_id}"
    try:
        executeCmd(cmd, logger=logger)
    except Exception as e:
        logger.warning(f"Failed to call variants on the filtered BAM file {output_bam}, because {e}")
        return phased_graph, output_bam, None
    
    # Annotate variants with HP tags
    output_vcf = annotate_vcf_HP_tag(output_vcf, 
                                     output_vcf, 
                                     output_bam, 
                                     "HP", 
                                     logger = logger)
    
    return phased_graph, output_bam, output_vcf



if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Realign and filter reads per coverage region")
    
    # Required parameters
    parser.add_argument("--bam", type=str, required=True, help="Path to input BAM file")
    
    # Optional parameters
    parser.add_argument("--output_bam", type=str, help="Path for output filtered BAM file (default: input.clean.bam)")
    parser.add_argument("--intrinsic_bam", type=str, help="Path to intrinsic BAM file")
    parser.add_argument("--bam_region_bed", type=str, help="Path to BED file defining covered regions")
    parser.add_argument("--recall_mq_cutoff", type=int, default=10, help="Mapping quality cutoff")
    parser.add_argument("--basequal_median_cutoff", type=int, default=15, help="Base quality median cutoff")
    parser.add_argument("--edge_weight_cutoff", type=float, default=0.201, help="Edge weight cutoff")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--tmp_dir", type=str, default="/tmp", help="Temporary directory")
    parser.add_argument("--ref_genome", type=str, help="Path to reference genome")
    parser.add_argument("--sample_id", type=str, help="Sample identifier")
    parser.add_argument("--log_level", type=str, default="INFO", 
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level")
    
    args = parser.parse_args()
 
    # Configure threading safely:
    # - Avoid changing NUMBA_NUM_THREADS via environment after Numba may have initialized.
    # - Prefer runtime API; if already initialized, respect existing value.
    try:
        import numba
        try:
            numba.set_num_threads(args.threads)
        except RuntimeError as e:
            # Threads already launched; keep current setting to avoid errors
            try:
                cur = numba.get_num_threads()
                logger.warning(f"Numba threads already initialized at {cur}; skipping change to {args.threads} ({e})")
            except Exception:
                logger.warning(f"Numba threads already initialized; skipping change to {args.threads} ({e})")
    except Exception:
        # numba not available or import failed; continue
        pass

    # Set BLAS/OpenMP vars if not preset
    os.environ.setdefault("MKL_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
 
    # Configure logging level if requested
    if args.log_level:
        import logging
        numeric_level = getattr(logging, args.log_level.upper(), None)
        if isinstance(numeric_level, int):
            logger.setLevel(numeric_level)
 
    # Run the main function
    phased_graph, output_bam, output_vcf = realign_filter_per_cov(
        bam=args.bam,
        output_bam=args.output_bam,
        intrinsic_bam=args.intrinsic_bam,
        bam_region_bed=args.bam_region_bed,
        recall_mq_cutoff=args.recall_mq_cutoff,
        basequal_median_cutoff=args.basequal_median_cutoff,
        edge_weight_cutoff=args.edge_weight_cutoff,
        threads=args.threads,
        tmp_dir=args.tmp_dir,
        ref_genome=args.ref_genome,
        sample_id=args.sample_id,
        logger=logger
    )
    
    # Print results
    logger.info(f"Processing complete:")
    logger.info(f"  Output BAM: {output_bam}")
    logger.info(f"  Output VCF: {output_vcf}")
    if phased_graph:
        logger.info(f"  Graph saved: {args.bam.replace('.bam', '.phased.graphml')}")
    
    # Return appropriate exit code
    import sys
    if output_bam is None:
        logger.error("Failed to generate output BAM")
        sys.exit(1)
    else:
        sys.exit(0)



