import os
import pysam
import logging
import traceback

from src.utils import executeCmd, prepare_tmp_file
from src.const import shell_utils
from src.log import logger


def imap_filter_out(args, log_dir=""):
    import sys
    """Worker function that processes a single region and returns results with log file path"""
    # Unpack arguments
    raw_bam, output_bam, intrinsic_bam, bam_region_bed, max_varno, mapq_cutoff, basequal_median_cutoff, edge_weight_cutoff, numba_threads, job_id = args
    import numba
    numba.set_num_threads(numba_threads)
    numba.config.THREADING_LAYER = 'omp'

    # Create unique log file for this subprocess
    log_file = os.path.join(log_dir, f"subprocess_{job_id}_{os.path.basename(raw_bam)}.log")
    
    # Configure file logger
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter("[%(asctime)s] [%(levelname)s] [%(pathname)s:%(funcName)s:%(lineno)d] %(message)s"))
    
    # Create logger that writes to file
    subprocess_logger = logging.getLogger(f"SubProcess-{job_id}")
    subprocess_logger.setLevel(logging.INFO)
    subprocess_logger.addHandler(file_handler)
    subprocess_logger.propagate = False  # Don't send logs to parent loggers
    
    try:
        # Log the start of processing with all parameters for debugging
        subprocess_logger.info(f"Processing region in {raw_bam}")
        subprocess_logger.info(f"Parameters: output_bam={output_bam}, intrinsic_bam={intrinsic_bam}, "
                               f"bam_region_bed={bam_region_bed}, max_varno={max_varno}, "
                               f"mapq_cutoff={mapq_cutoff}, basequal_median_cutoff={basequal_median_cutoff}, "
                               f"edge_weight_cutoff={edge_weight_cutoff}")
        
        # Call the main processing function
        phased_graph, output_bam = realign_filter_per_cov(
            bam=raw_bam,
            output_bam=output_bam,
            intrinsic_bam=intrinsic_bam,
            bam_region_bed=bam_region_bed,
            max_varno=max_varno,
            mapq_cutoff=mapq_cutoff,
            basequal_median_cutoff=basequal_median_cutoff,
            edge_weight_cutoff=edge_weight_cutoff,
            logger=subprocess_logger
        )
        
        # If successful, return result with comma-separated fields and log file path
        subprocess_logger.info(f"Successfully processed {raw_bam}")
        return f"{raw_bam},{output_bam},{log_file}"
        
    except Exception as e:
        # Log the full stack trace for debugging
        subprocess_logger.error(f"Error processing {raw_bam}: {str(e)}")
        subprocess_logger.error(traceback.format_exc())
        # Return error result with log file path
        return f"{raw_bam},NaN,{log_file}"


def realign_filter_per_cov(bam,
                           output_bam = None,
                           intrinsic_bam = None,
                           bam_region_bed = None,
                           max_varno = 5,
                           mapq_cutoff = 20,
                           basequal_median_cutoff = 15,
                           edge_weight_cutoff = 0.201,
                           threads = 4,
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
    - max_varno (float): Maximum variant number allowed
    - mapq_cutoff (int): Mapping quality cutoff to be included in the analysis
    - basequal_median_cutoff (int): Base quality median cutoff to be included in the analysis (if the median BQ of a read is lower than this cutoff, we will discard this read because it is too noisy)
    - edge_weight_cutoff (float): Edge weight cutoff separating two rounds of BK clique searches
    - logger (Logger object): Logger for output messages

    Returns:
    - phased_graph (Graph object): The constructed phasing graph

    This function integrates various analysis steps including BAM processing,
    graph construction, haplotype identification, and read filtering to improve
    the quality of genomic alignments and identify potential misalignments.
    '''

    from fp_control.bam_ncls import migrate_bam_to_ncls, calculate_mean_read_length
    from fp_control.graph_build import build_phasing_graph
    from fp_control.identify_misaligned_haps import inspect_by_haplotypes
    from fp_control.phasing import phasing_realigned_reads

    # Given the top 1% mismatch count per read (one structrual variant count as 1 mismatch)
    tmp_bam = prepare_tmp_file(suffix=".bam").name

    max_varno = float(max_varno)

    noisy_qnames = set([])
    mismap_qnames = set([])
    norm_qnames = {"":set([])}
    noisy_num = 0
    total_num = 0

    bam_ncls = migrate_bam_to_ncls(bam,
                                   mapq_filter = mapq_cutoff,
                                   basequal_median_filter = basequal_median_cutoff,
                                   logger=logger)

    # parse the results from the tuple returned by migrate_bam_to_ncls
    ncls_dict, read_dict, qname_dict, qname_idx_dict, total_lowqual_qnames = bam_ncls
    logger.info(f"Successfully migrated the BAM file {bam} to NCLS format\n\n")

    # Now migrate the intrinsic BAM file to NCLS format
    intrin_bam_ncls = migrate_bam_to_ncls(intrinsic_bam,
                                          mapq_filter = 0,
                                          basequal_median_filter = 0,
                                          paired = False,
                                          logger=logger)
    # Since intrinsic BAM reads are reference sequences, therefore there are no low quality reads
    intrin_bam_ncls = intrin_bam_ncls[:-1]
    logger.info(f"Successfully migrated the intrinsic BAM file {intrinsic_bam} to NCLS format\n")
    logger.info(f"Containing {len(intrin_bam_ncls[1])} reads in total.\n\n")
    bam_graph = bam.replace(".bam", ".phased.graphml")

    # Calculate the mean read length of the input bam file, which can be used for read pair similarity calculation
    mean_read_length = calculate_mean_read_length(bam)

    # Create the read-pair graph used for phasing
    # Detailed description of the graph construction can be found in the function docstring.
    phased_graph, weight_matrix, qname_to_node, total_readhap_vector, read_ref_pos_dict = build_phasing_graph(bam,
                                                                                                              ncls_dict,
                                                                                                              read_dict,
                                                                                                              qname_dict,
                                                                                                              mean_read_length,
                                                                                                              logger = logger)
    if phased_graph is None:
        return None

    logger.info(f"Now succesfully built the phasing graph with {phased_graph.num_vertices()} vertices and {phased_graph.num_edges()} edges. Save it to {bam_graph}\n\n")
    # Now we need to extract the components in the phased graph
    phased_graph.save(bam_graph)

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
    total_readerr_vector = {}
    compare_haplotype_meta_tab = bam.replace(".bam", ".haplotype_meta.tsv")

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
                        hap_id = f"{hap_id}_LOWQUAL"
                    # Filter out reads with oddly high editing distance that breakthrough the cutoff
                    if read.is_mapped and read.mapping_quality > 10:
                        gap_sizes = [t[1] for t in read.cigartuples if t[0] in [1,2] and t[1] > 1]
                        max_gap_size = max(gap_sizes) if len(gap_sizes) > 0 else 0
                        edit_dist = read.get_tag("NM")
                        scatter_edit_dist = edit_dist - max_gap_size
                        if qname in total_lowqual_qnames:
                            pass
                        elif qname in noisy_qnames:
                            # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is too high to be likely correctly mapped")
                            pass
                        elif qname in mismap_qnames:
                            # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is significantly higher than other reads so it is considered as noisy and misaligned")
                            pass
                        elif scatter_edit_dist > max_varno and qname not in correct_qnames:
                            noisy_num += 1
                            logger.info(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is so high that it is considered as noisy and misaligned")
                            noisy_qnames.add(qname)
                        elif not read.is_secondary and \
                            not read.is_supplementary and \
                            not read.is_duplicate and \
                            not read.is_qcfail:
                            # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {read.get_tag('NM')} which is considered as normal")
                            qname_pair = norm_qnames.get(qname, set([]))
                            qname_pair.add(read)
                            norm_qnames[qname] = qname_pair
                    else:
                        logger.info(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is not mapped or has low mapping quality {read.mapping_quality}, skip this read")

                    if qname in noisy_qnames:
                        hap_id = f"{hap_id}_HIGHVD"
                    read.set_tag('HP', f'HAP_{hap_id}')
                    tmp_handle.write(read)

            # logger.warning(f"Check {check_odd_num} reads to find oddly high editing distance reads, {zero_odd_num} reads found no odd editing distance read pairs")
            for qname, pair in norm_qnames.items():
                if (not qname in noisy_qnames) and (not qname in mismap_qnames):
                    for read in pair:
                        output_handle.write(read)

    logger.warning(f"Filtered out {len(noisy_qnames)} noisy read-pairs (Editing distance without the biggest gap > {max_varno}) and {len(mismap_qnames - noisy_qnames)} read-pairs with ODD high editing distance, remaining {len(set(norm_qnames.keys()) - noisy_qnames - mismap_qnames)} read-pairs from {bam} (with total {total_num} reads) and output to {output_bam}\n\n")

    # Replace the input BAM file with the tmp BAM file with modified RG tags for visualization of haplotype clusters
    executeCmd(f"samtools sort -O bam -o {bam} {tmp_bam} && \
                 samtools index {bam} && \
                 rm {tmp_bam}", logger=logger)
    # Sort and index the output bam file
    executeCmd(f"samtools sort -O bam -o {tmp_bam} {output_bam} && \
                 mv {tmp_bam} {output_bam} && \
                 samtools index {output_bam}", logger=logger)
    return phased_graph, output_bam







