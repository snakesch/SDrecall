import pysam
import numpy as np
from ncls import NCLS

from fp_control.numba_operators import fast_median, numba_sum
from src.log import logger
from src.utils import executeCmd
import shutil
import tempfile
import os


def overlapping_reads_iterator(ncls_dict, read_dict, chrom, start, end):
	"""
	Generator function to lazily yield overlapping read objects.

	Parameters:
	-----------
	ncls_dict : dict
		Dictionary of NCLS objects, keyed by chromosome names.
	read_dict : dict
		Dictionary of read objects, keyed by query name indices.
	qname_dict : dict
		Dictionary mapping query name indices to query names.
	chrom : str
		Chromosome name.
	start : int
		Start position of the interval (0-based, inclusive).
	end : int
		End position of the interval (0-based, exclusive).

	Yields:
	-------
	pysam.AlignedSegment
		Overlapping read object.

	Notes:
	------
	This function uses NCLS for efficient overlap queries and yields individual reads.
	It filters reads to ensure they actually overlap with the given interval.
	"""

	# Perform an overlap query on the NCLS tree
	if chrom not in ncls_dict or ncls_dict[chrom] is None:
		yield from ()
	else:
		_, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
		# overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
		for qname_idx in list(dict.fromkeys(overlapping_read_qnames)):
			for read in read_dict[qname_idx]:
				if read.reference_start < end and read.reference_end > start:
					yield read



def overlap_qname_idx_iterator(ncls_dict, chrom, start, end):
	"""
	Generator function to lazily yield query name indices of overlapping reads.

	Parameters:
	-----------
	ncls_dict : dict
		Dictionary of NCLS objects, keyed by chromosome names.
	chrom : str
		Chromosome name.
	start : int
		Start position of the interval (0-based, inclusive).
	end : int
		End position of the interval (0-based, exclusive).

	Yields:
	-------
	int
		Query name index of an overlapping read.

	Notes:
	------
	This function uses NCLS for efficient overlap queries.
	The start position is inclusive, and the end position is exclusive.
	"""

	# Perform an overlap query on the NCLS tree
	if chrom not in ncls_dict or ncls_dict[chrom] is None:
		yield from ()
	else:
		_, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
		# overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
		for qname_idx in dict.fromkeys(overlapping_read_qnames):
			yield qname_idx



def calculate_mean_read_length(bam_file_path, sample_size=100000):
	bamfile = pysam.AlignmentFile(bam_file_path, "rb")

	total_length = 0
	read_count = 0

	for read in bamfile:
		total_length += read.query_length
		read_count += 1
		if read_count >= sample_size:
			break

	bamfile.close()

	if read_count == 0:
		return 0  # Handle the case where there are no reads

	mean_length = total_length / read_count
	return mean_length



def is_read_noisy(read, paired, mapq_filter, basequal_median_filter=15, filter_noisy = True):
	"""Decide if a read should be treated as noisy.

	Rules:
	- Only primary alignments are evaluated for noisiness. Secondary or supplementary
	  alignments are treated as non-noisy and skipped.
	- filter_noisy controls ONLY quality/soft-clip based checks (preserves original behavior).
	- For other structural flags (unmapped, low MAPQ, etc.), the read is considered noisy
	  regardless of filter_noisy.
	"""
	# Only evaluate primary alignments for noise
	if read.is_secondary or read.is_supplementary:
		logger.debug(f"is_read_noisy: {read.query_name} skipped (secondary/supplementary alignment); not considered noisy")
		return False

	# Common fast checks for both paired and unpaired
	if read.is_unmapped:
		logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: unmapped read")
		return True

	if read.is_qcfail:
		logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: QC fail flag set")
		return True

	if read.mapping_quality < mapq_filter:
		logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: MAPQ {read.mapping_quality} < threshold {mapq_filter}")
		return True

	if read.reference_end is None:
		logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: reference_end is None")
		return True

	if read.query_sequence is None:
		logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: missing query_sequence")
		return True

	# Paired-end specific checks
	if paired:
		aln_len = (read.reference_end - read.reference_start) if (read.reference_end is not None and read.reference_start is not None) else 0
		if aln_len < 75:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: alignment span {aln_len} < 75")
			return True

		# Ensure mate on same reference for proper pairing in this pipeline
		if read.reference_name != read.next_reference_name:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: reference_name ({read.reference_name}) != next_reference_name ({read.next_reference_name})")
			return True

		if not read.is_proper_pair:
			logger.debug(f"is_read_noisy: warning:{read.query_name} is not a proper pair")
	else:
		# Single-end specific: drop optical/PCR duplicates if marked
		if read.is_duplicate:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: duplicate read (single-end mode)")
			return True

	# Base quality and soft-clip based checks (controlled by filter_noisy)
	if filter_noisy and read.query_qualities is not None:
		q = np.array(read.query_qualities, dtype=np.uint8)
		med_q = fast_median(q)
		if med_q <= basequal_median_filter:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: median baseQ {med_q} <= threshold {basequal_median_filter}")
			return True

		num_low = numba_sum(q < basequal_median_filter)
		if num_low >= 50:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: #bases with Q<{basequal_median_filter} is {num_low} >= 50")
			return True

		# Soft-clip length from CIGAR (op 4)
		softclip = 0
		if read.cigartuples is not None:
			softclip = numba_sum(np.array([t[1] for t in read.cigartuples if t[0] == 4], dtype=np.uint16))
		if softclip >= 20:
			logger.debug(f"is_read_noisy: {read.query_name} flagged noisy: total soft-clip length {softclip} >= 20")
			return True

	# If none of the noisy conditions triggered
	return False


def _collate_bam_file(bam_file, threads=4, logger=logger):
	"""Optionally collate BAM by qname using samtools collate.
	Returns path to collated BAM, or None if samtools unavailable or failed."""
	if shutil.which("samtools") is None:
		return None

	# Create a temp BAM file
	fd, tmp_path = tempfile.mkstemp(suffix=".bam")
	os.close(fd)
	cmd = f"samtools collate -f -@ {int(threads)} {bam_file} -o {tmp_path}"
	try:
		executeCmd(cmd, logger=logger)
		return tmp_path
	except Exception as e:
		logger.warning(f"samtools collate failed, using original BAM: {e}")
		# Cleanup if created
		try:
			os.remove(tmp_path)
		except Exception:
			pass
		return None


def _process_qname_group(qname, reads, read_dict, qname_interval_dict, 
						  qname_idx_dict, qname_dict, noisy_qnames, total_qnames,
						  qname_idx_counter, paired, mapq_filter, basequal_median_filter, 
						  filter_noisy, logger):
	"""
	Process all reads for a single qname, similar to Rust's process_qname_group.
	
	This function checks if ANY read in the group is noisy. If so, the entire
	qname is marked as noisy. Only records the qname if we have both read1 and read2
	for paired-end data.
	"""
	total_qnames.add(qname)
	
	# Check if any read in the group is noisy
	is_noisy = any(is_read_noisy(read, paired, mapq_filter, basequal_median_filter, filter_noisy) 
				   for read in reads)
	
	if is_noisy:
		logger.debug(f"This qname {qname} is noisy. Skip it.")
		noisy_qnames.add(qname)
		return qname_idx_counter
	
	# For paired-end mode, ensure we have both read1 and read2
	if paired:
		read1 = None
		read2 = None
		
		for read in reads:
			if read.is_read1:
				read1 = read
			elif read.is_read2:
				read2 = read
		
		# We need both reads for a complete pair
		if read1 is None or read2 is None:
			logger.debug(f"Incomplete pair for qname {qname}, skipping")
			return qname_idx_counter
	
	# Assign qname_idx and record reads
	if qname not in qname_idx_dict:
		qname_idx = qname_idx_counter
		qname_idx_dict[qname] = qname_idx
		qname_dict[qname_idx] = qname
		qname_idx_counter += 1
	else:
		qname_idx = qname_idx_dict[qname]
	
	# Store all reads for this qname
	read_dict[qname_idx] = reads.copy()
	
	# Update intervals for all reads
	for read in reads:
		chrom = read.reference_name
		start = read.reference_start
		end = read.reference_end
		
		if chrom is None or start is None or end is None:
			continue
		
		if qname_idx in qname_interval_dict[chrom]:
			interval = qname_interval_dict[chrom][qname_idx]
			updated_interval = (min(start, int(interval[0])), max(end, int(interval[1])))
		else:
			updated_interval = (int(start), int(end))
		qname_interval_dict[chrom][qname_idx] = updated_interval
	
	return qname_idx_counter


def migrate_bam_to_ncls(bam_file,
				mapq_filter = 25,
				basequal_median_filter = 20,
				paired = True,
				filter_noisy = True,
				logger = logger):
	"""
	Migrate BAM file data to NCLS (Nested Containment List) format for efficient interval querying.

	This function filters reads based on quality criteria, and organizes the data
	into NCLS structures for fast overlap queries.

	It completely put all alignment info into memory, thus faster than pysam disk IO based query.

	Parameters:
	-----------
	bam_file : str
		Path to the input BAM file.
	mapq_filter : int, optional
		Minimum mapping quality threshold for reads (default is 25).
	basequal_median_filter : int, optional
		Minimum median base quality threshold for reads (default is 20).
	paired : bool, optional
		If True, process as paired-end data; if False, as single-end (default is True).
	logger : logging.Logger, optional
		Logger object for output messages.

	Returns:
	--------
	tuple
		A tuple containing:
		- ncls_dict : dict
			Dictionary of NCLS objects, keyed by chromosome names.
		- read_dict : dict
			Dictionary of read objects, keyed by query name indices.
		- qname_dict : dict
			Dictionary mapping query name indices to query names. (qname_idx -> qname) The qname_idx is determined by the growing variable n with the iteration of all reads
		- qname_idx_dict : dict
			Dictionary mapping query names to query name indices. (qname -> qname_idx)
		- noisy_qnames : set
			Set of query names identified as noisy and filtered out.

	Within an NCLS object:
	(start, end) ---ncls_dict[chrom]---> qname_idx(interval_idx) ---qname_dict---> qname
	qname_idx(interval_idx) ---read_dict---> read objects(pysam.AlignedSegment)
	qname ---qname_idx_dict---> qname_idx(interval_idx) (which is also the interval index)

	Notes:
	------
	- The function uses pysam for BAM file processing and NCLS for interval data structure.
	- Reads are filtered based on mapping quality, median base quality, and other criteria.
	- For paired-end data, additional checks ensure proper pairing and orientation.
	- The NCLS structure allows for O(log n + m) time complexity for overlap queries, where n is the
	  number of intervals and m is the number of overlaps found.

	Example:
	--------
	>>> bam_file = "path/to/aligned_reads.bam"
	>>> ncls_dict, read_dict, qname_dict, qname_idx_dict, noisy_qnames = migrate_bam_to_ncls(bam_file)
	>>> # Now you can use ncls_dict for efficient overlap queries
	>>> chrom = "chr1"
	>>> start, end = 1000, 2000
	>>> overlapping_interval_indices = ncls_dict[chrom].find_overlap(start, end)
	>>> overlapping_reads = [read_dict[idx] for idx in overlapping_interval_indices]

	See Also:
	---------
	pysam.AlignmentFile : For BAM file handling
	ncls.NCLS : Nested Containment List Structure for interval queries
	"""

	# Collate input BAM by qname if samtools is available
	logger.info(f"Starting BAM processing: file={bam_file}, mapq_filter={mapq_filter}, basequal_median_filter={basequal_median_filter}, filter_noisy={filter_noisy}, paired={paired}")
	collated_path = _collate_bam_file(bam_file, threads=4, logger=logger) if paired else None
	use_bam = collated_path or bam_file

	with pysam.AlignmentFile(use_bam, "rb") as bam:
		chroms = bam.references
		logger.info(f"BAM file opened successfully, {len(chroms)} chromosomes found")
		read_dict = {}
		qname_interval_dict = {chrom: {} for chrom in chroms}
		ncls_dict = {chrom: None for chrom in chroms}
		qname_idx_dict = {}
		qname_dict = {}

		qname_idx_counter = 0
		noisy_qnames = set()
		total_qnames = set()
		total_reads_processed = 0
		skipped_alignments = 0

		# Buffer for collecting reads with the same qname (like Rust implementation)
		current_qname = None
		current_reads = []

		logger.info("Starting to process BAM records")
		for read in bam:
			total_reads_processed += 1
			# Skip secondary/supplementary/duplicate early
			if read.is_secondary or read.is_supplementary or read.is_duplicate:
				skipped_alignments += 1
				continue

			qname = read.query_name
			
			# Check if we've moved to a new qname
			if current_qname != qname:
				# Process the previous qname group if any
				if current_qname is not None and current_reads:
					qname_idx_counter = _process_qname_group(
						current_qname,
						current_reads,
						read_dict,
						qname_interval_dict,
						qname_idx_dict,
						qname_dict,
						noisy_qnames,
						total_qnames,
						qname_idx_counter,
						paired,
						mapq_filter,
						basequal_median_filter,
						filter_noisy,
						logger
					)
				
				current_qname = qname
				current_reads = []
			
			# Skip if this qname is already marked as noisy
			if qname in noisy_qnames:
				current_reads = []
				continue
			
			current_reads.append(read)

		# Process the last qname group
		if current_qname is not None and current_reads:
			_process_qname_group(
				current_qname,
				current_reads,
				read_dict,
				qname_interval_dict,
				qname_idx_dict,
				qname_dict,
				noisy_qnames,
				total_qnames,
				qname_idx_counter,
				paired,
				mapq_filter,
				basequal_median_filter,
				filter_noisy,
				logger
			)

	for chrom in chroms:
		qname_intervals = qname_interval_dict[chrom]
		if len(qname_intervals) > 0:
			qname_indices = np.fromiter(qname_intervals.keys(), dtype=np.int32)
			starts = np.array([int(qname_intervals[qname_idx][0]) for qname_idx in qname_indices])
			ends = np.array([int(qname_intervals[qname_idx][1]) for qname_idx in qname_indices])
			ncls_dict[chrom] = NCLS(starts, ends, qname_indices)

	logger.info(f"BAM processing complete: {total_reads_processed} total reads processed, {skipped_alignments} alignments skipped, {len(qname_idx_dict)} qnames retained, {len(noisy_qnames)} noisy qnames filtered")
	logger.info(f"Containing {len(qname_dict)} qnames in NCLS, left out {len(noisy_qnames)} noisy qnames, totally went through {len(total_qnames)} qnames")
	logger.info(f"Containing {len(qname_idx_dict)} key-value pairs in qname_idx_dict, the largest qname_idx is {max(qname_idx_dict.values()) if len(qname_idx_dict)>0 else -1}")

	if len(total_qnames) <= 2:
		return None

	# Cleanup temp collated file
	if collated_path:
		try:
			os.remove(collated_path)
		except Exception:
			pass

	return ncls_dict, read_dict, qname_dict, qname_idx_dict, noisy_qnames