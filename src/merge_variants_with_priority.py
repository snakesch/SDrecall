import pysam
import argparse
import multiprocessing
import subprocess
import logging
import uuid
from io import StringIO
import inspect
import time
import traceback
import sys

from src.log import logger, log_command
from src.const import shell_utils
from src.utils import executeCmd


class VariantRecordWrapper:
	def __init__(self, record):
		self.record = record
		
	def __getstate__(self):
		return self.record.to_dict()

	def __setstate__(self, state):
		self.record = pysam.VariantRecord.from_dict(state)

	def __hash__(self):
		# Define a custom hash function based on relevant attributes
		return hash((self.record.chrom, self.record.pos, tuple(self.record.alleles)))

	def __eq__(self, other):
		# Define equality comparison based on relevant attributes
		return (self.record.chrom == other.record.chrom and
				self.record.pos == other.record.pos and
				self.record.ref == other.record.ref and 
				self.record.alts == other.record.alts)
		
	def __repr__(self):
		return f"{self.record}"
	
	def pickable(self):
		return (self.record.chrom, 
				self.record.pos, 
				self.record.ref, 
				self.record.alts, 
				self.record.id, 
				self.record.qual, 
				tuple(list(self.record.filter)),
				tuple([(k, v) for k, v in self.record.info.items()]),
				tuple([(k, tuple(v.items())) for k, v in self.record.samples.items()]))
		

def extract_ad_info(record):
	# Extract the AD information from the record
	ad_info = record.samples[0]['AD']
	return ad_info


def modify_gt_based_on_ad_gq(record, rrecord):
	# Iterate over each sample in the record
	for sample_name in record.samples:
		rsample_info = rrecord.samples[sample_name]
		sample_info = record.samples[sample_name]
		# Extract AD and GQ information
		ref, alt = sample_info.get('AD', [0, 0]) # Default to [0, 0] if AD is not present
		gq = sample_info.get('GQ', 0) # Default to 0 if GQ is not present

		rref, ralt = rsample_info.get('AD', [0, 0]) # Default to [0, 0] if AD is not present
		rgq = rsample_info.get('GQ', 0) # Default to 0 if GQ is not present

		hps = rsample_info.get('HPSUP', ".") # Default to "." if HPSUP is not present
		# logger.info(f"The hps returned by get is {hps} and its type is {type(hps)}")
		num_hps = len(hps[0].split(";"))

		if rrecord.samples[sample_name]['GT'] == (1, 1):
			continue

		if ralt/(ralt + rref) >= 0.9:
			rrecord.samples[sample_name]['GT'] = (1, 1)
			logger.info(f"Setting the GT to (1, 1) due to the alt/dp ratio {ralt/(ralt + rref)} for the variant {record.chrom}:{record.pos}:{record.ref} -> {record.alts}")
			continue

		if num_hps >= 2 and ralt/(ralt + rref) >= 0.33:
			rrecord.samples[sample_name]['GT'] = (1, 1)
			logger.info(f"Setting the GT to (1, 1) due to the num_hps {num_hps} and alt/dp ratio {ralt/(ralt + rref)} for the variant {record.chrom}:{record.pos}:{record.ref} -> {record.alts}")
			continue

		if num_hps >= 3 and ralt/(ralt + rref) >= 0.30:
			rrecord.samples[sample_name]['GT'] = (1, 1)
			logger.info(f"Setting the GT to (1, 1) due to the num_hps {num_hps} and alt/dp ratio {ralt/(ralt + rref)} for the variant {record.chrom}:{record.pos}:{record.ref} -> {record.alts}")
			continue

		if num_hps >= 4 and ralt/(ralt + rref) >= 0.25:
			rrecord.samples[sample_name]['GT'] = (1, 1)
			logger.info(f"Setting the GT to (1, 1) due to the num_hps {num_hps} and alt/dp ratio {ralt/(ralt + rref)} for the variant {record.chrom}:{record.pos}:{record.ref} -> {record.alts}")
			continue

		if rgq < 5 and ralt/(ralt + rref) >= 0.5 and (ralt + rref) > 5:
			rrecord.samples[sample_name]['GT'] = (1, 1)
			logger.info(f"Setting the GT to (1, 1) due to the gq {rgq} and alt/dp ratio {ralt/(ralt + rref)} for the variant {record.chrom}:{record.pos}:{record.ref} -> {record.alts}")
			continue

	return rrecord


def merge_same_variant_rec( qrecord, 
							rrecord, 
							qv_tag=None,
							rv_tag=None,
							modify_gt = True ):
	qrecord = qrecord.record
	rrecord = rrecord.record

	if qv_tag:
		if qv_tag not in rrecord.filter:
			rrecord.filter.add(qv_tag)
	
	if rv_tag:
		if rv_tag not in rrecord.filter:
			rrecord.filter.add(rv_tag)

	if modify_gt:
		mrecord = modify_gt_based_on_ad_gq(qrecord, rrecord)
	else:
		mrecord = rrecord
	return VariantRecordWrapper(mrecord)



def check_empty_iterator(iterator):
	try:
		next(iterator)
	except StopIteration:
		return True
	else:
		return False
	
	
def compare_variant_positions(variant1, variant2):
	# -1 means not in the same chromosome
	# 0 means exact the same location
	# 1 means variant1 comes before variant2
	# 2 means variant1 comes after variant2
	if variant1.contig != variant2.contig:
		return -1
	if variant1.start < variant2.start:
		return 1
	elif variant1.start > variant2.start:
		return 2
	else:
		if variant1.stop < variant2.stop:
			return 1
		elif variant1.stop > variant2.stop:
			return 2
		else:
			return 0


def deal_with_same_loc_variants(var1: VariantRecordWrapper,
								var2: VariantRecordWrapper,
								var1_iterator, 
								var2_iterator,
								buffer_var1: list,
								buffer_var2: list,
								non_overlap_var1: set,
								non_overlap_var2: set,
								merged_records: set,
								merging_func,
								logger = logger, 
								**kwargs):
	# Considering the possibility that there might be multiple variant records at the same location
	# Solution is: iterate both iterators to extract all the records at the same location and compare two lists to find out matched records.
	current_loc_1recs = [var1]
	current_loc_2recs = [var2]
	# Iterate through complement_iterator for a few times until hitting the downstream variant records or bottom
	logger.debug(f"Run into a situation where the two comparing records are at the same location but with different alleles: variant record 1: {var1}, variant record 2: {var2}")
	while True:
		try:
			next_var1 = VariantRecordWrapper(next(var1_iterator))
			logger.debug(f"Trying to collect the next variant record 1 from iterator: {next_var1}")
		except StopIteration:
			# Hit the bottom of the query iterator but still captured all the records at this position to buffer_crecs and current_loc_crecs
			break
		else:
			if compare_variant_positions(var1.record, next_var1.record) == 0:
				current_loc_1recs.append(next_var1)
			elif compare_variant_positions(var1.record, next_var1.record) == 1:
				buffer_var1.append(next_var1)
				break
			else:
				raise ValueError(f'The variant 1 records are not sorted by position: \n{buffer_var1}\n')
			
	# Iterate through priority_iterator for a few times until hitting the downstream variant records or bottom
	while True:
		try:
			next_var2 = VariantRecordWrapper(next(var2_iterator))
			logger.debug(f"Trying to collect the next variant 2 record from iterator: {next_var2}")
		except StopIteration:
			# Hit the bottom of the reference iterator but still captured all the records at this position to buffer_precs and current_loc_precs
			break
		else:
			if compare_variant_positions(var2.record, next_var2.record) == 0:
				current_loc_2recs.append(next_var2)
			elif compare_variant_positions(var2.record, next_var2.record) == 1:
				buffer_var2.append(next_var2)
				break
			else:
				raise ValueError(f'The reference records are not sorted by position: \n{buffer_var2}\n')
			
	# Now we need to note that buffer_precs and buffer_crecs might contain the downstream variant records
	# But the current_loc_precs and current_loc_crecs are the records at the same location
	# So we need to compare the records between the two lists to find the matched records
	matching_var2 = []
	for v1 in current_loc_1recs:
		found_match = False
		for v2 in current_loc_2recs:
			if v1 == v2:
				matching_var2.append(v2)
				# Perform the binomial test
				merged_record = merging_func(v1, v2, **kwargs)
				merged_records.add(merged_record)
				found_match = True
				break
		if not found_match:
			# If no matching reference record is found, add the query record as is
			non_overlap_var1.add(v1)
	non_overlap_var2.update(set([v2 for v2 in current_loc_2recs if v2 not in matching_var2]))

	return merged_records, non_overlap_var1, non_overlap_var2, buffer_var1, buffer_var2



def imap_process_region(args):
	return process_region(*args)


@log_command
def process_region(region, 
				   query_vcf, 
				   reference_vcf, 
				   added_filter, 
				   qv_tag,
				   rv_tag,
				   modify_gt = True,
				   logger=logger):
	'''
	The function is not only used to return the modified overlapping records
	but also to return the records that are not overlapping with the reference VCF file
	'''
	
	query_vcf = pysam.VariantFile(query_vcf)
	reference_vcf = pysam.VariantFile(reference_vcf)
	
	# Fetch the records in the specified region from both VCF files
	query_iterator = query_vcf.fetch(region, reopen=True)
	reference_iterator = reference_vcf.fetch(region, reopen=True)

	# Process edge cases when eighter of the iterator returns no records
	if check_empty_iterator(query_iterator):
		return (frozenset([]), frozenset([]), frozenset([ VariantRecordWrapper(r).pickable() for r in reference_iterator ]))
	if check_empty_iterator(reference_iterator):
		query_iterator = query_vcf.fetch(region)
		return (frozenset([]), frozenset([ VariantRecordWrapper(r).pickable() for r in query_iterator ]), frozenset([]))
	
	# Now I ensure both iterators are not empty
	if added_filter:
		query_vcf.header.add_meta('FILTER', items=[('ID', added_filter), ('Description', f'Variant is likely to be {added_filter}')])
		reference_vcf.header.add_meta('FILTER', items=[('ID', added_filter), ('Description', f'Variant is likely to be {added_filter}')])
	if qv_tag:
		query_vcf.header.add_meta('FILTER', items=[('ID', qv_tag), ('Description', f'Variant is from {qv_tag}')])
		reference_vcf.header.add_meta('FILTER', items=[('ID', qv_tag), ('Description', f'Variant is from {qv_tag}')])
	if rv_tag:
		query_vcf.header.add_meta('FILTER', items=[('ID', rv_tag), ('Description', f'Variant is from {rv_tag}')])
		reference_vcf.header.add_meta('FILTER', items=[('ID', rv_tag), ('Description', f'Variant is from {rv_tag}')])
		
	query_iterator = query_vcf.fetch(region)
	reference_iterator = reference_vcf.fetch(region)
	
	# Initialize variables
	merged_records = set([])
	buffer_rrecs = []
	buffer_qrecs = []
	non_overlap_qrecs = set([])
	non_overlap_rrecs = set([])

	reference_record = None
	
	# Process records on the fly
	while True:
		try:
			if len(buffer_qrecs) > 0:
				query_record = buffer_qrecs.pop(0)
			else:
				query_record = VariantRecordWrapper(next(query_iterator))
		except StopIteration:
			# If we hit the bottom of the query iterator, we need to clear the buffer
			if reference_record is not None:
				non_overlap_rrecs.add(reference_record)
			
			if len(buffer_rrecs) > 0:
				non_overlap_rrecs.update(set(buffer_rrecs))
				
			while True:
				try:
					reference_record = VariantRecordWrapper(next(reference_iterator))
					logger.debug(f"Getting the next priority variant record from iterator: {reference_record}")
				except StopIteration:
					# If we hit the bottom of the reference iterator, we need to clear the buffer
					break
				else:
					non_overlap_rrecs.add(reference_record)
					
			break
		else:
			if reference_record is None:
				pass
			elif compare_variant_positions(query_record.record, reference_record.record) == 2:
				non_overlap_rrecs.add(reference_record)
			elif compare_variant_positions(query_record.record, reference_record.record) == 1:
				non_overlap_qrecs.add(query_record)
				continue
			elif compare_variant_positions(query_record.record, reference_record.record) == -1:
				raise ValueError('The chromosome of the query record and the reference record are not the same: {} vs {}'.format(query_record, reference_record))
			else:
				# The two records are at the same location
				merged_records, non_overlap_qrecs, non_overlap_rrecs, buffer_qrecs, buffer_rrecs = deal_with_same_loc_variants( query_record,
																																reference_record,
																																query_iterator,
																																reference_iterator,
																																buffer_qrecs,
																																buffer_rrecs,
																																non_overlap_qrecs,
																																non_overlap_rrecs,
																																merged_records,
																																merge_same_variant_rec,
																																modify_gt = modify_gt,
																																logger = logger,
																																qv_tag=qv_tag,
																																rv_tag=rv_tag)
				reference_record = None
				continue
			while True:
				try:
					if len(buffer_rrecs) > 0:
						reference_record = buffer_rrecs.pop(0)
					else:
						reference_record = VariantRecordWrapper(next(reference_iterator))
				except StopIteration:
					# If we hit the bottom of the reference iterator, we need to clear the buffer
					non_overlap_qrecs.add(query_record)
					break
				else:
					# Find the matching reference record
					if compare_variant_positions(query_record.record, reference_record.record) == 2:
						# Meaning the reference_record is upstream of the query_record, then we need to continue the inner while_loop
						non_overlap_rrecs.add(reference_record)
						continue
					elif compare_variant_positions(query_record.record, reference_record.record) == 1:
						# Meaning the reference_record is downstream of the query_record, then we need to break the inner while_loop and let the outer while_loop to move to the next variant record
						non_overlap_qrecs.add(query_record)
						break
					elif compare_variant_positions(query_record.record, reference_record.record) == -1:
						# This shouldn't happen because we are fetching records from the same chromosome region
						raise ValueError('The chromosome of the query record and the reference record are not the same: {} vs {}'.format(query_record, reference_record))
					else:
						# The two records are at the same location
						merged_records, non_overlap_qrecs, non_overlap_rrecs, buffer_qrecs, buffer_rrecs = deal_with_same_loc_variants( query_record, 
																																		reference_record, 
																																		query_iterator, 
																																		reference_iterator, 
																																		buffer_qrecs, 
																																		buffer_rrecs, 
																																		non_overlap_qrecs, 
																																		non_overlap_rrecs, 
																																		merged_records, 
																																		merge_same_variant_rec,
																																		logger = logger,
																																		qv_tag=qv_tag,
																																		rv_tag=rv_tag,
																																		modify_gt = modify_gt)
						reference_record = None
						break
	return (frozenset([ r.pickable() for r in merged_records ]), 
			frozenset([ r.pickable() for r in non_overlap_qrecs ]), 
			frozenset([ r.pickable() for r in non_overlap_rrecs ]))



def sort_vcf(vcf_file, ref_genome, output_vcf = None, logger = logger):
	# Sort the VCF file using bcftools
	if output_vcf is None:
		output_vcf = vcf_file.replace('.vcf', '.sorted.vcf')
	
	cmd = f'''bcftools norm -m -both -f {ref_genome} --multi-overlaps 0 -a -Ou {vcf_file} | \
			  bcftools norm -d exact -Ou - | \
			  bcftools filter -e 'ALT[0] == "*" || GT = "mis" || GT = "ref"' -Ou - | \
			  bcftools sort -Oz -o {output_vcf} - && \
			  bcftools index -f -t {output_vcf}'''
	executeCmd(cmd, logger=logger)
	return output_vcf



def merge_vcf_headers(header1, header2):
	merged_header = header1.copy()

	# Merge FILTER header records
	for record in header2.filters.values():
		if record.name not in merged_header.filters:
			merged_header.add_record(record.record)

	# Merge INFO header records
	for record in header2.info.values():
		if record.name not in merged_header.info:
			merged_header.add_record(record.record)

	# Merge FORMAT header records
	for record in header2.formats.values():
		if record.name not in merged_header.formats:
			merged_header.add_record(record.record)

	# Merge contig (SQ) header records
	for record in header2.contigs.values():
		if record.name not in merged_header.contigs:
			merged_header.add_record(record.record)

	# Merge other header records
	for record in header2.records:
		if record.key not in merged_header.records:
			merged_header.add_record(record)
			
	logger.info(f"The merged header looks like this: \n{merged_header}\n")

	return merged_header



def merge_with_priority(query_vcf = "", 
						reference_vcf = "", 
						output_vcf = "", 
						added_filter = None,
						qv_tag = None,
						rv_tag = None, 
						ref_genome = "",
						modify_gt = True,
						threads = 4):
	# Sort the query VCF file
	sorted_query_vcf = sort_vcf(query_vcf, ref_genome, logger = logger)

	# Sort the reference VCF file
	sorted_reference_vcf = sort_vcf(reference_vcf, ref_genome, logger = logger)

	# Open the sorted query VCF file
	bcf_query = pysam.VariantFile(sorted_query_vcf)
	bcf_reference = pysam.VariantFile(sorted_reference_vcf)

	# Open the output VCF file for writing
	tmp_output_vcf = output_vcf.replace(".vcf", ".tmp.vcf")
	# Add new FILTER tags to the header
	if added_filter:
		bcf_query.header.add_meta('FILTER', items=[('ID', added_filter), ('Description', f'Variant is likely to be {added_filter}')])
	if qv_tag:
		bcf_query.header.add_meta('FILTER', items=[('ID', qv_tag), ('Description', f'Variant is from {qv_tag}')])
	if rv_tag:
		bcf_query.header.add_meta('FILTER', items=[('ID', rv_tag), ('Description', f'Variant is from {rv_tag}')])

	merged_header = merge_vcf_headers(bcf_reference.header, bcf_query.header)
	bcf_output = pysam.VariantFile(tmp_output_vcf, 'w', header=merged_header)

	# Get the list of regions to process
	regions = list(dict.fromkeys(list(bcf_query.header.contigs.keys()) + list(bcf_reference.header.contigs.keys())))
	
	logger.info(f"Whether modify the GT field: {modify_gt}")
	logger.info(f"The query VCF file has records across {len(regions)} chromosomes, they are: {regions}")

	# Create a multiprocessing pool
	pool = multiprocessing.Pool(threads)

	# Process each region in parallel
	result_iterator = pool.imap_unordered(imap_process_region, [(region, 
																 sorted_query_vcf, 
																 sorted_reference_vcf, 
																 added_filter,
																 qv_tag,
																 rv_tag,
																 modify_gt) for region in regions])
	i=0
	for success, result, log_contents in result_iterator:
		i+=1
		print(f"\n************************************{i}_subprocess_start_for_compare_regional_variants************************************\n", file=sys.stderr)
		if success:
			merged_recs, non_overlap_qrecs, non_overlap_rrecs = result
			for rec_tup in merged_recs:
				chrom, pos, ref, alts, vid, qual, filters, info, samples = rec_tup
				if rv_tag:
					filters = [f for f in filters if f != rv_tag] + [rv_tag]
				if qv_tag:
					filters = [f for f in filters if f != qv_tag] + [qv_tag]
				rec = bcf_output.new_record( contig=chrom, 
											 start=pos-1, 
											 alleles=(ref,) + alts, 
											 id=vid, 
											 qual=qual, 
											 filter=filters, 
											 info=dict(info) )
				for sample, values in samples:
					for key, value in values:
						try:
							rec.samples[sample][key] = value
						except TypeError as te:
							logger.error(f"Failed to write value {value} to FORMAT field {key} at variant {chrom}:{pos}:{ref} -> {alts}, with error: {te}. Something wrong with the header definition or maybe your record in conventional VCF is missing genotype?")
							continue
				bcf_output.write(rec)
			for qrec_tup in non_overlap_qrecs:
				chrom, pos, ref, alts, vid, qual, filters, info, samples = qrec_tup
				if qv_tag:
					filters = [f for f in filters if f != qv_tag] + [qv_tag]
				if added_filter:
					filters = [f for f in filters if f != added_filter] + [added_filter]
				qrec = bcf_output.new_record(contig=chrom, 
											 start=pos-1, 
											 alleles=(ref,) + alts, 
											 id=vid, 
											 qual=qual, 
											 filter=filters, 
											 info=dict(info))
				for sample, values in samples:
					hps = [t[1] for t in values if t[0] == "HPSUP"]
					logger.debug(f"The hps returned by get is {hps} and its type is {type(hps)}")
					num_hps = len(hps[0][0].split(";")) if len(hps) > 0 else 0
					ad = [t[1] for t in values if t[0] == "AD"][0] # Default to [0, 0] if AD is not present
					if len(ad) > 2:
						logger.warning(f"This record is not biallelic, take a look at the record: {chrom}:{pos}:{ref} -> {alts}")
						ref_dp, alt_dp = ad[0], ad[1]
					else:
						ref_dp, alt_dp = ad
					
					gq = [t[1] for t in values if t[0] == "GQ"][0] # Default to 0 if GQ is not present
					gq = 0 if gq is None else gq
					ref_dp = 0 if ref_dp is None else ref_dp
					alt_dp = 0 if alt_dp is None else alt_dp
					for key, value in values:
						if key == "HPSUP":
							value = ";".join(value)
						if key == "GT":
							if num_hps >= 2 and alt_dp/(alt_dp + ref_dp) >= 0.55 and (alt_dp + ref_dp) >= 5 and modify_gt:
								logger.warning(f"Setting the GT to (1, 1) due to the num_hps {num_hps} and alt/dp ratio {alt_dp/(alt_dp + ref_dp)} for the variant {chrom}:{pos}:{ref} -> {alts}")
								value = (1, 1)
							if num_hps >= 4 and alt_dp >= ref_dp and (alt_dp + ref_dp) >= 5 and modify_gt:
								logger.warning(f"Setting the GT to (1, 1) due to the num_hps {num_hps} and alt/dp ratio {alt_dp/(alt_dp + ref_dp)} for the variant {chrom}:{pos}:{ref} -> {alts}")
								value = (1, 1)
							if alt_dp/(alt_dp + ref_dp) >= 0.9 and (alt_dp + ref_dp) >= 5 and modify_gt:
								logger.warning(f"Setting the GT to (1, 1) due to the alt/dp ratio {alt_dp/(alt_dp + ref_dp)} for the variant {chrom}:{pos}:{ref} -> {alts}")
								value = (1, 1)
						try:
							qrec.samples[sample][key] = value
						except TypeError as te:
							logger.error(f"Failed to write value {value} to FORMAT field {key} at variant {chrom}:{pos}:{ref} -> {alts}, with error: {te}")
							continue
				bcf_output.write(qrec)
			for rrec_tup in non_overlap_rrecs:
				chrom, pos, ref, alts, vid, qual, filters, info, samples = rrec_tup
				if rv_tag:
					filters = [f for f in filters if f != rv_tag] + [rv_tag]
				rrec = bcf_output.new_record(contig=chrom, 
											 start=pos-1, 
											 alleles=(ref,) + alts, 
											 id=vid, 
											 qual=qual, 
											 filter=filters, 
											 info=dict(info))
				for sample, values in samples:
					ad = [t[1] for t in values if t[0] == "AD"][0] # Default to [0, 0] if AD is not present
					if len(ad) > 2:
						logger.warning(f"This record is not biallelic, take a look at the record: {chrom}:{pos}:{ref} -> {alts}")
						ref_dp, alt_dp = ad[0], ad[1]
					else:
						ref_dp, alt_dp = ad
					gq = [t[1] for t in values if t[0] == "GQ"][0] # Default to 0 if GQ is not present
					gq = 0 if gq is None else gq
					ref_dp = 0 if ref_dp is None else ref_dp
					alt_dp = 0 if alt_dp is None else alt_dp
					for key, value in values:
						if key == "GT":
							if gq < 5 and alt_dp/(alt_dp + ref_dp) >= 0.7 and modify_gt:
								value = (1, 1)
						try:
							rrec.samples[sample][key] = value
						except TypeError as te:
							logger.error(f"Failed to write value {value} to FORMAT field {key} at variant {chrom}:{pos}:{ref} -> {alts}. The error message is {te}")
							continue
				bcf_output.write(rrec)
			print(f"Successfully processed the variant records in one chromosome. The log info are:\n{log_contents}\n", file=sys.stderr)
		else:
			error_mes, tb_str = result
			raise ValueError(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")
		print(f"\n************************************{i}_subprocess_end_for_compare_regional_variants************************************\n", file=sys.stderr)

	pool.close()

	# Close the VCF files
	bcf_query.close()
	bcf_output.close()

	# Index the output VCF file
	cmd = f"bcftools norm -d exact {tmp_output_vcf} | \
			bcftools sort -Oz -o {output_vcf} - && \
			bcftools index {output_vcf} && \
			rm {tmp_output_vcf} && \
			ls -lh {output_vcf}"
	executeCmd(cmd)
	
	
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compare and filter VCF files')
	parser.add_argument("-qv", '--query_vcf', help='Path to the query VCF file')
	parser.add_argument("-rv", '--reference_vcf', help='Path to the reference VCF file')
	parser.add_argument("-ov", '--output_vcf', help='Path to the output VCF file')
	parser.add_argument("-fl", "--filter_tag", type=str, help="The filter tag you want to added to the matched query VCF records", required=False, default="")
	parser.add_argument("-qt", "--query_vcf_tag", type=str, help="The tag added to the FILTERs to show the source of the variant from query_vcf_file", required=False, default="RAW")
	parser.add_argument("-rt", "--reference_tag", type=str, help="The tag added to the FILTERs to show the source of the variant from reference_vcf_file", required=False, default="CLEAN")
	parser.add_argument("-rg", "--ref_genome", type=str, help="The reference genome file", required=False, default="")
	parser.add_argument("-mg", "--modify_gt", action="store_false", help="Whether to modify the GT field", required=False, default=True)
	parser.add_argument('--threads', type=int, default=4, help='Number of threads for parallel processing', required=False)

	args = parser.parse_args()

	merge_with_priority(query_vcf = args.query_vcf, 
						reference_vcf = args.reference_vcf, 
						output_vcf = args.output_vcf, 
						added_filter = args.filter_tag, 
						qv_tag = args.query_vcf_tag,
						rv_tag = args.reference_tag,
						ref_genome = args.ref_genome,
						modify_gt = args.modify_gt,
						threads = args.threads)