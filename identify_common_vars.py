import pysam
import logging
import os
import re
import sys

import pandas as pd
import numpy as np
import multiprocessing as mp

from scipy.stats import binom

from src.log import logger, log_command
from src.merge_variants_with_priority import sort_vcf, compare_variant_positions
from src.utils import executeCmd, prepare_tmp_file, na_value


def determine_common_per_pysam_record(record, 
                                      AC_tag="AC", 
                                      AN_tag="AN", 
                                      inhouse_common_cutoff=0.01,
                                      conf_level=0.999):
    """
    Determine if a small variant is common in the inhouse cohort using pysam
    
    Args:
        record: A pysam.VariantRecord object
        AC_tag: Tag for allele count info field (e.g., "AC")
        AN_tag: Tag for allele number info field (e.g., "AN")
        inhouse_common_cutoff: Frequency threshold for common variants
        conf_level: Confidence level threshold (default: 0.999)
        
    Returns:
        Boolean indicating if the variant is common
    """
    # Get AC and AN values from INFO fields
    try:
        AC = record.info.get(AC_tag, None)
        AN = record.info.get(AN_tag, None)
    except KeyError:
        logger.debug(f"Missing AC/AN fields for {record.contig}:{record.pos}")
        return False
    
    # Handle list type for AC
    if isinstance(AC, (list, tuple)):
        AC = AC[0]
    
    logger.debug(f"AC is {AC} and AN is {AN}")
    
    # Handle missing values
    if na_value(AC) or na_value(AN):
        logger.debug(f"AC or AN is NA for record: {record.contig}:{record.pos}")
        return False
    
    # Calculate statistical power
    stat_power = binom.cdf(AC, AN, inhouse_common_cutoff)
    
    # Return True if variant should be considered common
    return stat_power > conf_level


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
        



def process_target_records( qrecord, 
                            crecord,
                            added_filter="INHOUSE_COMMON",
                            **kwargs ):
    qrecord = qrecord.record
    crecord = crecord.record

    is_common = determine_common_per_pysam_record(crecord, **kwargs)
    if is_common and "SDrecall" in qrecord.filter:
        qrecord.filter.add(added_filter)

    return VariantRecordWrapper(qrecord)



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
                                process_func,
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
            # Hit the bottom of the cohort iterator but still captured all the records at this position to buffer_precs and current_loc_precs
            break
        else:
            if compare_variant_positions(var2.record, next_var2.record) == 0:
                current_loc_2recs.append(next_var2)
            elif compare_variant_positions(var2.record, next_var2.record) == 1:
                buffer_var2.append(next_var2)
                break
            else:
                raise ValueError(f'The cohort records are not sorted by position: \n{buffer_var2}\n')
            
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
                merged_record = process_func(v1, v2, **kwargs)
                merged_records.add(merged_record)
                found_match = True
                break
        if not found_match:
            # If no matching cohort record is found, add the query record as is
            non_overlap_var1.add(v1)
    non_overlap_var2.update(set([v2 for v2 in current_loc_2recs if v2 not in matching_var2]))

    return merged_records, non_overlap_var1, non_overlap_var2, buffer_var1, buffer_var2



def imap_process_region(args):
    return process_region(*args)


@log_command
def process_region(region, 
                   query_vcf, 
                   cohort_vcf, 
                   added_filter = "INHOUSE_COMMON", 
                   qv_tag = None,
                   cv_tag = None,
                   inhouse_common_cutoff=0.05,
                   conf_level=0.999,
                   logger=logger):
    '''
    The function is not only used to return the modified overlapping records
    but also to return the records that are not overlapping with the cohort VCF file
    '''
    
    query_vcf = pysam.VariantFile(query_vcf)
    cohort_vcf = pysam.VariantFile(cohort_vcf)
    
    # Fetch the records in the specified region from both VCF files
    query_iterator = query_vcf.fetch(region, reopen=True)
    cohort_iterator = cohort_vcf.fetch(region, reopen=True)

    # Process edge cases when eighter of the iterator returns no records
    if check_empty_iterator(query_iterator):
        logger.warning(f"No records found in the query VCF for region {region}. All records from the cohort VCF will be returned as non-overlapping.")
        return (frozenset([]), frozenset([]), frozenset([ VariantRecordWrapper(r).pickable() for r in cohort_iterator ]))
    if check_empty_iterator(cohort_iterator):
        logger.warning(f"No records found in the cohort VCF for region {region}. All records from the query VCF will be returned as non-overlapping.")
        query_iterator = query_vcf.fetch(region)
        return (frozenset([]), frozenset([ VariantRecordWrapper(r).pickable() for r in query_iterator ]), frozenset([]))
    
    # Now I ensure both iterators are not empty
    if added_filter:
        query_vcf.header.add_meta('FILTER', items=[('ID', added_filter), ('Description', f'Variant is likely to be {added_filter}')])
    if qv_tag:
        query_vcf.header.add_meta('FILTER', items=[('ID', qv_tag), ('Description', f'Variant is from {qv_tag}')])
    if cv_tag:
        query_vcf.header.add_meta('FILTER', items=[('ID', cv_tag), ('Description', f'Variant is from {cv_tag}')])
        
    query_iterator = query_vcf.fetch(region)
    cohort_iterator = cohort_vcf.fetch(region)
    
    # Initialize variables
    merged_records = set([])
    buffer_crecs = []
    buffer_qrecs = []
    non_overlap_qrecs = set([])
    non_overlap_crecs = set([])

    cohort_record = None
    
    # Process records on the fly
    while True:
        try:
            if len(buffer_qrecs) > 0:
                query_record = buffer_qrecs.pop(0)
            else:
                query_record = VariantRecordWrapper(next(query_iterator))
        except StopIteration:
            # If we hit the bottom of the query iterator, we need to clear the buffer
            if cohort_record is not None:
                non_overlap_crecs.add(cohort_record)
            
            if len(buffer_crecs) > 0:
                non_overlap_crecs.update(set(buffer_crecs))
                
            while True:
                try:
                    cohort_record = VariantRecordWrapper(next(cohort_iterator))
                    logger.debug(f"Getting the next priority variant record from iterator: {cohort_record}")
                except StopIteration:
                    # If we hit the bottom of the cohort iterator, we need to clear the buffer
                    break
                else:
                    non_overlap_crecs.add(cohort_record)
                    
            break
        else:
            if cohort_record is None:
                pass
            elif compare_variant_positions(query_record.record, cohort_record.record) == 2:
                non_overlap_crecs.add(cohort_record)
            elif compare_variant_positions(query_record.record, cohort_record.record) == 1:
                non_overlap_qrecs.add(query_record)
                continue
            elif compare_variant_positions(query_record.record, cohort_record.record) == -1:
                raise ValueError('The chromosome of the query record and the cohort record are not the same: {} vs {}'.format(query_record, cohort_record))
            else:
                # The two records are at the same location
                merged_records, non_overlap_qrecs, non_overlap_crecs, buffer_qrecs, buffer_crecs = deal_with_same_loc_variants( query_record,
                                                                                                                                cohort_record,
                                                                                                                                query_iterator,
                                                                                                                                cohort_iterator,
                                                                                                                                buffer_qrecs,
                                                                                                                                buffer_crecs,
                                                                                                                                non_overlap_qrecs,
                                                                                                                                non_overlap_crecs,
                                                                                                                                merged_records,
                                                                                                                                process_target_records,
                                                                                                                                logger = logger,
                                                                                                                                added_filter=added_filter,
                                                                                                                                inhouse_common_cutoff=inhouse_common_cutoff,
                                                                                                                                conf_level=conf_level)
                cohort_record = None
                continue
            while True:
                try:
                    if len(buffer_crecs) > 0:
                        cohort_record = buffer_crecs.pop(0)
                    else:
                        cohort_record = VariantRecordWrapper(next(cohort_iterator))
                except StopIteration:
                    # If we hit the bottom of the cohort iterator, we need to clear the buffer
                    non_overlap_qrecs.add(query_record)
                    break
                else:
                    # Find the matching cohort record
                    if compare_variant_positions(query_record.record, cohort_record.record) == 2:
                        # Meaning the cohort_record is upstream of the query_record, then we need to continue the inner while_loop
                        non_overlap_crecs.add(cohort_record)
                        continue
                    elif compare_variant_positions(query_record.record, cohort_record.record) == 1:
                        # Meaning the cohort_record is downstream of the query_record, then we need to break the inner while_loop and let the outer while_loop to move to the next variant record
                        non_overlap_qrecs.add(query_record)
                        break
                    elif compare_variant_positions(query_record.record, cohort_record.record) == -1:
                        # This shouldn't happen because we are fetching records from the same chromosome region
                        raise ValueError('The chromosome of the query record and the cohort record are not the same: {} vs {}'.format(query_record, cohort_record))
                    else:
                        # The two records are at the same location
                        merged_records, non_overlap_qrecs, non_overlap_crecs, buffer_qrecs, buffer_crecs = deal_with_same_loc_variants( query_record, 
                                                                                                                                        cohort_record, 
                                                                                                                                        query_iterator, 
                                                                                                                                        cohort_iterator, 
                                                                                                                                        buffer_qrecs, 
                                                                                                                                        buffer_crecs, 
                                                                                                                                        non_overlap_qrecs, 
                                                                                                                                        non_overlap_crecs, 
                                                                                                                                        merged_records, 
                                                                                                                                        process_target_records,
                                                                                                                                        logger = logger,
                                                                                                                                        added_filter=added_filter,
                                                                                                                                        inhouse_common_cutoff=inhouse_common_cutoff,
                                                                                                                                        conf_level=conf_level)
                        cohort_record = None
                        break
    return (frozenset([ r.pickable() for r in merged_records ]), 
            frozenset([ r.pickable() for r in non_overlap_qrecs ]), 
            frozenset([ r.pickable() for r in non_overlap_crecs ]))


def annotate_inhouse_common(query_vcf = "", 
                            cohort_vcf = "", 
                            output_vcf = "", 
                            ref_genome = "",
                            added_filter = "INHOUSE_COMMON",
                            qv_tag = None,
                            cv_tag = None,
                            inhouse_common_cutoff = 0.05,
                            conf_level = 0.999,
                            threads = 4):
    # Sort the query VCF file
    sorted_query_vcf = sort_vcf(query_vcf, ref_genome, logger = logger)

    # Sort the cohort VCF file
    sorted_cohort_vcf = sort_vcf(cohort_vcf, ref_genome, logger = logger)

    # Open the sorted query VCF file
    bcf_query = pysam.VariantFile(sorted_query_vcf)
    bcf_cohort = pysam.VariantFile(sorted_cohort_vcf)

    # Open the output VCF file for writing
    tmp_output_vcf = prepare_tmp_file(suffix=".vcf.gz").name
    # Add new FILTER tags to the header
    if added_filter:
        bcf_query.header.add_meta('FILTER', items=[('ID', added_filter), ('Description', f'Variant is likely to be {added_filter}')])
    if qv_tag:
        bcf_query.header.add_meta('FILTER', items=[('ID', qv_tag), ('Description', f'Variant is from {qv_tag}')])
    if cv_tag:
        bcf_cohort.header.add_meta('FILTER', items=[('ID', cv_tag), ('Description', f'Variant is from {cv_tag}')])

    # Get the list of regions to process
    regions = list(dict.fromkeys(list(bcf_query.header.contigs.keys()) + list(bcf_cohort.header.contigs.keys())))
    
    logger.info(f"The query VCF file has records across {len(regions)} chromosomes, they are: {regions}")

    # Create a multiprocessing pool
    with mp.Pool(threads) as pool:
        # Process each region in parallel
        result_iterator = pool.imap_unordered(imap_process_region, [(region, 
                                                                    sorted_query_vcf, 
                                                                    sorted_cohort_vcf, 
                                                                    added_filter,
                                                                    qv_tag,
                                                                    cv_tag,
                                                                    inhouse_common_cutoff,
                                                                    conf_level) for region in regions if re.match(r"^chr[0-9MTXY]+$", region)])

        with pysam.VariantFile(tmp_output_vcf, 'w', header=bcf_query.header) as bcf_output:
            i=0
            for success, result, log_contents in result_iterator:
                i+=1
                print(f"\n************************************{i}_subprocess_start_for_compare_regional_variants************************************\n", file=sys.stderr)
                if success:
                    merged_recs, non_overlap_qrecs, non_overlap_crecs = result
                    for rec_tup in merged_recs:
                        chrom, pos, ref, alts, vid, qual, filters, info, samples = rec_tup
                        if cv_tag:
                            filters = [f for f in filters if f != cv_tag] + [cv_tag]
                        if qv_tag:
                            filters = [f for f in filters if f != qv_tag] + [qv_tag]
                        rec = bcf_output.new_record(contig=chrom, 
                                                    start=pos-1, # the new_record treats start as 0-based coordinates 
                                                    alleles=(ref,) + alts, 
                                                    id=vid, 
                                                    qual=qual, 
                                                    filter=filters, 
                                                    info=dict(info) )
                        for sample, values in samples:
                            for key, value in values:
                                # logger.info(f"Setting sample {sample} field {key} to value {value}")
                                try:
                                    rec.samples[sample][key] = value
                                except TypeError as te:
                                    # logger.warning(f"TypeError encountered when setting sample {sample} field {key} to value {value}: {te}.")
                                    continue

                        bcf_output.write(rec)
                    for qrec_tup in non_overlap_qrecs:
                        chrom, pos, ref, alts, vid, qual, filters, info, samples = qrec_tup
                        if qv_tag:
                            filters = [f for f in filters if f != qv_tag] + [qv_tag]
                        qrec = bcf_output.new_record(contig=chrom, 
                                                     start=pos-1, 
                                                     alleles=(ref,) + alts, 
                                                     id=vid, 
                                                     qual=qual, 
                                                     filter=filters, 
                                                     info=dict(info))
                        for sample, values in samples:
                            for key, value in values:
                                # logger.info(f"Setting sample {sample} field {key} to value {value}")
                                try:
                                    qrec.samples[sample][key] = value
                                except TypeError as te:
                                    # logger.warning(f"TypeError encountered when setting sample {sample} field {key} to value {value}: {te}.")
                                    continue
                        bcf_output.write(qrec)
                    
                    print(f"Successfully processed the variant records in one chromosome. The log info are:\n{log_contents}\n", file=sys.stderr)
                else:
                    error_mes, tb_str = result
                    raise ValueError(f"An error occurred: {error_mes}\nTraceback: {tb_str}\nThe error message is :\n{log_contents}\n")
                print(f"\n************************************{i}_subprocess_end_for_compare_regional_variants************************************\n", file=sys.stderr)

    # Close the VCF files
    bcf_query.close()
    bcf_cohort.close()

    # Sort and index the output VCF file
    sort_vcf(tmp_output_vcf, 
             ref_genome, 
             output_vcf = output_vcf)

    return output_vcf
    
    
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Compare and filter VCF files')
    parser.add_argument("-qv", '--query_vcf', help='Path to the query VCF file')
    parser.add_argument("-rv", '--cohort_vcf', help='Path to the cohort VCF file')
    parser.add_argument("-ov", '--output_vcf', help='Path to the output VCF file')
    parser.add_argument("-rg", '--ref_genome', help='Path to the reference genome file corresponding to the VCF file')
    parser.add_argument("-fl", "--filter_tag", type=str, help="The filter tag you want to added to the matched query VCF records", required=False, default="INHOUSE_COMMON")
    parser.add_argument("-qt", "--query_vcf_tag", type=str, help="The tag you want to added to the matched query VCF records", required=False, default=None)
    parser.add_argument("-ct", "--cohort_vcf_tag", type=str, help="The tag you want to added to the matched cohort VCF records", required=False, default=None)
    parser.add_argument("-ic", "--inhouse_common_cutoff", type=float, help="The cutoff for the inhouse common variants", required=False, default=0.01)
    parser.add_argument("-cl", "--confidence_level", type=float, help="The confidence level for the inhouse common variants", required=False, default=0.999)
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for parallel processing', required=False)

    args = parser.parse_args()

    annotate_inhouse_common(query_vcf = args.query_vcf, 
                            cohort_vcf = args.cohort_vcf, 
                            output_vcf = args.output_vcf, 
                            ref_genome = args.ref_genome,
                            added_filter = args.filter_tag,
                            qv_tag = args.query_vcf_tag,
                            cv_tag = args.cohort_vcf_tag,
                            threads = args.threads)