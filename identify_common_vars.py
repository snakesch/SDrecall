import pandas as pd
import numpy as np
from scipy.stats import binom
import pysam
import logging
import os
import re
from python_utils import na_value
import argparse


from src.log import logger
from src.utils import executeCmd, prepare_tmp_file


def determine_common_per_pysam_record(record, 
                                     AC_tag="AC", 
                                     AN_tag="AN", 
                                     inhouse_common_cutoff=0.05,
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


def identify_inhouse_common(query_vcf,
                           cohort_vcf,
                           output_vcf="",
                           inhouse_common_cutoff=0.05,
                           conf_level=0.999):
    """
    Identify common variants in a cohort VCF and annotate them in a query VCF
    
    Args:
        query_vcf: The VCF file to be annotated
        cohort_vcf: The cohort-level VCF with population data
        output_vcf: Output VCF file path (if empty, will modify input)
        inhouse_common_cutoff: Cutoff for common variant frequency
        conf_level: Confidence level threshold (default: 0.999)
        
    Returns:
        Path to the annotated VCF file
    """
    logger.info(f"Processing cohort VCF: {cohort_vcf}")
    logger.info(f"Query VCF to annotate: {query_vcf}")
    
    # Step 1: Create a VCF file of common variants from the cohort VCF
    common_variants_vcf_tmp = prepare_tmp_file(suffix='.vcf')
    common_variants_file = common_variants_vcf_tmp.name
    common_variants_vcf = f"{common_variants_file}.gz"
    
    # Open cohort VCF
    cohort_vcf_file = pysam.VariantFile(cohort_vcf)
    
    # Add inhouse_common filter to header
    header = cohort_vcf_file.header.copy()
    if "inhouse_common" not in header.filters:
        header.add_filter("inhouse_common", 
                            description="Variant common in inhouse cohort",
                            source="identify_inhouse_common")
    
    # Create output VCF for common variants
    with pysam.VariantFile(common_variants_file, 'w', header=header) as common_out:
        # Process each record
        common_count = 0
        total_count = 0
        
        for record in cohort_vcf_file:
            total_count += 1
            
            # Check if this variant is common based on our criteria
            is_common = determine_common_per_pysam_record(
                record,
                AC_tag="AC",
                AN_tag="AN",
                inhouse_common_cutoff=inhouse_common_cutoff,
                conf_level=conf_level
            )
            
            if is_common:
                common_count += 1
                # Add the inhouse_common filter
                record.filter.add("inhouse_common")
                # Write to common variants VCF
                common_out.write(record)
    
    cohort_vcf_file.close()
    logger.info(f"Found {common_count} common variants out of {total_count} total variants")
    
    # Compress and index the common variants VCF
    cmd = f"bgzip -f {common_variants_file} && tabix -f -p vcf {common_variants_vcf}"
    executeCmd(cmd, logger=logger)
    
    # Step 2: Set up output file
    if not output_vcf:
        final_output = prepare_tmp_file(suffix='.vcf.gz').name
    else:
        final_output = output_vcf
    
    # Step 3: Create a header file with the inhouse_common filter definition
    header_tmp = prepare_tmp_file(suffix='.txt')
    header_file = header_tmp.name
    with open(header_file, 'w') as hf:
        hf.write('##FILTER=<ID=inhouse_common,Description="Variant common in inhouse cohort">\n')
    
    # Step 4: Use bcftools to annotate the query VCF
    logger.info(f"Annotating query VCF with common variants")
    
    # First, index the query VCF if needed
    if not os.path.exists(f"{query_vcf}.tbi"):
        cmd = f"bcftools index --tbi {query_vcf}"
        executeCmd(cmd, logger=logger)
    
    # Now annotate the VCF with corrected syntax for VCF-based annotation
    # Using +inhouse_common to ensure it's added to the FILTER field
    cmd = (f"bcftools annotate \
             -a {common_variants_vcf} \
             --columns CHROM,POS,REF,ALT \
             --mark-sites FILTER:+inhouse_common \
             -h {header_file} \
             -Ou {query_vcf} | \
             bcftools sort -Oz -o {final_output}")
    
    logger.info(f"Running command: {cmd}")
    executeCmd(cmd, logger=logger)
    
    # Step 5: Index the output file
    cmd = f"bcftools index --tbi {final_output}"
    executeCmd(cmd, logger=logger)
    
    # Step 6: If modifying the input file, replace it
    if not output_vcf:
        cmd = f"mv {final_output} {query_vcf} && mv {final_output}.tbi {query_vcf}.tbi"
        executeCmd(cmd, logger=logger)
        logger.info(f"Inhouse common filters annotated to {query_vcf}")
        return query_vcf
    else:
        logger.info(f"Inhouse common filters annotated to {output_vcf}")
        return output_vcf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify common small variants in inhouse cohort")
    parser.add_argument("-q", "--query_vcf", required=True, help="Query VCF file to annotate")
    parser.add_argument("-c", "--cohort_vcf", required=True, help="Cohort-level VCF with population data")
    parser.add_argument("-o", "--output", default="", help="Output VCF file (default: modify input)")
    parser.add_argument("-t", "--cutoff", type=float, default=0.05, help="Frequency cutoff for common variants")
    parser.add_argument("-p", "--conf_level", type=float, default=0.999, help="Confidence level threshold")
    
    args = parser.parse_args()
    
    identify_inhouse_common(
        query_vcf=args.query_vcf,
        cohort_vcf=args.cohort_vcf,
        output_vcf=args.output,
        inhouse_common_cutoff=args.cutoff,
        conf_level=args.conf_level
    )