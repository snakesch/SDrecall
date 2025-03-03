import pysam
import numba
import numpy as np
import pandas as pd
import argparse as ap
import logging
import uuid
import os
import subprocess
import sys
import tempfile
import glob
import re
from subprocess import PIPE
bash_utils_hub = "/paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh"



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def prepare_tmp_file(tmp_dir="/paedyl01/disk1/yangyxt/test_tmp", **kwargs):
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass
    
    return tempfile.NamedTemporaryFile(dir = "/paedyl01/disk1/yangyxt/test_tmp", delete = False, **kwargs)




def executeCmd(cmd, stdout_only = False, shell="/home/yangyxt/miniforge3/envs/ngs_pipeline/bin/bash", logger=logger):
    logger.info("About to run this command in shell invoked within python: \n{}\n".format(cmd))

    if stdout_only:
        result = subprocess.run(cmd, shell=True, executable=shell, capture_output=True)
        
    else:
        result = subprocess.run(cmd, shell=True, executable=shell, stderr=subprocess.STDOUT, stdout=PIPE)

    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        if stdout_only:
            logger.error("Error in \n{}\nAnd the output goes like:\n{}\n".format(" ".join(cmd_lst), result.stderr.decode()))
        else:
            logger.error("Error in \n{}\nAnd the output goes like:\n{}\n".format(" ".join(cmd_lst), result.stdout.decode()))
        if cmd_lst[1][0] != "-":
            raise RuntimeError
        else:
            raise RuntimeError

    if stdout_only:
        logger.info(f"Ran the following shell command inside python:\n{cmd}\nAnd it receives a return code of {code}\nThe process looks like this:\n{result.stderr.decode()}\n\n")
    else:
        logger.info(f"Ran the following shell command inside python:\n{cmd}\nAnd it receives a return code of {code}, the output goes like this:\n{result.stdout.decode()}\n\n")
    return result.stdout.decode()



def extract_missing_tps(golden_vcf, test_vcf, output_dir):
    """
    Compare two VCF files using bcftools isec and return variants unique to the golden VCF.
    
    :param golden_vcf: Path to the golden VCF file
    :param test_vcf: Path to the test VCF file
    :param output_dir: Path to the output directory for bcftools isec results
    :return: Path to the VCF file containing variants unique to the golden VCF
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the bcftools isec command
    cmd = f"bcftools isec -p {output_dir} -n=1 -c all {golden_vcf} {test_vcf} && bcftools view -Oz -o {output_dir}/0000.vcf.gz {output_dir}/0000.vcf && bcftools index -t -f {output_dir}/0000.vcf.gz"
    
    # Execute the command using the existing executeCmd function
    executeCmd(cmd)
    
    # The file containing variants unique to the golden VCF will be named 0000.vcf in the output directory
    unique_variants_file = os.path.join(output_dir, "0000.vcf.gz")
    
    if not os.path.exists(unique_variants_file):
        logger.error(f"Expected output file {unique_variants_file} not found")
        raise FileNotFoundError(f"bcftools isec did not produce the expected output file: {unique_variants_file}")
    
    return unique_variants_file



def extract_allele_depths(variant, bam_300x, bam_30x, reference_genome):
    """
    Extract allele depths for a variant from two BAM files and a realigned BAM.
    
    :param variant: pysam VariantRecord object
    :param bam_300x: Path to the 300x coverage BAM file
    :param bam_30x: Path to the 30x coverage BAM file
    :param reference_genome: Path to the reference genome file
    :return: Dictionary containing allele depths for both BAM files and realigned BAM
    """
    realigned_bam = locate_chunk_raw_realigned_bam(reference_genome, bam_30x, variant)
    # Consider the edge case where realigned BAM is not found
    # Just skip the AD extraction in realigned BAM and fill NA to the corresponding fields

    result = {'300x': {'ref':np.nan, 'alt_1':np.nan}, '30x': {'ref':np.nan, 'alt_1':np.nan}, 'realigned': {'ref': 0, 'alt': 0, 'misaligned_alt': 0}}
    logger.info(f"Extracting allele depths for variant {variant.chrom}:{variant.start}-{variant.stop} among {bam_300x}, {bam_30x}, and {realigned_bam}")
    if reference_genome == "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta":
        nochr_genome = "/paedyl01/disk1/yangyxt/indexed_genome/GRCh37/human_g1k_v37.fasta"
    else:
        nochr_genome = reference_genome

    try:
        with pysam.AlignmentFile(bam_300x, "rb", reference_filename=nochr_genome) as bam_300x_file, \
             pysam.AlignmentFile(bam_30x, "rb", reference_filename=reference_genome) as bam_30x_file:
            
            for bam_file, coverage in [(bam_300x_file, '300x'), (bam_30x_file, '30x')]:
                if coverage == '300x' and nochr_genome != reference_genome:
                    chrom = variant.chrom.replace("chr", "")
                else:
                    chrom = variant.chrom

                pileup = bam_file.pileup(chrom, variant.pos - 1, variant.pos)
                for pileupcolumn in pileup:
                    if pileupcolumn.pos == variant.pos - 1:
                        bases = pileupcolumn.get_query_sequences(add_indels=True)
                        ref_count = bases.count(variant.ref)
                        alt_counts = [bases.count(alt) for alt in variant.alts]
                        
                        result[coverage]['ref'] = ref_count
                        for i, alt in enumerate(variant.alts):
                            result[coverage][f'alt_{i+1}'] = alt_counts[i]
            
        # Process realigned BAM
        if not realigned_bam:
            result['realigned']['ref'] = np.nan
            result['realigned']['alt'] = np.nan
            result['realigned']['misaligned_alt'] = np.nan
            return result

        with pysam.AlignmentFile(realigned_bam, "rb", reference_filename=reference_genome) as realigned_bam_file:
            for pileupcolumn in realigned_bam_file.pileup(variant.chrom, variant.pos - 1, variant.pos, truncate=True):
                if pileupcolumn.pos == variant.pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            hp_tag = pileupread.alignment.get_tag('HP') if pileupread.alignment.has_tag('HP') else None
                            
                            is_misaligned = hp_tag and ('HIGHVD' in hp_tag or 'LOWQUAL' in hp_tag)
                            
                            if base == variant.ref:
                                result['realigned']['ref'] += 1
                            elif base in variant.alts:
                                result['realigned']['alt'] += 1
                                if is_misaligned:
                                    result['realigned']['misaligned_alt'] += 1
    
    except Exception as e:
        logger.error(f"Error extracting allele depths for variant {variant.chrom}:{variant.pos}: {str(e)}")
        raise
    
    return result




def locate_chunk_raw_realigned_bam(ref_genome, bam_30x, pysam_variant):
    assembly = os.path.basename(ref_genome).split(".")[1]
    sample_ID = os.path.basename(bam_30x).split(".")[0]
    region_str = f"{pysam_variant.chrom}:{pysam_variant.start}-{pysam_variant.stop}"

    """
    Locate the chunk ID by region.
    
    :param region_str: String representing the genomic region
    :param sample_ID: Sample ID
    :param assembly: Genome assembly (default: "hg19")
    :return: Chunk ID if found, None otherwise
    """
    # Create a temporary BED file
    with tempfile.NamedTemporaryFile(mode='w+t', delete=False, suffix='.bed') as temp_bed:
        region_bed = temp_bed.name
        
        # Parse the region string
        region_parts = region_str.split(':')
        region_chr = region_parts[0]
        region_start, region_end = map(int, region_parts[1].split('-'))
        
        # Adjust the region
        region_start -= 20
        region_end += 20
        
        # Write to the temporary BED file
        temp_bed.write(f"{region_chr}\t{region_start}\t{region_end}\n")
    
    # Find chunk BED files
    chunk_beds_pattern = f"/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{assembly}/{sample_ID}.SD.deduped.pooled.raw.deduped.chunk*.bed"
    logger.info(f"The chunk bed pattern is {chunk_beds_pattern}")
    chunk_beds = glob.glob(chunk_beds_pattern)
    
    chunk_bed_overlap = None
    for chunk_bed in chunk_beds:
        if re.search(r'\.chunk(\d+)\.bed', os.path.basename(chunk_bed)) is None:
            continue
        # Use bedtools intersect
        cmd = f"bedtools intersect -a {region_bed} -b {chunk_bed} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        intersect_interval_count = int(result.stdout.strip())
        
        if intersect_interval_count > 0:
            logger.info(f"{chunk_bed} is overlapping with the specified region {region_str}")
            chunk_bed_overlap = chunk_bed
            break
    
    # Clean up the temporary file
    os.unlink(region_bed)
    
    if not chunk_bed_overlap:
        logger.warning(f"No chunk bed file is overlapping with the specified region {region_str}")
        return None
    
    # Extract chunk ID from the filename
    chunk_id = re.search(r'\.chunk(\d+)\.bed', os.path.basename(chunk_bed_overlap)).group(1)
    temp = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/{refgen_tag}/{sampID}.SD.deduped.pooled.raw.deduped.chunk{chunkID}.bam"
    return temp.format(refgen_tag=assembly, sampID=sample_ID, chunkID=chunk_id)



def main_function(golden_vcf, test_vcf, bam_300x, bam_30x, reference_genome):
    tmp_tag = str(uuid.uuid4())
    gv_uniq_dir = os.path.join(os.path.dirname(golden_vcf), os.path.basename(golden_vcf).replace(".vcf.gz", f".uniq_vars.{tmp_tag}"))
    # Extract variants unique to the golden VCF
    unique_variants_file = extract_missing_tps(golden_vcf, test_vcf, gv_uniq_dir)
    logger.info(f"The unique variants in {golden_vcf} is stored in {unique_variants_file}\n")
    
    # Load the unique variants into a pysam VCF file
    unique_variants = pysam.VariantFile(unique_variants_file)
    
    # Prepare the output DataFrame
    output_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', '300x_REF', '300x_ALT', '30x_REF', '30x_ALT', 'REALIGNED_REF', 'REALIGNED_ALT', 'MISALIGNED_ALT'])
    
    # Iterate over the unique variants and extract allele depths
    idx = 0
    for variant in unique_variants:
        allele_depths = extract_allele_depths(variant, bam_300x, bam_30x, reference_genome)
        
        output_df.loc[idx] = [variant.chrom, variant.pos, variant.ref, variant.alts[0], 
                              allele_depths['300x']['ref'], allele_depths['300x']['alt_1'], 
                              allele_depths['30x']['ref'], allele_depths['30x']['alt_1'], 
                              allele_depths['realigned']['ref'], allele_depths['realigned']['alt'], 
                              allele_depths['realigned']['misaligned_alt']]
        idx += 1
    
    # Use test_vcf to offer a column of GT
    output_meta_table = test_vcf.replace(".vcf.gz", ".missing_TP_meta.tsv")
    output_df.to_csv(output_meta_table, sep="\t", index=False)
    logger.info(f"The final meta table is stored at {output_meta_table}. And it looks like:\n{output_df[:10].to_string(index=False)}\n")

    # Now we can use the table to identify the variants that are missing because lost of alternative alleles in the down-sample process
    lost_alt_variant_bools = (output_df['30x_ALT'] == 0) | \
                             (output_df['REALIGNED_ALT'] == 0) | \
                             (output_df['30x_ALT'].isna())
    lost_alt_variants = output_df[lost_alt_variant_bools]

    # Now we need to remove those variants from the golden_vcf and output a new golden vcf without them
    lost_alt_variants_df = lost_alt_variants[['CHROM', 'POS', 'REF', 'ALT']]
    # Use pysam to iterate through the golden vcf. And write the record to a new golden vcf file if the variant is not contained in the lost_alt_variants_vcf
    new_golden_vcf = golden_vcf.replace(".vcf.gz", ".no_lost_alt.vcf.gz")
    with pysam.VariantFile(golden_vcf) as golden_vcf_file, pysam.VariantFile(new_golden_vcf, 'w', header=golden_vcf_file.header) as new_golden_vcf_file:
        for record in golden_vcf_file:
            if lost_alt_variants_df.loc[(lost_alt_variants_df['CHROM'] == record.chrom) & \
                                        (lost_alt_variants_df['POS'] == record.pos) & \
                                        (lost_alt_variants_df['REF'] == record.ref) & \
                                        (lost_alt_variants_df['ALT'] == record.alts[0]), :].empty:
                new_golden_vcf_file.write(record)

    logger.info(f"The new golden vcf without the variants lost alternative alleles in the down-sample process is stored at {new_golden_vcf}\n")
    cmd = "bcftools index -t " + new_golden_vcf
    executeCmd(cmd)

    return new_golden_vcf




if __name__ == "__main__":
    # Just use argparse to parse the command line arguments
    # Allow 6 arguments: 
    # 1. The path to the golden_vcf file
    # 2. The path to the vcf file to be tested
    # 3. The path to the output directory of the vcf file uniq to golden_vcf (the output of bcftools isec)
    # 4. The path to the BAM file with 300x coverage
    # 5. The path to the BAM file with downsampled 30x coverage
    # 6. The path to the reference genome
    # Thats all
    parser = ap.ArgumentParser(description="Identify allele bias in a VCF file")
    parser.add_argument("-gv", "--golden_vcf", type=str, help="Path to the golden VCF file")
    parser.add_argument("-tv", "--test_vcf", type=str, help="Path to the VCF file to be tested")
    parser.add_argument("--bam_300x", type=str, help="Path to the BAM file with 300x coverage")
    parser.add_argument("--bam_30x", type=str, help="Path to the BAM file with 30x coverage")
    parser.add_argument("-rg", "--reference_genome", type=str, help="Path to the reference genome")
    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main_function(args.golden_vcf, args.test_vcf, args.bam_300x, args.bam_30x, args.reference_genome)


