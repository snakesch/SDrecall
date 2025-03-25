import pysam
import numpy as np
import pandas as pd
import argparse as ap
import logging
import uuid
import os
import subprocess
import tempfile
from subprocess import PIPE


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



def stat_ad_to_dict(bam_file, ref_genome, region = None, empty_dict={}, logger = logger):
    # Build up an AD query dict by bcftools mpileup
    
    try:
        if region is None:
            bam_ad_file = f"{bam_file}.ad"
            # If no region is specified, use the original approach
            cmd = f"""bcftools mpileup -Ou --fasta-ref {ref_genome} -a FORMAT/AD --indels-2.0 -q 0 -Q 15 {bam_file} | \
                  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n' - > {bam_ad_file}"""
        else:
            bam_ad_file = f"{bam_file}.{region.replace(':', '_').replace('-', '_')}.ad"
            # Use sambamba slice to extract the region first, then pipe to bcftools
            cmd = f"""sambamba slice -q {bam_file} {region} | \
                  bcftools mpileup -Ou --fasta-ref {ref_genome} -a FORMAT/AD --indels-2.0 -q 0 -Q 15 - | \
                  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n' - > {bam_ad_file}"""
        
        executeCmd(cmd, logger = logger)
        
        # Use pandas read_csv with optimized parameters
        ad_table = pd.read_csv(
            bam_ad_file, 
            sep="\t", 
            names=["chrom", "pos", "ref", "alt", "ad"],
            dtype={"chrom": str, "pos": int, "ref": str, "alt": str, "ad": str},
            na_values=["", "<*>"]
        ).dropna(subset=["alt"])
        
        ad_table.loc[:, "alt"] = ad_table["alt"].str.rstrip(",<*>").str.upper()  # Convert all the bases to uppercase
        ad_table.loc[:, "ref"] = ad_table["ref"].str.upper()  # Convert the reference allele to uppercase

        ad_expanded = ad_table["ad"].str.split(",", expand=True).replace({None: np.nan, "": np.nan, "0": np.nan}).astype(float).dropna(axis=1, how="all")
        alt_expanded = ad_table["alt"].str.split(",", expand=True).replace({None: np.nan, "": np.nan}).dropna(axis=1, how="all")

        if alt_expanded.shape[1] < 1 or ad_expanded.shape[1] <= 1:
            logger.warning(f"No ALT allele found in this BAM file: {bam_file} at region {region}. Skip this entire script")
            logger.warning(f"Look at the ad table: \n{ad_table.to_string(index=False)}\n")
            return None

        logger.info("For bam {}, the original AD table looks like \n{}\nThe AD expanded table looks like: \n{}\nThe ALT expanded table looks like: \n{}\n".format(bam_file,
                                                                                                                ad_table[:10].to_string(index=False),
                                                                                                                ad_expanded[:10].to_string(index=False),
                                                                                                                alt_expanded[:10].to_string(index=False)))

        # Initialize the nested dictionary
        nested_ad_dict = {chrom: {} for chrom in ad_table["chrom"].unique()}

        # Iterate over the rows of dfA and dfB
        for i in range(len(ad_table)):
            chrom = ad_table.iloc[i, 0]
            pos = ad_table.iloc[i, 1] # position
            ref = ad_table.iloc[i, 2] # reference allele
            alts = alt_expanded.iloc[i, :].dropna().unique().tolist() # alternative alleles
            allele_depths = ad_expanded.iloc[i, 1:].dropna().astype(int).tolist() # allele depths

            if pd.isna(allele_depths).all():
                continue

            if len(alts) == 0:
                logger.warning(f"The ALT alleles are null for the row {i} of the ad table: \n{ad_table.iloc[i, :].to_string(index=False)}\n")
                continue

            total_dp = sum(allele_depths)

            # Initialize the inner dictionary if the outer key is not present
            if pos not in nested_ad_dict[chrom]:
                nested_ad_dict[chrom][pos] = empty_dict
            
            assert len(alts) == len(allele_depths), f"The number of alternative alleles and the number of allele depths are not the same at {chrom}:{pos} in bam file {bam_file}"
            for alt, allele_depth in zip(alts, allele_depths):
                # Now we need to handle insertions
                if len(alt) > len(ref) and len(ref) > 1:
                    assert ref[0] == alt[0], f"The reference allele and the alternative allele have different first base: {ref} and {alt} at {chrom}:{pos} in bam file {bam_file}"
                    # This is an insertion, we need to crop out the reference allele to have normalized format
                    ref_pos = alt.find(ref)
                    ins_pos = ref_pos + len(ref)
                    alt = alt[ins_pos-1:]
                    ref = ref[-1]
                elif len(ref) > len(alt) and len(alt) > 1:
                    assert ref[0] == alt[0], f"The reference allele and the alternative allele have different first base: {ref} and {alt} at {chrom}:{pos} in bam file {bam_file}"
                    # This is a deletion, we need to crop out the alternative allele to have normalized format
                    del_seq = ref[1:]
                    alt_crop = alt[1:] # We need to crop out the sequence from the del_seq in the right most position
                    crop_pos = del_seq.rfind(alt_crop)
                    del_seq = del_seq[:crop_pos]
                    ref = ref[0] + del_seq
                    alt = alt[0]

                # Add the pair of inner key-value
                nested_ad_dict[chrom][pos][alt] = allele_depth

            nested_ad_dict[chrom][pos]["DP"] = total_dp
        
        return nested_ad_dict
    
    finally:
        # Clean up the temporary file
        if os.path.exists(bam_ad_file):
            os.remove(bam_ad_file)



def extract_allele_depths(variant, bam_300x, bam_30x, realigned_bam, reference_genome):
    """
    Extract allele depths for a variant from two BAM files and a realigned BAM.
    
    :param variant: pysam VariantRecord object
    :param bam_300x: Path to the 300x coverage BAM file
    :param bam_30x: Path to the 30x coverage BAM file
    :param reference_genome: Path to the reference genome file
    :return: Dictionary containing allele depths for both BAM files and realigned BAM
    """
    # Consider the edge case where realigned BAM is not found
    # Just skip the AD extraction in realigned BAM and fill NA to the corresponding fields

    result = {'300x': {'DP':np.nan}, '30x': {'DP':np.nan}, 'realigned': {'DP':np.nan}}
    logger.info(f"Extracting allele depths for variant {variant.chrom}:{variant.start}-{variant.stop} among {bam_300x}, {bam_30x}, and {realigned_bam}")
    region_str = f"{variant.chrom}:{variant.start - 10}-{variant.stop + 10}"

    ad_300x_dict = stat_ad_to_dict(bam_300x, reference_genome, region = region_str)
    ad_30x_dict = stat_ad_to_dict(bam_30x, reference_genome, region = region_str)
    ad_realigned_dict = stat_ad_to_dict(realigned_bam, reference_genome, region = region_str)

    if ad_300x_dict is not None:
        dp = ad_300x_dict[variant.chrom].get(variant.pos, {'DP': np.nan})['DP']
        result['300x']['DP'] = dp
        for alt in variant.alts:
            ad = ad_300x_dict[variant.chrom].get(variant.pos, {}).get(alt, 0)
            if dp > 0:
                result['300x'][alt] = ad
            else:
                result['300x'][alt] = np.nan

    if ad_30x_dict is not None:
        dp = ad_30x_dict[variant.chrom].get(variant.pos, {'DP': np.nan})['DP']
        result['30x']['DP'] = dp
        for alt in variant.alts:
            ad = ad_30x_dict[variant.chrom].get(variant.pos, {}).get(alt, 0)
            if dp > 0:
                result['30x'][alt] = ad
            else:
                result['30x'][alt] = np.nan

    if ad_realigned_dict is not None:
        dp = ad_realigned_dict[variant.chrom].get(variant.pos, {'DP': np.nan})['DP']
        result['realigned']['DP'] = dp
        for alt in variant.alts:
            ad = ad_realigned_dict[variant.chrom].get(variant.pos, {}).get(alt, 0)
            if dp > 0:
                result['realigned'][alt] = ad
            else:
                result['realigned'][alt] = np.nan
            
    return result




def main_function(golden_vcf, test_vcf, bam_300x, bam_30x, realigned_bam, reference_genome):
    tmp_tag = str(uuid.uuid4())
    gv_uniq_dir = os.path.join(os.path.dirname(golden_vcf), os.path.basename(golden_vcf).replace(".vcf.gz", f".uniq_vars.{tmp_tag}"))
    # Extract variants unique to the golden VCF
    unique_variants_file = extract_missing_tps(golden_vcf, test_vcf, gv_uniq_dir)
    logger.info(f"The unique variants in {golden_vcf} is stored in {unique_variants_file}\n")

    if not os.path.exists(bam_300x):
        raise FileNotFoundError(f"The 300x BAM file {bam_300x} does not exist")
    
    if not os.path.exists(bam_30x):
        raise FileNotFoundError(f"The 30x BAM file {bam_30x} does not exist")
    
    if not os.path.exists(realigned_bam):
        raise FileNotFoundError(f"The realigned BAM file {realigned_bam} does not exist")

    if not os.path.exists(reference_genome):
        raise FileNotFoundError(f"The reference genome {reference_genome} does not exist")
    
    # Load the unique variants into a pysam VCF file
    unique_variants = pysam.VariantFile(unique_variants_file)
    
    # Prepare the output DataFrame
    output_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', '300x_DP', '300x_ALT', '30x_DP', '30x_ALT', 'REALIGNED_DP', 'REALIGNED_ALT'])
    
    # Iterate over the unique variants and extract allele depths
    idx = 0
    for variant in unique_variants:
        allele_depths = extract_allele_depths(variant, bam_300x, bam_30x, realigned_bam, reference_genome)
        
        for alt in variant.alts:
            output_df.loc[idx] = [variant.chrom, variant.pos, variant.ref, alt, 
                                  allele_depths['300x']['DP'], allele_depths['300x'].get(alt, 0) if allele_depths['300x']['DP'] > 0 else np.nan, 
                                  allele_depths['30x']['DP'], allele_depths['30x'].get(alt, 0) if allele_depths['30x']['DP'] > 0 else np.nan, 
                                  allele_depths['realigned']['DP'], allele_depths['realigned'].get(alt, 0) if allele_depths['realigned']['DP'] > 0 else np.nan]
            idx += 1
    
    # Use test_vcf to offer a column of GT
    output_meta_table = test_vcf.replace(".vcf.gz", ".missing_TP_meta.tsv")
    output_df.to_csv(output_meta_table, sep="\t", index=False)
    logger.info(f"The final meta table is stored at {output_meta_table}. And it looks like:\n{output_df[:10].to_string(index=False)}\n")

    # Now we can use the table to identify the variants that are missing because lost of alternative alleles in the down-sample process
    true_allele_absent = np.logical_and(output_df['REALIGNED_ALT'] <= 1, output_df['300x_ALT'] <= 1) | \
                         np.logical_and(output_df['REALIGNED_ALT'].isna(), output_df['300x_ALT'].isna())
    
    true_allele_lost = np.logical_and(output_df['30x_ALT'] == 0, output_df['300x_ALT'] > 0) | \
                       np.logical_and(output_df['30x_ALT'].isna(), output_df['300x_ALT'] > 0) | \
                       np.logical_or(output_df['REALIGNED_ALT'] == 0, output_df['REALIGNED_ALT'].isna())

    serious_allele_bias = np.logical_and((output_df['30x_ALT']/(output_df['30x_DP']) <= 0.1), (output_df['REALIGNED_ALT']/(output_df['REALIGNED_DP']) <= 0.1))
    excl_tp_bools = true_allele_absent | true_allele_lost | serious_allele_bias
    lost_alt_variants = output_df[excl_tp_bools]

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
    # Allow 7 arguments: 
    # 1. The path to the golden_vcf file
    # 2. The path to the vcf file to be tested
    # 3. The path to the output directory of the vcf file uniq to golden_vcf (the output of bcftools isec)
    # 4. The path to the BAM file with 300x coverage
    # 5. The path to the BAM file with downsampled 30x coverage
    # 6. The path to the realigned BAM file
    # 7. The path to the reference genome
    # Thats all
    parser = ap.ArgumentParser(description="Identify allele bias in a VCF file")
    parser.add_argument("-gv", "--golden_vcf", type=str, help="Path to the golden VCF file")
    parser.add_argument("-tv", "--test_vcf", type=str, help="Path to the VCF file to be tested")
    parser.add_argument("--bam_300x", type=str, help="Path to the BAM file with 300x coverage")
    parser.add_argument("--bam_30x", type=str, help="Path to the BAM file with 30x coverage")
    parser.add_argument("--realigned_bam", type=str, help="Path to the realigned BAM file")
    parser.add_argument("-rg", "--reference_genome", type=str, help="Path to the reference genome")
    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main_function(args.golden_vcf, args.test_vcf, args.bam_300x, args.bam_30x, args.realigned_bam, args.reference_genome)


