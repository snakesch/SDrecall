import pandas as pd
from Bio import bgzf
import gzip
import re
from io import StringIO
import subprocess
import logging
import os
import sys
import argparse as ap
from utils import read_vcf_gz, read_vcf_gz_chunkwise
            
def convert_vcf_to_diploidy_from_high_ploidy(vcf_df=pd.core.frame.DataFrame):
    vcf_body_dip = vcf_df.copy()
    
    for i in range(0, vcf_body_dip.shape[0]):
        '''
        # Deal with INFO fields
        info_fields = vcf_body_dip.iloc[i,7].split(";")
        info_dict = dict(zip([x.split("=")[0] for x in info_fields],[x.split("=")[-1] for x in info_fields]))
        if int(info_dict['AC']) >= 2:
            info_dict['AC'] = 2
        if int(info_dict['AN']) >= 2:
            info_dict['AN'] = 2
        info_dict['AF'] = float(info_dict['AC']/info_dict['AN'])
        # It seems we can leave AC and AN untouched
        '''
        # Deal with FORMAT fields
        format_fields = vcf_body_dip.iat[i, 8].split(":")
        vcf_headers = vcf_body_dip.columns.tolist()
        samples = [ vcf_headers[i] for i in range(9, len(vcf_headers)) ]
        for sample in samples:
            format_values = vcf_body_dip.iat[i,vcf_headers.index(sample)].split(":")
            format_dict = dict(zip(format_fields, format_values))

            # Dealing with GT first
            try:
                gt_val_list = [str(x) for x in str(format_dict['GT']).split(str(format_dict['GT'])[1])]
            except IndexError:
                logging.warning("It seems this sample {} in row {} does not have a valid GT value, lets take a look: \n{}".format(sample, i+1, vcf_body_dip.iloc[i,:].to_string()))
            else:
                if gt_val_list.count("1") >= 2:
                    format_dict['GT'] = "1/1"
                elif gt_val_list.count("1") == 1:
                    format_dict['GT'] = "0/1"
                elif gt_val_list.count("1") == 0 and "." in str(format_dict['GT']):
                    format_dict['GT'] = "."
                elif gt_val_list.count("1") == 0:
                    format_dict['GT'] = "0/0"

            # Dealing with PL & GQ if there is one key named PL
            if 'PL' in format_dict:
                pl_val_list = str(format_dict['PL']).split(",")
                if len(str(format_dict['GT'])) < 3:
                    format_dict['PL'] = "."
                    format_dict['GQ'] = "."
                elif gt_val_list.count("1") <= 1:
                    new_pl_vals = pl_val_list[:2] + [str(min([int(n) for n in pl_val_list[2:]]))]
                    format_dict['PL'] = ",".join([str(x) for x in new_pl_vals])
                    format_dict['GQ'] = sorted(set([int(x) for x in new_pl_vals]))[1]
                    if int(format_dict['GQ']) >= 99:
                        format_dict['GQ'] = 99
                elif gt_val_list.count("1") >= 2:
                    new_pl_vals = pl_val_list[:2] + ["0"]
                    format_dict['PL'] = ",".join([str(x) for x in new_pl_vals])
                    format_dict['GQ'] = sorted(set([int(x) for x in new_pl_vals]))[1]
                    if int(format_dict['GQ']) >= 99:
                        format_dict['GQ'] = 99
            
            # Check the results
            logging.debug(format_dict)
            vcf_body_dip.iloc[i,vcf_headers.index(sample)] = ":".join([str(v) for k,v in format_dict.items()])
    
    return vcf_body_dip
    
def convert_haploid_vcf_to_diploid(vcf_df=pd.core.frame.DataFrame):
    vcf_df_dip = vcf_df.copy()
    for i in range(0, vcf_df_dip.shape[0]): 
        # Deal with FORMAT fields
        format_fields = vcf_df_dip.iloc[i, 8].split(":")
        vcf_headers = vcf_df_dip.columns.tolist()
        samples = [ vcf_headers[i] for i in range(9, len(vcf_headers)) ]
        for sample in samples:
            format_values = vcf_df_dip.iloc[i,vcf_headers.index(sample)].split(":")
            format_dict = dict(zip(format_fields, format_values))

            # Dealing with GT first
            if str(format_dict['GT']) == "1":
                format_dict['GT'] = "1/1"
            elif str(format_dict['GT']) == "0":
                format_dict['GT'] = "0/0"
            elif str(format_dict['GT']) == ".":
                format_dict['GT'] = "./."

            # Dealing with PL & GQ if there is one key named PL
            if 'PL' in format_dict:
                pl_val_list = str(format_dict['PL']).split(",")
                print(pl_val_list)
                if str(format_dict['GT']) == "1/1" :
                    new_pl_vals = [ pl_val_list[0], pl_val_list[0], pl_val_list[-1] ]
                    print(new_pl_vals)
                    format_dict['PL'] = ",".join([str(x) for x in new_pl_vals])
                    format_dict['GQ'] = sorted(set([int(x) for x in new_pl_vals]))[1]
                    if int(format_dict['GQ']) >= 99:
                        format_dict['GQ'] = 99 
                if str(format_dict['GT']) == "0/0" or str(format_dict['GT']) == "./.":
                    new_pl_vals = [ pl_val_list[0], pl_val_list[-1], pl_val_list[-1] ]
                    format_dict['PL'] = ",".join([str(x) for x in new_pl_vals])
                    format_dict['GQ'] = sorted(set([int(x) for x in new_pl_vals]))[1]
                    if int(format_dict['GQ']) >= 99:
                        format_dict['GQ'] = 99
                        
    return vcf_df_dip
    

def main_convert(vcf_file="", chunksize=0, replace=False):    
    vcf_handle = read_vcf_gz(vcf_file)
    vcf_head = "".join(vcf_handle[0])
    vcf_body = vcf_handle[1]
    
    # Set the output path 
    output_vcf_file = ".dip.vcf".join(vcf_file.rsplit(".vcf", 1)).replace(".gz", "")
    logging.info("Writing diploid vcf to {}".format(output_vcf_file))
    
    # Output the vcf headers 
    with open(output_vcf_file, "w") as ov:
        ov.write(vcf_head)
    
    # Detect the vcf ploidy
    gt_list = str(vcf_body.iat[0, 9]).split(str(vcf_body.iat[0, 9])[1]) # choose the first row
    ploidy = len(gt_list)
    
    # Convert vcf ploidy
    if ploidy > 2:
        dip_vcf_body = convert_vcf_to_diploidy_from_high_ploidy(vcf_body)
    elif ploidy == 1:
        dip_vcf_body = convert_haploid_vcf_to_diploid(vcf_body)
    else:
        dip_vcf_body = vcf_body
    
    # Append the vcf body content to the output file
    dip_vcf_body.to_csv(output_vcf_file, index=False, sep='\t', mode='a')
    if replace:
        # mv the file to the original file
        code = subprocess.run("source miscellaneous.sh && isValidVCF {}".format(output_vcf_file), shell=True, capture_output=True).returncode
        if code:
            logging.warning("Invalid output VCF ({}). Original VCF kept.".format(output_vcf_file))
            cmd = 'bgzip -f ' + output_vcf_file + ' && bcftools sort -Oz -o {gz_vcf} {gz_vcf} && tabix -f -p vcf {gz_vcf} && ls -lh {gz_vcf}'.format(gz_vcf=output_vcf_file + ".gz")
        else:
            cmd = 'bgzip -f ' + output_vcf_file + ' && bcftools sort -Oz -o {gz_vcf} {gz_vcf} && tabix -f -p vcf '.format(gz_vcf=output_vcf_file + ".gz") + output_vcf_file + '.gz && mv ' + output_vcf_file + '.gz ' + vcf_file + \
                    ' && mv ' + output_vcf_file + '.gz.tbi ' + vcf_file + '.tbi && ls -lh ' + vcf_file
    else:
        # Bgzip the plain vcf file and create index of it.
        cmd = 'bgzip -f ' + output_vcf_file + ' && bcftools sort -Oz -o {gz_vcf} {gz_vcf} && tabix -f -p vcf {gz_vcf} && ls -lh {gz_vcf}'.format(gz_vcf=output_vcf_file + ".gz")
    logging.debug("Command stack: \n{}\n".format(cmd))
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    print(res.stdout, file=sys.stderr)
    

if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-vp", "--vcf_path", type=str, help="The path to the vcf files", required=True)
    parser.add_argument("-cs", "--chunk_size", type=int, help="The chunk size if vcf file is too big", required=False, default=0)
    parser.add_argument("-rp", "--replace", dest="replace", help="Whether to replace the original file.", action="store_true")
    
    args=parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    main_convert(vcf_file=args.vcf_path, chunksize=args.chunk_size, replace=args.replace)
    
    
    
        
    