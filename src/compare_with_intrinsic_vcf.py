import vcfpy
import argparse as ap
import subprocess
import os
import re
import logging
from collections import OrderedDict

def main_process(query_vcf_path, intrin_vcf_path, lowerbound, output_vcf):
    intrin_vcf = vcfpy.Reader.from_path(intrin_vcf_path)
    intrin_recs = [r for r in intrin_vcf]
    query_vcf = vcfpy.Reader.from_path(query_vcf_path)
    query_recs = [r for r in query_vcf]
    query_head = query_vcf.header
    
    # Add filter label in header
    mapping = OrderedDict({'ID':"LIKELY_INTRINSIC",'Description':"The variant called is likely due to intrinsic difference between homologous sequences in the ref genome."})
    query_head.add_filter_line(mapping)
    mapping = OrderedDict({'ID':"UNLIKELY_INTRINSIC",'Description':"The variant called is unlikely due to intrinsic difference between homologous sequences in the ref genome. ALT/REF ratio is higher than theoretical value for wildtype."})
    query_head.add_filter_line(mapping)
    
    overlapped_recs = [ (r,rec) for r in query_recs for rec in intrin_recs if if_same_vcfrec(r, rec) ]
    non_overlap_recs = [ r for r in query_recs if r not in [ t[0] for t in overlapped_recs ] ]
    
    filtered_recs = [ modify_intrin_vars(t, lowerbound) for t in overlapped_recs ]
    final_recs = [ r for r in filtered_recs if r != None ] + non_overlap_recs
    
    if len(output_vcf) == 0:
        output_vcf = ".filtered.vcf".join(query_vcf_path.rsplit(".vcf", 1)).replace(".gz", "")
    else:
        ori_opv = output_vcf
        output_vcf = output_vcf.replace(".gz", "")

    cmd = 'if [[ -f {pv} ]]; then rm -f {pv}; fi'.format(pv = output_vcf)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    
    vw = vcfpy.Writer.from_path(path = output_vcf, header = query_head)
    [ vw.write_record(line) for line in final_recs ]
    vw.close()
    
    if re.search(r'\.gz$', ori_opv):
        cmd = 'bgzip -f {pv} && bcftools sort -Oz -o {gz} {gz} && tabix -f -p vcf {gz}'.format(pv = output_vcf, gz = re.sub(r'\.vcf$', '.vcf.gz', output_vcf))
    elif re.search(r'\.vcf$', ori_opv):
        cmd = 'bgzip -f {pv} && bcftools sort -Ov -o {pv} {gz}'.format(pv = output_vcf, gz = re.sub(r'\.vcf$', '.vcf.gz', output_vcf))
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    print(result.stdout)
    print(result.stderr)
    
    
def if_same_vcfrec(rec1, rec2):
    '''
    Input rec1 and rec2 are vcfpy.Record object
    '''
    if rec1.CHROM == rec2.CHROM and rec1.POS == rec2.POS and rec1.REF == rec2.REF and rec1.ALT[0] == rec2.ALT[0]:
        return True
    else:
        return False
    

def modify_intrin_vars(rec_tuple, lowerbound=1.1):
    '''
    Input rec is vcfpy.Record object
    '''
    rec = rec_tuple[0]
    irec = rec_tuple[1]
    
    rec.ID=["hom_intrinsic"]
    AD_values = [int(x) for x in rec.calls[0].data["AD"]]
    ref_dp = AD_values[0]
    alt_dp = AD_values[1]
    ra_ratio = ref_dp/alt_dp if alt_dp > 0 else 0
    if alt_dp == 0:
        logging.warning("This vcf record has 0 read carrying ALT allele? Take a look: {}".format(rec))
    
    iAD_values = [int(x) for x in irec.INFO['AD']]
    iref_dp = iAD_values[0] + 1
    ialt_dp = iAD_values[1]
    ira_ratio = iref_dp/ialt_dp if alt_dp > 0 else 0
    
    if ira_ratio == 0:
        rec.add_filter("UNLIKELY_INTRINSIC")
        return rec
    elif ra_ratio/ira_ratio <= lowerbound:
        rec.add_filter("LIKELY_INTRINSIC")
        return rec
    else:
        rec.add_filter("UNLIKELY_INTRINSIC")
        return rec
    
    

if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-qv", "--query_vcf", type=str, help="The path to the vcf file containing short var callings on homo regions", required=True)
    parser.add_argument("-iv", "--intrinsic_vcf", type=str, help="The path to the vcf file directly storing intrinsic vars", required=True)
    parser.add_argument("-lb", "--lowerbound", type=float, help="lowerbound of the REF/ALT depth ratio", required=False, default=1.1)
    parser.add_argument("-ov", "--output_vcf", type=str, help="output vcf path", required=True, default="")

    args=parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    main_process(query_vcf_path=args.query_vcf, intrin_vcf_path=args.intrinsic_vcf, lowerbound=args.lowerbound, output_vcf=args.output_vcf)