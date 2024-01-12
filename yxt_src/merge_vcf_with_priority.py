import vcfpy
import argparse as ap
import numpy as np
import subprocess
import os
import io
import re
import logging
from collections import OrderedDict
from scipy.stats import binom_test
from postprocess_mergesvvcf import add_header_line
from python_utils import na_value



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def merge_vcf_heads(head_list, proband=None, ori_genome_dict="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.dict"):
    merged_lines = []
    samples = []
    for vcf_head in head_list:
        samples += vcf_head.samples.names
        for line in vcf_head.lines:
            try:
                if any(True if line.key == l.key and line.mapping['ID'] == l.mapping['ID'] else False for l in merged_lines):
                    continue
                else:
                    merged_lines.append(line)
            except AttributeError as e:
                # print(line)
                if line in merged_lines:
                    logger.info("Header line {} already included".format(line))
                    continue
                else:
                    merged_lines.append(line)
                    # print(e)
        
    samples = list(dict.fromkeys(samples))
    if proband: samples = [proband] + [s for s in samples if s != proband]
    sampleinfo = vcfpy.SamplesInfos(sample_names=samples)
    return vcfpy.Header(lines = merged_lines, samples = sampleinfo)



def prioritize_homozygous(rec, gq_cutoff=30):
    rec_call = rec.calls[0]
    ref_depth = int(rec_call.data["AD"][0])
    alt_depth = int(rec_call.data["AD"][1])
    if int(rec_call.data["GQ"]) <= gq_cutoff and alt_depth >= ref_depth:
        logger.warning(f"This variant has ambiguity between heterozygous call and homozygous call. So we just pick homozygous form for the sake of sensitivity: {rec}")
        return "1/1"
    else:
        return rec_call.data["GT"]


def main_process(priority_vcf_path, 
                 ori_vcf_path,
                 pv_tag = "GATK",
                 ov_tag = "SDrecall",
                 output_vcf="",
                 homo_gq_cutoff=None):
    
    ori_vfile = io.StringIO(subprocess.run("bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp -Ov {}".format(ori_vcf_path), shell=True, stdout=subprocess.PIPE, encoding="utf-8").stdout)
    priority_vfile = io.StringIO(subprocess.run("bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp -Ov {}".format(priority_vcf_path), shell=True, stdout=subprocess.PIPE, encoding="utf-8").stdout)
    
    ori_vcf = vcfpy.Reader.from_stream(ori_vfile)
    ori_head = ori_vcf.header
    priority_vcf = vcfpy.Reader.from_stream(priority_vfile)
    priority_head = priority_vcf.header


    merged_header = merge_vcf_heads([ori_head, priority_head])
    merged_header = add_header_line(merged_header, "FILTER", ov_tag, description="variants called from {}".format(ov_tag))
    merged_header = add_header_line(merged_header, "FILTER", pv_tag, description="variants called from {}".format(pv_tag))

    ori_recs = [ VCF_Rec(r) for r in ori_vcf ]
    all_ori_chrs = set([ r.__chrom__ for r in ori_recs ])

    priority_recs = [ VCF_Rec(r) for r in priority_vcf ]
    all_priority_chrs = set([ r.__chrom__ for r in priority_recs ])

    union_chrs = all_ori_chrs.union(all_priority_chrs)
    
    orec_chr_map = { chrom: [r for r in ori_recs if r.__chrom__ == chrom] for chrom in union_chrs }
    prec_chr_map = { chrom: [r for r in priority_recs if r.__chrom__ == chrom] for chrom in union_chrs }
    non_overlap_precs = set([])
    non_overlap_orecs = []
    match_precs = set([])
    buffer_precs = { chrom: [] for chrom in union_chrs }

    for chrom, chr_orecs in orec_chr_map.items():
        if len(chr_orecs) == 0:
            non_overlap_precs = non_overlap_precs.union(set(prec_chr_map[chrom]))
            continue
        
        chr_precs = prec_chr_map[chrom]
        if len(chr_precs) == 0:
            non_overlap_orecs = non_overlap_orecs + chr_orecs
            continue
        
        qindex = 0
        for iindex, chr_orec in enumerate(chr_orecs):
            continue_orec = False
            non_overlap_orecs.append(chr_orec)
            found_match = False

            # Inspect the buffer precs
            while True:
                if len(buffer_precs[chrom]) == 0:
                    continue_orec = False
                    break
                # The buffer is not empty, the first prec in the buffer is downstream than the last irec
                chr_prec = buffer_precs[chrom][0]
                if chr_prec > chr_orec:
                    # The first prec in the buffer is still downstream than the newest irec
                    continue_orec = True
                    break
                elif chr_prec < chr_orec:
                    # The first prec in the buffer is upstream than the newest irec
                    chr_prec = buffer_precs[chrom].pop(0)
                    non_overlap_precs.add(chr_prec)
                else:
                    chr_prec = buffer_precs[chrom].pop(0)
                    if chr_prec == chr_orec:
                        logger.debug("Found a match:\n")
                        logger.debug("The iterating query rec: {}".format(chr_prec))
                        logger.debug('The iterating itrinsic rec: {}\n'.format(chr_orec))
                        match_prec = chr_prec
                        match_prec.__orirecord__.FILTER = [f for f in list(dict.fromkeys(match_prec.__orirecord__.FILTER + chr_orec.__orirecord__.FILTER))]
                        orec_call = chr_orec.__orirecord__.calls[0]
                        prec_call = chr_prec.__orirecord__.calls[0]
                        
                        orec_call.data["GT"] = prioritize_homozygous(chr_orec.__orirecord__)
                        prec_call.data["GT"] = prioritize_homozygous(chr_prec.__orirecord__)

                        try:
                            match_prec.__orirecord__.calls[0].data["GT"] = prec_call.data["GT"] if prec_call.data["GT"].count("1") >= orec_call.data["GT"].count("1") else orec_call.data["GT"]
                        except TypeError:
                            logger.error("One of the records containing NoneType GT info. Check them out:\nReference record:{},{}\nQuery record:{},{}\n".format(orec_call.data, orec_call.gt_type, prec_call.data, prec_call.gt_type))
                            raise TypeError
                        else:
                            logger.debug("We found a match at {}:{}:{}->{}, ref rec call is {} and query rec call is {}.".format(chr_orec.__chrom__, chr_orec.__start__, chr_orec.__ref__, chr_orec.__alt__, orec_call.data["GT"].count("1"), prec_call.data["GT"].count("1")))
                            logger.debug("Now the new merged record's call is {}".format(match_prec.__orirecord__.calls[0].data["GT"]))
                        
                        match_precs.add(match_prec)
                        found_match = True
                    else:
                        non_overlap_precs.add(chr_prec)
                    continue

                if len(buffer_precs[chrom]) == 0:
                    continue_orec = False
                    break

            if continue_orec and iindex < (len(chr_orecs) - 1):
                if found_match:
                    non_overlap_orecs.pop(-1)
                continue
            
            # Iterate through the non-buffer precs
            while qindex < len(chr_precs):
                chr_prec = chr_precs[qindex]
                qindex += 1
                if chr_prec < chr_orec:
                    non_overlap_precs.add(chr_prec)
                elif chr_prec > chr_orec:
                    if iindex < (len(chr_orecs) - 1):
                        buffer_precs[chrom].append(chr_prec)
                        break
                    else:
                        non_overlap_precs.add(chr_prec)
                else:
                    if chr_prec == chr_orec:
                        logger.debug("Found a match:\n")
                        logger.debug("The iterating query rec: {}".format(chr_prec))
                        logger.debug('The iterating reference rec: {}\n'.format(chr_orec))
                        match_prec = chr_prec
                        match_prec.__orirecord__.FILTER = [f for f in list(dict.fromkeys(match_prec.__orirecord__.FILTER + chr_orec.__orirecord__.FILTER))]
                        orec_call = chr_orec.__orirecord__.calls[0]
                        prec_call = chr_prec.__orirecord__.calls[0]

                        orec_call.data["GT"] = prioritize_homozygous(chr_orec.__orirecord__)
                        prec_call.data["GT"] = prioritize_homozygous(chr_prec.__orirecord__)
                        
                        try:
                            match_prec.__orirecord__.calls[0].data["GT"] = prec_call.data["GT"] if prec_call.data["GT"].count("1") >= orec_call.data["GT"].count("1") else orec_call.data["GT"]
                        except TypeError:
                            logger.error("One of the records containing NoneType GT info. Check them out:\nReference record:{},{}\nQuery record:{},{}\n".format(orec_call.data, orec_call.gt_type, prec_call.data, prec_call.gt_type))
                            raise TypeError
                        else:
                            logger.debug("We found a match at {}:{}:{}->{}, ref rec call is {} and query rec call is {}.".format(chr_orec.__chrom__, chr_orec.__start__, chr_orec.__ref__, chr_orec.__alt__, orec_call.data["GT"].count("1"), prec_call.data["GT"].count("1")))
                            logger.debug("Now the new merged record's call is {}".format(match_prec.__orirecord__.calls[0].data["GT"]))
                        
                        match_precs.add(match_prec)
                        found_match = True
                        if iindex < (len(chr_orecs) - 1):
                            break
                    else:
                        non_overlap_precs.add(chr_prec)
            
            if found_match:
                non_overlap_orecs.pop(-1)
        
            non_overlap_precs = non_overlap_precs.union(set(buffer_precs[chrom]))
            

    non_overlap_orecs = set(non_overlap_orecs)

    for prec in match_precs:
        prec.__orirecord__.FILTER = [f for f in list(dict.fromkeys(prec.__orirecord__.FILTER)) if f not in [pv_tag, ov_tag]] + [pv_tag, ov_tag]
    for prec in non_overlap_precs:
        prec.__orirecord__.FILTER = [f for f in list(dict.fromkeys(prec.__orirecord__.FILTER)) if f not in [pv_tag]] + [pv_tag]
    for orec in non_overlap_orecs:
        orec.__orirecord__.FILTER = [f for f in list(dict.fromkeys(orec.__orirecord__.FILTER)) if f not in [ov_tag]] + [ov_tag]
    

    final_recs = match_precs.union(non_overlap_precs).union(non_overlap_orecs)
    final_recs = list(dict.fromkeys(final_recs))

    if len(output_vcf) == 0:
        output_vcf = ".merged.vcf".join(priority_vcf_path.rsplit(".vcf", 1)).replace(".gz", "")

    cmd = 'if [[ -f {pv} ]]; then rm -f {pv}; fi'.format(pv = output_vcf)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    
    with vcfpy.Writer.from_path(path = output_vcf, header = merged_header) as vw:
        rec_hashes = set([])
        for rec in final_recs:
            not_ref = rec.genotype != 0
            not_mis = rec.genotype != "miss"
            not_na = not na_value(rec.__alt__)
            rec.__orirecord__.REF = rec.__orirecord__.REF.upper()
            for alt in rec.__orirecord__.ALT:
                alt.value = alt.value.upper()
            rec.__orirecord__.ID = [e for e in rec.__orirecord__.ID if len(e) > 0]
            if all([len(l) == 1 for l in rec.__orirecord__.ID]):
                rec.__orirecord__.ID = [ "".join(rec.__orirecord__.ID) ]
            if not_ref and not_mis and not_na:
                if not rec.__hash__ in rec_hashes:
                    if homo_gq_cutoff:
                        ref_depth = int(rec.__orirecord__.calls[0].data["AD"][0])
                        alt_depth = int(rec.__orirecord__.calls[0].data["AD"][1])
                        if int(rec.__orirecord__.calls[0].data.get("GQ", 99)) <= homo_gq_cutoff and alt_depth >= ref_depth:
                            logger.warning(f"This variant has ambiguity between heterozygous call and homozygous call. So we just pick homozygous form for the sake of sensitivity: {rec}")
                            rec.__orirecord__.calls[0].data["GT"] = "1/1"
                    vw.write_record(rec.__orirecord__)
                    rec_hashes.add(rec.__hash__)
                else:
                    logger.warning("This record has been written before, so skip the current one:\n{}".format(rec.__orirecord__))
            
    
    if re.search(r'\.gz$', output_vcf):
        cmd = '/home/yangyxt/anaconda3/envs/SDrecall/bin/bcftools norm -a --atom-overlaps . -m-both --multi-overlaps 0 -Ou {opv} | bcftools norm -d exact -Ou - | bcftools filter -e \'ALT == "*"\' -Ou - | bcftools sort -Oz -o {tmp_opv} - && mv {tmp_opv} {opv} && tabix -f -p vcf {opv}'.format(opv = output_vcf, tmp_opv = re.sub(r'\.vcf\.gz$', '.tmp.vcf.gz', output_vcf))
    elif re.search(r'\.vcf$', output_vcf):
        cmd = '/home/yangyxt/anaconda3/envs/SDrecall/bin/bcftools norm -a --atom-overlaps . -m-both --multi-overlaps 0 -Ou {opv} | bcftools norm -d exact -Ou - | bcftools filter -e \'ALT == "*"\' -Ou - | bcftools sort -Ov -o {tmp_opv} - && mv {tmp_opv} {opv}'.format(opv = output_vcf, tmp_opv = re.sub(r'\.vcf\.gz$', '.tmp.vcf.gz', output_vcf))
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    print(result.stdout)
    print(result.stderr)
    


class VCF_Rec(object):
    '''
    This class is used as a node object in Graph 
    '''
    def __init__(self, vcfpy_record, caller=".", rough = True, onlysite=True):
        self.__chrom__ = str(vcfpy_record.CHROM)
        self.__start__ = int(vcfpy_record.POS)
        try:
            self.__end__ = int(vcfpy_record.INFO['END'])
        except:
            self.__end__ = int(vcfpy_record.POS)
        self.__orirecord__ = vcfpy_record
        self.infos = vcfpy_record.INFO
        self.__ref__ = str(vcfpy_record.REF).upper()
        try:
            self.depth = float(vcfpy_record.calls[0].data["DP"])
        except IndexError:
            self.depth = 0
        except KeyError:
            self.depth = np.nan
        except TypeError:
            self.depth = np.nan
        try:
            self.genotype = vcfpy_record.calls[0].gt_type
            if self.depth == 0:
                self.genotype = "miss"
            if self.depth in [np.nan]:
                self.genotype = "miss"
        except IndexError:
            self.genotype = "miss"
        else:
            if "GT" in vcfpy_record.calls[0].data:
                proband_gt = vcfpy_record.calls[0].data["GT"]
                self.genotype = proband_gt.count("1") if proband_gt is not None else None
            elif "AD" in vcfpy_record.calls[0].data:
                ad_values = vcfpy_record.calls[0].data["AD"]
                if None in ad_values: 
                    logger.warning(f"This record does not have a valid GT field and a None value is in its AD values: {vcfpy_record}")
                    ad_values = [ 0 if ad is None else ad for ad in ad_values ]
                if len(ad_values) < 2:
                    raise ValueError(f"This record does not have a valid GT field and it only has one AD values: {vcfpy_record}")                   
                if ad_values[0] > 0 and ad_values[1] == 0:
                    self.genotype = 0
                else:
                    self.genotype = 1
            if self.genotype is None:
                self.genotype = "miss"
        try:
            self.__alt__ = str(vcfpy_record.ALT[0].value).upper()
        except IndexError as ie:
            logger.debug("This vcfpy record seems being a ref record. Let's take a look at it: \n{}".format(vcfpy_record))
            self.__alt__ = self.__ref__
            self.genotype = 0
        self.__caller__ = str(caller)
        self.__samples__ = list(vcfpy_record.call_for_sample.keys())
        self.filters = vcfpy_record.FILTER
        self.ids = vcfpy_record.ID
        self.rough = rough
        self.onlysite = onlysite

    def __hash__(self):
        return hash(self.__chrom__) ^ hash(self.__start__) ^ hash(self.__ref__) ^ hash(self.__alt__) ^ hash(self.genotype)

    def __eq__(self, other):
        '''
        Test if this object is equal to other object
        '''
        if self.__chrom__ != other.__chrom__:
            return False
        elif self.onlysite:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__:
                return True
            else:
                return False
        elif self.rough:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__ and self.genotype != 0 and other.genotype != 0 and self.genotype != "miss" and other.genotype != "miss":
                return True
            else:
                return False
        else:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__ and self.genotype == other.genotype and self.genotype != 0 and other.genotype != 0 and self.genotype != "miss" and other.genotype != "miss":
                return True
            else:
                return False
    
    def __lt__(self, other):
        if self.__chrom__ == other.__chrom__:
            if self.__start__ < other.__start__:
                return True
            elif self.__start__ == other.__start__:
                if self.__end__ < other.__end__:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
        
    def __gt__(self, other):
        if self.__chrom__ == other.__chrom__:
            if self.__start__ < other.__start__:
                return False
            elif self.__start__ == other.__start__:
                if self.__end__ > other.__end__:
                    return True
                else:
                    return False
            else:
                return True
        else:
            return True

    def __str__(self):
        '''
        Convert this object to string
        '''
        return "chr:{},start:{},end:{},ref:{},alt:{},caller:{},genotype:{},samples{}".format(self.__chrom__, self.__start__, self.__end__, self.__ref__, self.__alt__, self.__caller__, self.genotype, self.__samples__)

    def __repr__(self):
        '''
        Representation for printing
        '''
        return "chr:{},start:{},end:{},ref:{},alt:{},caller:{},genotype:{},samples{}".format(self.__chrom__, self.__start__, self.__end__, self.__ref__, self.__alt__, self.__caller__, self.genotype, self.__samples__)
    


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-pv", "--priority_vcf", type=str, help="The path to the vcf file containing short var callings on homo regions", required=True)
    parser.add_argument("-ov", "--original_vcf", type=str, help="The path to the vcf file directly storing reference vars", required=True)
    parser.add_argument("-pt", "--pri_vcf_tag", type=str, help="The filter tag you want to added to the matched query VCF records", required=False)
    parser.add_argument("-ot", "--ori_vcf_tag", type=str, help="The filter tag you want to added to the matched query VCF records", required=False, default = "SDrecall")
    parser.add_argument("-hc", "--homozygous_gq_cutoff", type=int, help="The GQ cutoff to force homozygous call from heterozygous variants", required=False, default = None)
    parser.add_argument("-mv", "--merged_vcf", type=str, help="output vcf path", required=True, default="")

    args=parser.parse_args()
    main_process(priority_vcf_path=args.priority_vcf, 
                 ori_vcf_path=args.original_vcf, 
                 pv_tag=args.pri_vcf_tag,
                 ov_tag=args.ori_vcf_tag,
                 output_vcf=args.merged_vcf,
                 homo_gq_cutoff=args.homozygous_gq_cutoff)
