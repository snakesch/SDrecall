#!/usr/bin/env python3

# 3_maskedAlignPolyVar.py
# Description: Masked alignment and multi-ploidy variant calling.
# Author: Yang XT, She CH (2022)
# Contact: Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

import os
import pandas as pd
import glob
import argparse
import subprocess
from collections import Counter
import time
from src.utils import *
import logging

def maskedAlign(BAM_FILE, BED_DIR, MGENOME_DIR, FASTQ_DIR, REF_GENOME, NTHREADS):

    logging.info(f"Taking up {NTHREADS} threads.")
    sample_ID = os.path.basename(BAM_FILE).split(".")[0]

    # Available homologous regions coordinates
    os.chdir(os.path.join(BED_DIR, "homologous_regions"))
    bed_fns = (glob.glob("*_related_homo_regions.bed"))
    bed_regions = [fn.rstrip("_related_homo_regions.bed") for fn in bed_fns]

    # Available masked genomes
    os.chdir(MGENOME_DIR)
    fns = (glob.glob("*[0-9A-Z]_masked.fasta"))
    if any(os.path.getsize(f)/1024**3 > 2 for f in fns):
        logging.warning("Some masked genomes have size > 2Gb. BWA index IS algorithm may not be used.")
    mgenomes = ["_".join(mg[:-13].split("_")[1:]) for mg in fns]

    for key in Counter(mgenomes).keys():
        num = Counter(mgenomes)[key]
        if num > 1:
            raise KeyError(f"Multiple masked genomes found for region {key}")

    # Analytic regions
    regions = set(bed_regions) & set(mgenomes)
    logging.info(f"Realigning {len(regions)} regions ...")

    # Index reference genome
    if not os.path.exists(".".join(REF_GENOME.split(".")[:-1]) + ".dict"):
        cmd = f"gatk CreateSequenceDictionary -R {REF_GENOME}"
        executeCmd(cmd)
        cmd = f"samtools faidx {REF_GENOME}"
        executeCmd(cmd)

    # Masked alignment for each region
    high_ploidy = {}
    for region in regions:
        # Check if we have precisely 2 fastq files
        logging.debug(f"Processing {region}")
        os.chdir(FASTQ_DIR)
        fastq = (glob.glob(f"{sample_ID}_{region}-XA*.fq.gz"))
        if len(fastq) > 2:
            raise ValueError(f"More than 2 fastq files found for region {region}. ")
        elif len(fastq) == 0:
            logging.warning(f"No FASTQ files for region {region}. ")
            continue
        fastq = sorted([os.path.join(FASTQ_DIR, f) for f in fastq])

        # Check for multi-aligned reads first
        bedf = os.path.join(BED_DIR, "homologous_regions", f"{region}_related_homo_regions.bed")
        original_hq_depth = getDepth(BAM_FILE, bedf, NTHREADS, mode="HQ")
        logging.debug(f"Average depth before realignment (without multi-aligned reads): {original_hq_depth}")
        if original_hq_depth == 0.0:
            logging.info(f"No high-quality reads in {region}")
            continue

        # Calculate sequence depth
        original_depth = getDepth(BAM_FILE, bedf, NTHREADS)
        logging.debug(f"Average depth before realignment: {original_depth}")
        if original_depth == 0.0:
            logging.info(f"Zero coverage for region {region}")
            continue

        os.chdir(MGENOME_DIR)
        mg = (glob.glob(f"*{region}_masked.fasta"))[0]
        cmd = f"bwa index {mg}"
        executeCmd(cmd)

        output_fn = os.path.join(os.path.dirname(BAM_FILE), sample_ID + "_only_" + region + ".bam")
        cmd = fr"bwa mem -t {NTHREADS} -M -R '@RG\tID:{sample_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:{sample_ID}' {mg} {fastq[0]} {fastq[1]} | samtools sort -o bam -@ {NTHREADS} -o {output_fn}"
        executeCmd(cmd)

        XA_depth = getDepth(output_fn, bedf, NTHREADS)
        logging.debug(f"XA_depth = {XA_depth}")
        if XA_depth == 0.0:
            continue

        post_depth = XA_depth - original_hq_depth
        ploidy = 2*post_depth / original_depth
        logging.info(f"Read depth of {region} increased by {post_depth/original_hq_depth} fold.")
        logging.info(f"Estimated ploidy for region {region} : {ploidy}")

        if ploidy < 2.0:
            os.remove(output_fn)
        else:
            high_ploidy[region] = ploidy

    logging.info(f"Regions with ploidy >= 2 : {list(high_ploidy.keys())}")

    return high_ploidy

def multiploidVariantCaller(BAM_FILE, BED_DIR, MGENOME_DIR, REF_GENOME, NTHREADS, VCF_DIR, high_ploidy: dict, keep_vcf = False):

    # Create VCF directory for storing GATK outputs
    os.makedirs(VCF_DIR, exist_ok=True)

    sample_ID = os.path.basename(BAM_FILE).split(".")[0]

    processed = []
    for region, ploidy in high_ploidy.items():

        # Get path of masked genome
        os.chdir(MGENOME_DIR)
        mg = glob.glob(f"*{region}_masked.fasta")[0]
        # Get dict index if not found
        if not os.path.exists(mg[:-5] + 'dict'):
            cmd = f"gatk CreateSequenceDictionary -R {mg} -O {mg[:-5] + 'dict'} "
            executeCmd(cmd)

        # Get fai index if not found
        if not os.path.exists(mg + '.fai'):
            cmd = f"samtools faidx {mg} "
            executeCmd(cmd)

        # Get path of masked alignment BAM file
        masked_bamf = os.path.join(os.path.dirname(BAM_FILE), sample_ID + "_only_" + region + ".bam")

        # Index masked alignments
        if not os.path.exists(masked_bamf + ".bai"):
            cmd = f"samtools index {masked_bamf} "
            executeCmd(cmd)

        # Get BED path
        bedf = os.path.join(BED_DIR, "homologous_regions", f"{region}_related_homo_regions.bed")

        # Get principal component region
        pcbedf = os.path.join(BED_DIR, "principal_components", region + ".bed")

        # Output path
        HC_out = os.path.join(VCF_DIR, sample_ID + "_only_" + region + ".HC.g.vcf.gz")

        # HaplotypeCaller (Add --bamout for debugging)
        cmd = f"""gatk --java-options \"-Xmx8G\" HaplotypeCaller """ \
                    f"""-R {mg} """ \
                    f"""-I {masked_bamf} """ \
                    f"""-O {HC_out} """ \
                    f"""-L {pcbedf} """ \
                    f"""--ploidy {round(ploidy)} """ \
                    f"""--min-base-quality-score 15 """ \
                    f"""--force-active """ \
                    f"""--allow-non-unique-kmers-in-ref """ \
                    f"""--assembly-region-padding 150 """ \
                    f"""--annotate-with-num-discovered-alleles """ \
                    f"""--output-mode EMIT_ALL_CONFIDENT_SITES """ \
                    f"""--emit-ref-confidence GVCF """ \
                    f"""--min-pruning 1 """ \
                    f"""--recover-all-dangling-branches """ \
                    f"""--min-dangling-branch-length 2 """ \
                    f"""--allele-informative-reads-overlap-margin 5 """ \
                    f"""--pair-hmm-implementation EXACT """ \
                    f"""--max-reads-per-alignment-start 0 """ \
                    f"""--max-unpruned-variants 100 """ \
                    f"""--min-assembly-region-size 20 """ \
                    f"""--soft-clip-low-quality-ends """ \
                    f"""--native-pair-hmm-threads {NTHREADS} """ \
                    f"""--enable-dynamic-read-disqualification-for-genotyping """ \
                    f"""-DF WellformedReadFilter """ \
                    f"""-DF GoodCigarReadFilter """ \
                    f"""-DF NotSecondaryAlignmentReadFilter """ \
                    f"""-DF PassesVendorQualityCheckReadFilter """ \
                    f"""-DF MappedReadFilter """ \
                    f"""-ip 1000 """ \
                    f"""--debug-assembly """ \
                    f"""--kmer-size 41 """ \
                    f"""--max-num-haplotypes-in-population 1280 """ \
                    f"""--bam-writer-type CALLED_HAPLOTYPES """ \
                    f"""--G AS_StandardAnnotation """ \
                    f"""-L {pcbedf} """ \
                    f"""-imr OVERLAPPING_ONLY """ \
                    f"""-OBI true """ \
                    f"""--verbosity WARNING """
        executeCmd(cmd)
        vcf_header, subject_ID, records = loadVCF(HC_out)
        if records.shape[0] == 0:
            os.remove(HC_out)
            logging.warning(f"No variants remained after GATK HaplotypeCaller for region {region}")
            continue

        # GenotypeGVCFs
        gvcf_out = os.path.join(VCF_DIR, sample_ID + "_only_" + region + ".g.vcf.gz")
        cmd = f"""gatk --java-options \"-Xmx8G -XX:ConcGCThreads=2 -XX:ParallelGCThreads=5\" GenotypeGVCFs """ \
                    f"""-R {REF_GENOME} """ \
                    f"""-V {HC_out} """ \
                    f"""-O {gvcf_out} """ \
                    f"""--force-output-intervals {HC_out} """ \
                    f"""-G StandardAnnotation """ \
                    f"""-G AS_StandardAnnotation """
        executeCmd(cmd)
        vcf_header, subject_ID, records = loadVCF(gvcf_out)
        if records.shape[0] == 0:
            os.remove(HC_out)
            os.remove(gvcf_out)
            logging.warning(f"No variants remained after GATK GenotypeGVCFs for region {region}")
            continue

        _filter = records.apply(checkDP, axis=1, threshold=0)
        records = records[_filter]
        os.remove(gvcf_out)
        writeVCF(vcf_header, records, gvcf_out)

        # GenotypeGVCF gives malformed VCF files (Need to recompress)
        if os.path.exists(gvcf_out + '.tbi'):
            os.remove(gvcf_out + '.tbi')
        cmd = f"gunzip {gvcf_out} && bgzip {gvcf_out[:-3]} && tabix -f -p vcf {gvcf_out} "
        executeCmd(cmd)

        # LeftAlignAndTrimVariants
        vcf_out = os.path.join(VCF_DIR, sample_ID + "_only_" + region + ".vcf.gz")

        cmd = f"""gatk LeftAlignAndTrimVariants """ \
                    f"""-R {REF_GENOME} """ \
                    f"""-V {gvcf_out} """ \
                    f"""-O {vcf_out} """ \
                    f"""--max-leading-bases 2000 """ \
                    f"""--max-indel-length 2000 """ \
                    f"""--split-multi-allelics """ \
                    f"""--dont-trim-alleles """
        executeCmd(cmd)
        vcf_header, subject_ID, records = loadVCF(vcf_out)
        records = records[records["ALT"] != "*"]
        records = records.apply(convert_ploidy, axis=1)

        # Reconstruct header with contigs from reference genome
        other_header = [header for header in vcf_header if "##contig" not in header]
        fai_path = REF_GENOME + ".fai"
        fai = pd.read_csv(fai_path, sep="\t", header=None)
        fai[5] = "##contig=<ID=" + fai[0] + ",length=" + fai[1].astype(str) +">"
        contig_header = fai[5].tolist()
        header = other_header + contig_header
        logging.debug("Original VCF : \n{}".format(records.to_string()))
        os.remove(vcf_out)
        records = records.apply(convert_ploidy, axis=1)
        _filter = records.apply(checkDP, axis=1, threshold=5, how="first")
        records = records[_filter]

        dip_vcf_out = os.path.join(VCF_DIR, sample_ID + "_only_" + region + ".dip.vcf.gz")

        writeVCF(header, records, dip_vcf_out)
        logging.info(f"Written {region} VCF to {dip_vcf_out}")

        # Extract low coverage regions from VCF
        logging.debug(f"Processing BAM for {region} ...")
        cmd = "module load samtools; "
        cmd = cmd + f"samtools view -c -@ {NTHREADS} -L {bedf} "
        _cnt = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout

        if _cnt == 0:
            logging.warning(f"No reads in {region}.")
            continue

        # Descend to BAM directory
        os.chdir(os.path.dirname(BAM_FILE))
        mos_out_prefix = BAM_FILE.removesuffix("bam") + "highMQ"
        low_coverage_fn = BAM_FILE.removesuffix("bam") + region + ".bed"

        # Extract low coverage regions from BAM file
        logging.debug(f"Extracting low coverage regions for {region} ...")
        cmd = f"mosdepth -Q 10 -b {bedf} {mos_out_prefix} {BAM_FILE}"
        executeCmd(cmd)
        mos_bed = pd.read_csv(mos_out_prefix + ".per-base.bed.gz", compression="gzip", sep="\t", header=None)
        mos_bed = mos_bed[mos_bed[3] <= 15]
        mos_bed.to_csv(low_coverage_fn, sep="\t", index=False, header=False)

        # Remove files named *highMQ*
        mos_fns = (glob.glob("*highMQ"))
        for file in mos_fns:
            os.remove(file)

        # Extract low coverage regions from VCF
        final_vcf_out = os.path.join(VCF_DIR, sample_ID + "_only_" + region + "_extracted.vcf.gz")
        _final_vcf_out = final_vcf_out + ".tmp"

        dip_vcf_out_bgz = os.path.join(VCF_DIR, sample_ID + "_only_" + region + ".dip.vcf.bgz")
        cmd = f"gunzip -c {dip_vcf_out} | bgzip > {dip_vcf_out_bgz} "
        executeCmd(cmd)
        cmd = f"tabix -p vcf {dip_vcf_out_bgz} "
        executeCmd(cmd)

        cmd = f"bcftools view -R {low_coverage_fn} -Oz -o {_final_vcf_out} {dip_vcf_out_bgz}"
        executeCmd(cmd)
        cmd = f"bcftools sort -Oz -o {final_vcf_out} {_final_vcf_out}"
        executeCmd(cmd)
        os.remove(_final_vcf_out)
        os.remove(low_coverage_fn)
        os.remove(dip_vcf_out)

        cmd = f"tabix -p vcf {final_vcf_out} "
        executeCmd(cmd)

        if not keep_vcf:
            os.remove(HC_out)
            os.remove(HC_out + ".tbi")
            os.remove(gvcf_out)
            os.remove(gvcf_out + ".tbi")
            os.remove(dip_vcf_out_bgz)
            os.remove(dip_vcf_out_bgz + ".tbi")

        vcf_header, subject_ID, records = loadVCF(final_vcf_out)
        if records.shape[0] == 0:
            logging.error(f"No variants left for region {region}")
        else:
            logging.info(f"Multiploid variant calling completed for {region}")
            processed.append(region)
    logging.info(f"Processed regions : {processed}")

def getRegions(BED_DIR, VCF_DIR) -> str:

    """
    This function determines from VCF_DIR analytic regions to compare with intrinsic VCF.

    Analytic VCFs are concatenated and sorted.

    """
    # Get analytic regions for comparison
    bed_fns = (glob.glob(os.path.join(BED_DIR, "homologous_regions", "[A-Z0-9]*.bed")))
    vcf_fns = (glob.glob(os.path.join(VCF_DIR, "*_extracted.vcf.gz")))

    all_regions = {os.path.basename(f).split("_")[0] for f in bed_fns}
    logging.debug(f"All regions with BED coordinates : {all_regions}")

    vcf_regions = {os.path.basename(f)[:-17].split("_")[-1] for f in vcf_fns}
    logging.info(f"Regions with valid VCF : {vcf_regions}")
    sample_ID = os.path.basename(vcf_fns[0]).split("_")[0]

    if len(vcf_regions) == 0:
        raise ValueError("No valid regions extracted! ")
    vcf_fl = os.path.join(VCF_DIR, "region_filelist.tmp")

    with open(vcf_fl, "w") as f:
        for region in vcf_fns:
            f.write(region)
            f.write("\n")

    total_vcf = os.path.join(VCF_DIR, sample_ID + "_homo_regions.vcf.gz")
    cmd = f"bcftools concat -a --no-version -Oz -f {vcf_fl} -o {total_vcf} >/dev/null 2>&1 "
    executeCmd(cmd)

    vcf_header, subject_IDs, records = loadVCF(total_vcf)
    if records.shape[0] == 0:
        raise ValueError(f"No variants found in all VCFs! Check file list {vcf_fl}")

    query_vcf = os.path.join(VCF_DIR, sample_ID + "_extracted_homo_regions.vcf.gz")
    cmd = f"bcftools sort -Oz {total_vcf} -o {query_vcf} >/dev/null 2>&1 && tabix -f -p vcf {query_vcf} "
    executeCmd(cmd)

    os.remove(vcf_fl)
    os.remove(total_vcf)

    return query_vcf

def compareIntrinsic(query_vcf, intrinsic_vcf, lower_limit=1.1):
    """
    This function compares query VCF with intrinsic VCF to determine likely intrinsic and unlikely intrinsic variants.

    A lower limit can be specified to control variant filtering.

    """
    # Load query VCF and intrinsic VCF
    qvcf_header, qsubject_IDs, qrecords = loadVCF(query_vcf)
    ivcf_header, isubject_IDs, irecords = loadVCF(intrinsic_vcf)

    # Add filters to query VCF header
    likely_filter_head = "##FILTER=<ID=LIKELY_INTRINSIC,Description=\"The variant called is likely due to intrinsic difference between homologous sequences in the ref genome.\">"
    unlikely_filter_head = "##FILTER=<ID=UNLIKELY_INTRINSIC,Desciption=\"The variant called is unlikely due to intrinsic difference between homologous sequences in the ref genome. ALT/REF ratio is higher than theoretical value for wildtype.\">"

    for idx, header in enumerate(qvcf_header):
        if header.startswith("##FILTER"):
            qvcf_header.insert(idx, likely_filter_head)
            qvcf_header.insert(idx, unlikely_filter_head)
            break

    _nsubjs = qrecords.shape[1] - 9
    irecords["ALT_tmp"] = irecords["ALT"].str[0]
    qrecords["ALT_tmp"] = qrecords["ALT"].str[0]
    _all_ret = qrecords.merge(irecords, how="left", on=["#CHROM", "POS", "REF", "ALT_tmp"], suffixes=(None, "_i"))
    _all_ret = _all_ret.drop("ALT_tmp", axis=1)
    non_overlap_rec = _all_ret.loc[_all_ret.isnull().any(axis=1), :]
    overlap_rec = _all_ret.loc[~_all_ret.isnull().any(axis=1), :] # Likely intrinsic variants
    non_overlap_rec = non_overlap_rec.dropna(axis=1)
    overlap_rec.loc[:, "ID"] = "hom_intrinsic"

    # Intrinsic variant filter
    overlap_rec = overlap_rec.apply(intrinsicFilter, axis=1, args=(_nsubjs, lower_limit))
    final_rec = pd.concat([overlap_rec, non_overlap_rec], axis=0)
    final_vcf_out_tmp = query_vcf[:-6] + "filtered.tmp.vcf.gz"
    final_vcf_out = query_vcf[:-6] + "filtered.vcf.gz"
    writeVCF(qvcf_header, final_rec, final_vcf_out_tmp)

    # Output
    cmd = f"bcftools sort -Oz {final_vcf_out_tmp} -o {final_vcf_out} >/dev/null 2>&1  && tabix -p vcf {final_vcf_out} "
    executeCmd(cmd)
    os.remove(final_vcf_out_tmp)
    logging.info(f"Filtered VCF written to {final_vcf_out}")
    logging.info("****************** Intrinsic variant filtering completed ******************")

def intrinsicFilter(row, nsubjs: int, lower_limit: float):
    """

    A pandas apply method that computes qra_ratio and ira_ratio to determine intrinsic variants with a given lower limit.

    Input should contain AD values from both query VCF and intrinsic VCF.
    Output only returns query columns.

    """
    qinfo_lst = [row[8].split(":")] + [row[col].split(":") for col in range(9, nsubjs + 9, 1)]
    qtag = {info[0]: list(info[1:]) for info in zip(*qinfo_lst)}

    try:
        ad = list(map(int, qtag["AD"][0].split(","))) # Only query the first call
    except KeyError:
        return row[:nsubjs + 9]
    ref_dp, alt_dp = ad[0], ad[1]
    qra_ratio = ref_dp / alt_dp if alt_dp > 0 else 0
    if alt_dp == 0:
        logging.warning(f"Record {record} has 0 read with ALT allele.")

    iinfo_lst = [row["FORMAT_i"].split(":")] + [row[col].split(":") for col in range(row.index.get_loc("FORMAT_i") + 1, row.shape[0], 1)]
    itag = {info[0]: list(info[1:]) for info in zip(*iinfo_lst)}

    iad = list(map(int, itag["AD"][0].split(","))) # Only query the first call
    iref_dp = iad[0] + 1
    ialt_dp = iad[1]
    ira_ratio = iref_dp / ialt_dp if alt_dp > 0 else 0

    if ira_ratio == 0:
        row["FILTER"] = "UNLIKELY_INTRINSIC"
    elif qra_ratio / ira_ratio <= lower_limit:
        row["FILTER"] = "LIKELY_INTRINSIC"
    else:
        row["FILTER"] = "UNLIKELY_INTRINSIC"

    return row[:nsubjs + 9]

def cleanup(bamfp, outpath, vcfp="", refp=""):

    """
    This function takes multiple inputs of directories and clean up intermediate files generated in SDrecall. vcfp and refp are optional as out/vcf/ and out/ref/ by default.)
    """
    import shutil

    # Descend to BAM path
    os.chdir(os.path.dirname(bamfp))

    # Remove temporary masked alignment files
    for tmp_bamf in glob.glob("*_only_*"):
        os.remove(tmp_bamf)

    # Descend to out/
    os.chdir(outpath)

    # Remove multi-aligned reads
    if os.path.isdir("fastq"):
        shutil.rmtree("fastq")

    # Remove large masked genomes
    if os.path.isdir("masked_genome"):
        shutil.rmtree("masked_genome")

    # Remove VCFs if specified (Not recommended; currently not implemented)
    if vcfp and os.path.isdir(vcfp):

        # Descend to VCF directory
        os.chdir(vcfp)

        # Remove VCFs of individual regions
        for regionf in glob.glob("*_only_*"):
            os.remove(regionf)

        # Remove intrinic VCFs
        for ivcf in glob.glob("all_pc*"):
            os.remove(ivcf)

        # Remove raw extracted VCF
        for rawf in glob.glob("*extracted_homo_regions.vcf.gz*"):
            os.remove(rawf)

    # Remove reference files
    if refp and os.path.isdir(refp):

        # Descend to ref directory
        os.chdir(refp)

        # Remove homologous region BED files
        if os.path.isdir("homologous_regions"):
            shutil.rmtree("homologous_regions")

        # Remove principal component BED files
        if os.path.isdir("principal_components"):
            shutil.rmtree("principal_components")

if __name__ == "__main__":

    # Argparse setup
    parser = argparse.ArgumentParser(description = "Masked alignment and multi-ploidy variant calling.")
    parser._optionals.title = "Options"
    ROOT = os.path.dirname(__file__)
    parser.add_argument("--input_bam", type = str, help = "input BAM file", required = True)
    parser.add_argument("--bed_dir", type = str, help = "directory of BED files", default = None)
    parser.add_argument("--intrinsic_vcf", type = str, help = "path of intrinsic VCF", required = True)
                        # default = f"{ROOT}/ref/all_pc.realign.250.trim.vcf.gz")
    parser.add_argument("--lower", type = float, help = "lower bound for unlikely intrinsic variants", default = 1.1)
    parser.add_argument("--masked_genomes", type = str, help = "directory of masked genomes", default = None)
    parser.add_argument("--fastq_dir", type = str, help = "directory of FASTQ files", default = None)
    parser.add_argument("--ref_genome", type = str, help = "path of reference genome", default = None)
    parser.add_argument("--output_vcf", type = str, help = "directory to write VCF", default = None)
    parser.add_argument("--thread", type = int, help = "number of threads", default = 8)
    parser.add_argument("--keep_vcf", action = "store_true", help = "keep intermediate VCFs")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    logging.debug(f"Working in {ROOT}")

    BAM_FILE = args.input_bam
    BED_DIR = args.bed_dir or f"{ROOT}/ref/" # homologous_regions / principal_components
    MGENOME_DIR = args.masked_genomes or f"{ROOT}/out/masked_genome/"
    FASTQ_DIR = args.fastq_dir or f"{ROOT}/out/fastq/"
    REF_GENOME = args.ref_genome or f"{ROOT}/ref/ucsc.hg19.fasta"
    VCF_DIR = args.output_vcf or f"{ROOT}/out/vcf"
    NTHREADS = args.thread

    start = time.time()
    high_ploidy = maskedAlign(BAM_FILE, BED_DIR, MGENOME_DIR, FASTQ_DIR, REF_GENOME, NTHREADS)
    logging.info(f"****************** Masked alignment completed in {time.time() - start:.2f} seconds ******************")

    start = time.time()
    multiploidVariantCaller(BAM_FILE, BED_DIR, MGENOME_DIR, REF_GENOME, NTHREADS, VCF_DIR, high_ploidy, keep_vcf = args.keep_vcf)
    logging.info(f"****************** Multiploid variant calling completed in {time.time() - start:.2f} seconds ******************")

    start = time.time()
    if not os.path.exists(args.intrinsic_vcf):
        raise FileNotFoundError(f"Intrinsic VCF not found in path {args.intrinsic_vcf}")
    else:
        logging.info(f"Detected intrinsic VCF : {args.intrinsic_vcf}")

    query_vcf = getRegions(args.bed_dir, VCF_DIR)
    compareIntrinsic(query_vcf, args.intrinsic_vcf, lower_limit=args.lower)
    logging.info(f"****************** Comparison with intrinsic VCF completed in {time.time() - start:.2f} seconds ******************")
    if args.verbose.upper() != "DEBUG":
        logging.info("Cleaning up ... ")
        cleanup(BAM_FILE, os.path.dirname(VCF_DIR), refp=os.path.dirname(REF_GENOME))

