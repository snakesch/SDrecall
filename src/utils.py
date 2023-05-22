def timing(f):
    
    from functools import wraps
    from time import time
    import logging
    
    logger = logging.getLogger("root")
    
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        logger.info('Process :%r (args:[%r, %r]) completed in %2.2f min' % \
          (f.__name__, args, kw, (te-ts)/60))
        return result
    return wrap

def getIntersect(all_blocks: list):
    '''
    This function takes in a list of lists and groups sub-lists if their intersects is not NULL.
    '''
    # Base case
    if len(all_blocks) == 0:
        return None
    
    intersect = []
    tmp = []
    base = set()
    for cur in range(len(all_blocks)):
        # Does it intersect?
        if base.intersection(all_blocks[cur]) == set() and cur != 0:
            intersect.append(tuple(tmp))
            tmp.clear()
            tmp.append(all_blocks[cur])
            base = set(all_blocks[cur])
        else:
            tmp.append(all_blocks[cur])
            # Find new members
            for _idx in all_blocks[cur]:
                if _idx not in base:
                    base.add(_idx)
    intersect.append(tuple(tmp))
    
    return intersect

def loadVCF(path, omit_record=False):
    """
    Input a VCF file and returns:
    
    1. header (lines begin with #)
    2. subject ID list
    3. record dataframe
    
    """
    from gzip import GzipFile
    import pandas as pd
    
    header = []
    subjects = []
    with GzipFile(path) as f:
        line = next(f)
        while line:
            if line.decode()[1] == "#":
                header.append(line.decode().strip())
            elif line.decode()[0] == "#":
                col_names = line.decode().strip().split("\t")
                subjects = col_names[9:]
            else:
                break
            try:
                line = next(f)
            except StopIteration:
                break
    if omit_record:
        return header, subjects
    
    vcf = pd.read_csv(path, sep="\t", na_filter=False, engine="c", comment="#", header=None, compression="gzip")
    vcf.columns = col_names
    if header == None or subjects == None:
        raise RuntimeError("Incorrect VCF format.")
    return header, subjects, vcf

def writeVCF(header, records, outpath) -> None:
    """Writes to VCF with header, records and output path"""
    import gzip
    import pandas
    
    with gzip.open(outpath, "wb") as f:
        f.write("\n".join(header).encode())
        f.write("\n".encode())
        f.write("\t".join(records.columns.tolist()).encode())
        f.write("\n".encode())
  
    records.to_csv(outpath, mode="a", compression="infer", index=False, sep="\t", header=False)
    
    return

def executeCmd(cmd) -> None:
    
    import subprocess
    
    code = subprocess.run(cmd, shell=True, capture_output=False).returncode
    cmd_lst = cmd.split(" ")   
    if code != 0:
        if cmd_lst[1][0] != "-":
            raise RuntimeError("Error in " + " ".join(cmd_lst[0:1]))
        else:
            raise RuntimeError("Error in {}".format(cmd_lst[0]))
    return

def getDepth(bamf, target_bed, nthreads, mode="normal"):
    """
    Utility to compute average read depth for a given alignment file
    
    Recommend to use nthreads = 8 to enhance speed.
    
    """
    import subprocess, os
    if not os.path.exists(bamf):
        raise FileNotFoundError(f"File not found - {bamf}")
    elif not os.path.exists(target_bed):
        raise FileNotFoundError(f"File not found - {target_bed}")
    elif mode not in ["normal", "HQ"]:
        raise ValueError(f"Unknown mode - {mode}. Please select from normal / HQ.")
        
    if mode == "normal":
        cmd = f"""samtools depth -@ {nthreads} -a -b {target_bed} """\
               f"""-g DUP,UNMAP,QCFAIL,SECONDARY {bamf} | """\
               """awk '{sum+=$3} END{ if (NR != 0) { printf "%s",sum/NR;} else {print 0;}}'"""
    else:
        cmd = f"""samtools view -h -@ {nthreads} {bamf} | grep -v XA:Z | """\
               f"""samtools depth -@ {nthreads} -a -b {target_bed} -Q 30 """\
               """-g DUP,UNMAP,QCFAIL,SECONDARY - | """\
               """awk '{sum+=$3} END{ if (NR != 0) {printf "%s",sum/NR;} else {print 0;}}'"""
    p = subprocess.run([cmd], shell=True, capture_output = True)
    
    return float(p.stdout)

def convert_ploidy(row):
    """
    This function converts variants from higher diploidy / haploidy to diploidy by changing GT, PL and GQ fields.
    """

    info_lst = [row[8].split(":")] + [row[col].split(":") for col in range(9, row.shape[0], 1)]
    tag = {info[0]: list(info[1:]) for info in zip(*info_lst)}
    try:
        _len = len(tag["GT"][0])
    except KeyError:
        return row
    
    if _len > 3:
        # Multiploidy
        ol = tag["GT"]
        for idx, gt in enumerate(ol):
            if gt.count("1") >= 2:
                ol[idx] = "1/1"
            elif gt.count("1") == 1:
                ol[idx] = "0/1"
            elif gt.count("1") == 0 and "." in gt:
                ol[idx] = "."
            else:
                ol[idx] = "0/0"
        tag["GT"] = ol
        
        if "PL" in tag.keys():
            pl_lst, gq_lst = [], []
            for gt, pl in zip(tag["GT"], tag["PL"]):
                if len(gt) < 3:
                    gq_lst.append(".")
                    pl_lst.append(".")
                else:
                    old_pl = list(map(int, pl.split(",")))
                    if gt.count("1") <= 1:
                        new_pl = old_pl[0:2] + [min(old_pl[2:])]
                    elif gt.count("1") >= 2:
                        new_pl = old_pl[0:2] + [0]

                    pl_lst.append(",".join(map(str, new_pl)))
                    second_smallest = sorted(set(new_pl))[1]
                    if second_smallest > 99:
                        gq_lst.append("99")
                    else:
                        gq_lst.append(str(second_smallest))
            # Mutate original tag (GT is static in the second part)
            tag["PL"], tag["GQ"] = pl_lst, gq_lst
            
        # Mutate original VCF row
        row[8] = ":".join(tag.keys())     
        for i in range(9, row.shape[0], 1):
            info = list(map(list, zip(*tag.values())))[i-9]
            row[i] = ":".join(info)
    elif _len == 1:
        # haploid
        ol = tag["GT"]
        for idx, gt in enumerate(ol):
            if gt == "1":
                ol[idx] = "1/1"
            elif gt == "0":
                ol[idx] = "0/1"
            elif gt == ".":
                ol[idx] = "."
            else:
                raise ValueError(f"Unrecognized haploid genotype - {gt}")
        tag["GT"] = ol
        
        if "PL" in tag.keys():
            pl_lst, gq_lst = [], []
            for gt, pl in zip(tag["GT"], tag["PL"]):
                old_pl = list(map(int, pl.split(",")))
                if gt == "1/1":
                    new_pl = [ old_pl[0], old_pl[0], old_pl[-1] ]
                elif gt == "0/0":
                    new_pl = [ old_pl[0], old_pl[-1], old_pl[-1] ]
                else:
                    continue
                
                pl_lst.append(",".join(new_pl))
                second_smallest = sorted(set(new_pl))[1]
                if second_smallest > 99:
                    gq_lst.append("99")
                else:
                    gq_lst.append(str(second_smallest))
            # Mutate original tag (GT is static in the second part)
            tag["PL"], tag["GQ"] = pl_lst, gq_lst
    
        # Mutate original VCF row
        row[8] = ":".join(tag.keys())     
        for i in range(9, row.shape[0], 1):
            info = list(map(list, zip(*tag.values())))[i-9]
            row[i] = ":".join(info)
    
    return row

def checkDP(record, threshold=0, how="all"):
    """
    This function removes variant(s) with DP == 0 or DP < threshold. 
    
    It checks all subjects or just the first one.
    
    Require a valid VCF file.
    
    """
    import pandas as pd 
        
    info_lst = [record["FORMAT"].split(":")] + [record[col].split(":") for col in range(9, record.shape[0], 1)]
    tag = {info[0]: list(info[1:]) for info in zip(*info_lst)}
    
    if "DP" not in tag:
        return False
    
    if how == "all":
        return False if any(int(x) <= threshold for x in tag["DP"]) else True

    elif how == "first":
        return False if int(tag["DP"][0]) < threshold else True

        
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