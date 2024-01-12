import logging
from multiprocessing.sharedctypes import Value
from random import sample
import time
import uuid
import os
import numpy as np
import re
import pandas as pd
import sys
import argparse as ap
import math
import sqlite3
import subprocess
import io
import json
import requests
from subprocess import PIPE


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def executeCmd(cmd, stdout_only = False, logger = logger) -> None:
    if stdout_only:
        result = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    else:
        result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=PIPE)
        
    # logger.info(f"Running the following shell command inside python:\n{cmd}\nAnd the output goes like this:\n{result.stdout.decode()}\n\n")
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        if cmd_lst[1][0] != "-":
            raise RuntimeError("Error in {}:\n{}\n".format(" ".join(cmd_lst), result.stdout.decode()))
        else:
            raise RuntimeError("Error in {}:\n{}\n".format(cmd_lst, result.stdout.decode()))
    
    return result.stdout.decode()


def chi2_goodness_of_fit(data, dist=""):
    import importlib
    fit_dist = getattr(importlib.import_module("scipy.stats"), dist)
    

def na_value(v):
    # logger.info("Input value is {}, and its type is {}, whether its NAN? {}, {}".format(v, type(v), v in [np.nan], math.isnan(v)))
    if type(v) == str:
        if re.search(r"^[Nn][Aa][Nn]*$", v):
            return True
        elif re.search(r"^[;,\|]*[Nn][Aa][Nn]*[;,\|]*$", v):
            return True
        elif re.search(r"^([Nn][Aa][Nn]*[;,\|_\- ]+)+([Nn][Aa][Nn]*)*$", v):
            return True
        elif re.search(r"^[\.\-\*_ ]*$", v):
            return True
        else:
            return False
    elif v is None:
        return True
    elif v in [np.nan]:
        return True
    elif type(v) == pd._libs.tslibs.timestamps.Timestamp:
        return False
    elif math.isnan(v):
        return True
    else:
        return False



def drop_na_rows(df, subset=[], how="any"):
    # First determine the input data type of the df
    if isinstance(df, str):
        path = df
        df = pd.read_table(df, low_memory=False, encoding="utf-8")
    else:
        path = None
    
    # Then determine the columns of df
    if len(subset) == 0:
        subset = df.columns.tolist()
    elif isinstance(subset, str):
        subset = subset.strip().split(",")
        
    # Then get a boolean column determine whether a row contains all
    if how == "any":
        bool_series = df.loc[:, subset].applymap(na_value).any(axis="columns")
     
    if how == "all":
        bool_series = df.loc[:, subset].applymap(na_value).all(axis="columns")
    
    # Then drop the rows
    if isinstance(path, str):
        df.loc[np.logical_not(bool_series), :].to_csv(path, sep="\t", index=False)
    else:
        return df
     
    
    
def check_null_groupkeys(keys):
    if type(keys) == tuple or type(keys) == list:
        return all([ na_value(k) for k in keys ])
    else:
        return na_value(keys)
    
    
def deal_with_nullkey_group_return(func):
    def check_name_wrapper(*args, **kwargs):
        first_arg = [a for a in args][0]
        res = func(*args, **kwargs)
        if type(first_arg) == pd.core.frame.DataFrame:
            try:
                group_name = first_arg.name
            except AttributeError:
                logger.debug("WARNING: This function {} is supposed to be used as a groupby apply function, while we cannot access the groupdf name".format(func.__name__))
                return res
            else:  
                if check_null_groupkeys(group_name):
                    print("WARNING: Input group df to function {} has all NA group keys: {}, key type {}. The group_df has shape of {}, check the head part:\n{}".format(func.__name__, 
                                                                                                                                                            group_name, 
                                                                                                                                                            type(group_name),
                                                                                                                                                            first_arg.shape, 
                                                                                                                                                            first_arg[:3].to_string(index=False)), file=sys.stderr)
                    if type(res) == pd.core.frame.DataFrame:
                        input_cols = first_arg.columns.tolist()
                        output_cols = res.columns.tolist()
                        new_cols = [ l for l in output_cols if l not in input_cols ]
                        if len(new_cols) == 0:
                            return res
                        else:
                            res.loc[:, new_cols] = np.nan
                            return res.drop_duplicates()
                    elif res is None:
                        pass
                    else:
                        return res
                else:
                    return res
        else:
            print("WARNING: This function {} is supposed to be used as a groupby apply function, while the first arg is not a dataframe, instead its a {}".format(func.__name__, type(first_arg)), file=sys.stderr)
            return res
    return check_name_wrapper
            
        
def drop_tobe_anno_cols(df, anno_cols=[], tobe_dropped=[]):
    for col in df.columns.tolist():
        for ac in anno_cols:
            if re.search(r"^{}_*[a-z]*$".format(ac), col):
                tobe_dropped.append(col)
    
    tobe_dropped = list(dict.fromkeys(tobe_dropped))
    df.drop(columns=tobe_dropped, errors="ignore", inplace=True)
    return df            
    

def log(_func=None, *, my_logger=logger):
    import functools
    def decorator_log(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger = my_logger
            args_repr = [repr(a) for a in args]
            kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
            signature = ", ".join(args_repr + kwargs_repr)
            logger.info(f"function {func.__name__} called with args {signature}")
            start_time = time.time()
            try:
                result = func(*args, **kwargs)
                end_time = time.time()
                logger.info(f"Used {end_time - start_time} secs to run function {func.__name__}")
                return result
            except Exception as e:
                logger.exception(f"Exception raised in {func.__name__}. exception: {str(e)}")
                raise e
        return wrapper

    if _func is None:
        return decorator_log
    else:
        return decorator_log(_func)
    
    
def strip_control_characters(input_str):
    if input_str:
        import re
        # unicode invalid characters
        RE_XML_ILLEGAL = u'([\u0000-\u0008\u000b-\u000c\u000e-\u001f\ufffe-\uffff])' + \
                         u'|' + \
                         u'([%s-%s][^%s-%s])|([^%s-%s][%s-%s])|([%s-%s]$)|(^[%s-%s])' % \
                          (chr(0xd800),chr(0xdbff),chr(0xdc00),chr(0xdfff),
                           chr(0xd800),chr(0xdbff),chr(0xdc00),chr(0xdfff),
                           chr(0xd800),chr(0xdbff),chr(0xdc00),chr(0xdfff),
                           )
        input_str = re.sub(RE_XML_ILLEGAL, "", input_str)
        # ascii control characters
        input_str = re.sub(r"[\x01-\x08\x0B\x0C\x0E-\x1B]", "", input_str)
        input_str = re.sub(r"[\x0D]", "\n", input_str)
    return input_str


def batch_map_transcripts(transcript_ids):
    import requests
    import json
    from Bio import SeqIO
    from io import StringIO
    import numpy as np
    import time

    ensembl_rest_url = "http://rest.ensembl.org"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    # Attempt to fetch sequences for all transcripts in a batch query, retrying up to 10 times on network errors
    max_retries = 10
    retry_delay = 1  # seconds
    sequences = {}
    for attempt in range(max_retries + 1):
        try:
            data = json.dumps({"ids": transcript_ids})
            sequence_endpoint = f"{ensembl_rest_url}/sequence/id"
            sequence_response = requests.post(sequence_endpoint, headers=headers, data=data)
            sequence_response.raise_for_status()
            sequence_results = sequence_response.json()
            sequences = {result["id"]: result["seq"] for result in sequence_results}
            break  # break out of the loop on success
        except requests.exceptions.RequestException as e:
            logger.info(f"Error fetching sequences for batch query, attempt {attempt + 1}/{max_retries + 1}: {e}")
            if attempt < max_retries:
                time.sleep(retry_delay)

    # Initialize the dictionary for all transcripts
    transcript_data = {}
    for transcript_id in transcript_ids:
        transcript_data[transcript_id] = {
            "sequence": np.nan,
            "cds_sequence": np.nan,
            "reference_assembly": np.nan,
            "gene_symbol": np.nan,
            "ncbi_transcript_id": np.nan
        }

    # Populate the dictionary with sequences and metadata
    for transcript_id in transcript_ids:
        if transcript_id in sequences:
            transcript_data[transcript_id]["sequence"] = sequences[transcript_id]

        lookup_endpoint = f"{ensembl_rest_url}/lookup/id/{transcript_id}?expand=1;include=cds,refseq,display_name"
        try:
            lookup_response = requests.get(lookup_endpoint, headers=headers)
            lookup_response.raise_for_status()
            lookup_data = lookup_response.json()

            if "transcript_support_level" in lookup_data and lookup_data["transcript_support_level"] != "NA":
                transcript_data[transcript_id]["reference_assembly"] = lookup_data["assembly_name"]
                transcript_data[transcript_id]["gene_symbol"] = lookup_data["display_name"].split("-")[0]

            if "Transcript" in lookup_data and "Translation" in lookup_data["Transcript"]:
                cds_seq = lookup_data["Transcript"]["Translation"]["seq"]
                transcript_data[transcript_id]["cds_sequence"] = cds_seq

            if "xrefs" in lookup_data:
                for xref in lookup_data["xrefs"]:
                    if xref["dbname"] in ["RefSeq_mRNA", "RefSeq_ncRNA"]:
                        transcript_data[transcript_id]["ncbi_transcript_id"] = xref["primary_id"]
                        break
        except requests.exceptions.RequestException as e:
            logger.info(f"Error fetching metadata for {transcript_id}: {e}")
            # For missing metadata items, set the corresponding values to np.nan

    return transcript_data



def get_cds_len_from_ensembl_tranxid(tranxid="ENST00000498094", 
                                     cache_dir="/paedyl01/disk1/yangyxt/public_data/ensembl_resources/pyensembl",
                                     release=75):
    os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir
    from pyensembl import EnsemblRelease
    er = EnsemblRelease.cached(release=release)
    try:
        tranx_obj = er.transcript_by_id(tranxid)
    except ValueError as ve:
        logger.warning("Transcript {} cannot be matched".format(tranxid))
        return np.nan
    else:
        try:
            return len(tranx_obj.coding_sequence)
        except ValueError as ve:
            logger.debug("Transcript {} seems to be a non-coding transcript".format(tranxid))
            return np.nan
        
        

def prepare_ensembl_cache(cache_dir="/paedyl01/disk1/yangyxt/public_data/ensembl_resources/pyensembl", 
                          assembly = "hg19", 
                          assembly_release_map = {"hg19": 75,
                                                  "hg38": 104}):
    from pyensembl import EnsemblRelease
    os.environ["PYENSEMBL_CACHE_DIR"] = cache_dir
    
    ensembl = EnsemblRelease(assembly_release_map[assembly])
    ensembl.download()
    ensembl.index()
    
    

def cdna_to_genomic_coordinate(cdna_coordinate,
                               assembly = "hg19",
                               cache_dir = "/paedyl01/disk1/yangyxt/public_data/ensembl_resources/pyensembl",
                               assembly_release_map = {"hg19": 75,
                                                       "hg38": 104}):
    # Should input the cdna_coordinate following this format: NM_002524.5:c.35G>A
    
    from pyensembl import EnsemblRelease
    os.environ["PYENSEMBL_CACHE_DIR"] = cache_dir
    
    ensembl = EnsemblRelease.cached(assembly_release_map[assembly])

    # Parse the input cDNA coordinate
    transcript_id, _, coord = cdna_coordinate.partition(":")

    # Get the transcript object
    transcript = ensembl.transcript_by_id(transcript_id)
    print(transcript, type(transcript))
    coding_sequence = transcript.coding_sequence()
    
    # get the cdna coordinate
    ref, pos, alt = coord[0], int(coord[1:-1]), coord[-1]

    # Convert cDNA position to genomic position
    genomic_pos = transcript.spliced_offset(pos - 1)
    chrom = transcript.contig
    genomic_coordinate = f"{chrom}:{genomic_pos}:{ref}>{alt}"

    return genomic_coordinate


def create_gffutils_db(gff_file="/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.gz", 
                       db_file = None, **kwargs):
    if not db_file:
        db_file = re.sub("(\.g[tf]f[3]*)\.gz$", "\1.db", gff_file)
        logger.info(f"The output database file is {db_file}")
        assert db_file != gff_file, "The input gff file should not be compressed"
    
    if os.path.getmtime(gff_file) < os.path.getmtime(db_file):
        logger.info(f"The gff file {gff_file} is older than the db file {db_file}, skip the db creation")
        db = gffutils.FeatureDB(db_file)
    else:
        import gffutils
        # kwargs here can have merge_strategy = "merge" or merge_strategy = "create_unique"
        db = gffutils.create_db(gff_file, 
                                dbfn=db_file, 
                                force=True, 
                                keep_order=True,  
                                sort_attribute_values=False,
                                merge_strategy = "merge",
                                **kwargs)
        logger.info("Finished the db creation at {}\n".format(db_file))

    return db


def find_features_by_id_pattern(gff3_db=None, query_id=""):
    id_pattern=r"^.*{id}.*$".format(id = query_id)
    matching_features = [feature for feature in gff3_db.all_features() if re.search(id_pattern, feature.attributes["ID"][0])]
    return matching_features


def find_cds_features_by_id_pattern(gff3_db=None, query_id=""):
    id_pattern=r"^.*{id}.*$".format(id = query_id)
    matching_features = [feature for feature in gff3_db.all_features() if feature.featuretype == "CDS" and re.search(id_pattern, feature.attributes.get("Parent", ["None",])[0])]
    return matching_features


def find_cds_features_by_gene_symbol(gff3_db=None, query_id=""):
    id_pattern=r"^.*{id}.*$".format(id = query_id)
    matching_features = [feature for feature in gff3_db.all_features() if feature.featuretype == "CDS" and re.search(id_pattern, feature.attributes["gene"][0])]
    return matching_features



def extract_pseudogene_symbols(gff3_db = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db",
                               output_lst = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.pseudogenes.lst"):
    import gffutils
    db = gffutils.FeatureDB(gff3_db)
    all_pseudogenes = [(rec.attributes.get("Name", "")[0], rec.attributes.get("description", "")[0]) for rec in db.all_features() if rec.featuretype == "pseudogene"]
    all_pseudogenes = list(dict.fromkeys([x for x in all_pseudogenes if len(x[0]) > 0]))

    if output_lst:
        with open(output_lst, "w") as ol:
            ol.write("\n".join(t[0] + "\t" + t[1] for t in all_pseudogenes))
    else:
        return all_pseudogenes


        
def convert_gff_to_bed_coding(gff3_db = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db",
                              output_bed = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.coding.func.bed"):
    import gffutils
    db = gffutils.FeatureDB(gff3_db)
    # Extract mRNA transcript IDs of all the CDS regions
    all_cds_transcript_ids = [rec.attributes.get("ID", "")[0] for rec in db.all_features() if rec.featuretype == "mRNA" and re.search(r"^NC_", rec.chrom) ]
    all_cds_transcript_ids = set(list(dict.fromkeys([x.split("-")[1] for x in all_cds_transcript_ids if len(x) > 0 and re.search(r"NM_", x)])))
    # use the transcript IDs to extract the exon regions (exon includes UTR and CDS)
    all_cds_transcript_exons = [rec for rec in db.all_features() if rec.featuretype == "exon" and re.search(r"^NC_", rec.chrom) and rec.attributes.get("ID","none-none-none")[0].split("-")[1] in all_cds_transcript_ids ]
    with open(output_bed, "w") as ob:
        for exon_rec in all_cds_transcript_exons:
            try:
                chr_ind = re.search(r"^NC_0+([1-9]+[0-9]*)\.[0-9]+$", exon_rec.chrom).group(1)
            except AttributeError:
                logger.warning("This record cannot match the chromosome format: {}\n".format(exon_rec))
                raise ValueError
            sex_chr_map = {"23": "X", "24": "Y", "12920": "M"}
            chr_ind = chr_ind if chr_ind not in ["23", "24", "12920"] else sex_chr_map[chr_ind]
            chrom = "chr" + chr_ind
            start = exon_rec.start
            end = exon_rec.end
            strand = exon_rec.strand
            transcript_id = exon_rec.attributes.get("ID","none-none-none")[0].split("-")[1]
            exon_id = "exon_" + exon_rec.attributes.get("ID","none-none-none")[0].split("-")[2]
            gene_symbol = exon_rec.attributes.get("gene","")[0]
            ob.write(f"{chrom}\t{start}\t{end}\t{gene_symbol}\t{transcript_id}\t{exon_id}\t{strand}\n")
            
            

def coding_bed_of_canonical_tranx(input_bed = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.coding.func.bed",
                                  tranx_map = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.funcoding_genes.lst",
                                  output_bed = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.coding.func.canonical.bed"):
    input_df = pd.read_table(input_bed, header=None, names=['chr_exon', 'start_exon', 'end_exon',
                                                            'symbol', 'transcript_ID', 'exon_number', 'strand'])
    tranx_map_df = pd.read_table(tranx_map, header=None, names=['symbol', 'canonical_tranx_id'])
    can_tranx = set(tranx_map_df.loc[:, "canonical_tranx_id"].drop_duplicates().tolist())
    
    output_df = input_df.loc[ input_df['transcript_ID'].isin(can_tranx), : ]
    output_df.to_csv(output_bed, sep="\t", index=False, header=False)


def extract_canonical_tranx_id(gff3="/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                               output_tab="/paedyl01/disk1/yangyxt/public_data/gene_annotation/RefSeq_canonical_tranx_id.tsv"):
    db = create_gffutils_db(gff3)
    all_mrnas = [rec for rec in db.all_features() if rec.featuretype == "mRNA"]
    all_canonical_tranx = [ rec for rec in all_mrnas if "Select" in rec.attributes.get("tag", [""])[0] ]
    with open(output_tab, "w") as ot:
        for rec in all_canonical_tranx:
            symbol = rec.attributes.get("gene", [""])[0]
            tranx_id = rec.attributes.get("ID", [""])[0].split("-")[1]
            ot.write(f"{symbol}\t{tranx_id}\n")


def extract_gene_anno_for_msa(gff3_db="/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db",
                              gene_name="",
                              output_gff=""):
    import gffutils
    db = gffutils.FeatureDB(gff3_db)
    assert len(gene_name) > 0, "Please at least input a valid HGNC gene symbol name"
    
    # Build a records iterator
    rec_iterator = db.all_features()
    
    # Then start to look for the gene record
    while True:
        rec = next(rec_iterator)
        if rec.featuretype == "gene" or rec.featuretype == "pseudogene" and rec.attributes.get("Name", [None])[0] == gene_name:
            gene_rec = rec
            break
        
    # Then start to extract the mRNA record
    while True:
        rec = next(rec_iterator)
        if rec.featuretype == "mRNA" and rec.attributes.get("gene", [None])[0] == gene_name:
            mRNA_rec = rec
            break

    mRNA_name = mRNA_rec.attributes.get("Name", [None])[0]
    assert mRNA_name, f"Failed to extract the mRNA name for this record: {mRNA_rec}"

    # Then extract the exons record
    exon_recs = []
    cds_recs = []
    while True:
        rec = next(rec_iterator)
        if rec.featuretype == "exon" and rec.attributes.get("Parent", [None])[0] == f"rna-{mRNA_name}":
            exon_recs.append(rec)
        if rec.featuretype == "CDS" and rec.attributes.get("Parent", [None])[0] == f"rna-{mRNA_name}":
            cds_recs.append(rec)
        if rec.attributes.get("Parent", [None])[0] != f"rna-{mRNA_name}":
            break
    
    # Modify sequence name and coordinates
    gene_start = gene_rec.start
    for feature in [gene_rec, mRNA_rec] + exon_recs + cds_recs:
        feature.seqid = gene_name
        feature.start = feature.start - gene_start + 1
        feature.end = feature.end - gene_start + 1

    with open(output_gff, 'w') as out_gff:
        for feature in [gene_rec, mRNA_rec] + exon_recs + cds_recs:
            out_gff.write(str(feature) + '\n')    
    
    
def extract_functional_coding_symbols(gff3_db = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db",
                                      output_lst = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.funcoding_genes.lst"):
    import gffutils
    db = gffutils.FeatureDB(gff3_db)
    all_cds_symbols = [(rec.attributes.get("gene", " ")[0], rec.attributes.get("product", " ")[0]) for rec in db.all_features() if rec.featuretype == "CDS"]
    all_cds_symbols = list(dict.fromkeys([x for x in all_cds_symbols if len(x[0]) > 1]))
    
    if output_lst:
        with open(output_lst, "w") as ol:
            ol.write("\n".join([f"{t[0]}\t{t[1]}" for t in all_cds_symbols]))
    else:
        return all_cds_symbols
        

def extract_exon_regions_per_gene(gff3_db = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db",
                                chr_alias_map = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/chromosome.alias.tsv",
                                meta_table = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/all_exon_beds/GCF_000001405.25_GRCh37.p13_genomic.exon.meta.tsv"):
    import gffutils
    import pybedtools
    db = gffutils.FeatureDB(gff3_db)
    
    chr_alias_df = pd.read_table(chr_alias_map, low_memory=False)
    logger.info("The chromosome name alias map table looks like:\n{}\n".format(chr_alias_df[:5].to_string(index=False)))
    # header is ucsc, assembly, genbank, refseq

    # Extract their canonial transcript IDs
    rec_iterator = db.all_features()

    with open(meta_table, "w") as mt:
        while True:
            try:
                rec = next(rec_iterator)
            except StopIteration:
                break
            else:
                if rec.featuretype == "mRNA":
                    symbol = rec.attributes.get("gene", "")[0]
                    if len(symbol) > 0:
                        tranx_id = rec.attributes.get("ID", "-")[0].split("-")[1]
                        tranx_bed = meta_table.replace(".exon.meta.tsv", f".{tranx_id}.exon.bed")
                        intervals = []
                        while True:
                            try:
                                exon_rec = next(rec_iterator)
                            except StopIteration:
                                break
                            else:
                                if exon_rec.featuretype != "exon":
                                    break
                                else:
                                    seq_ID = exon_rec.seqid
                                    if seq_ID in chr_alias_df.loc[:, "refseq"].values:
                                        chrom = chr_alias_df.loc[chr_alias_df.loc[:, "refseq"] == seq_ID, "ucsc"].values[0]
                                    elif chr_alias_df.loc[:, "refseq"].str.contains(seq_ID.split(".")[0]).sum() > 0:
                                        chrom = chr_alias_df.loc[chr_alias_df.loc[:, "refseq"].str.contains(seq_ID.split(".")[0]), "ucsc"].values[0]
                                    else:
                                        logger.warning("Cannot find any matches in the chromosome alias map table, check the seqID: {}\n".format(seq_ID))
                                        break
                                    
                                    start = exon_rec.start - 1
                                    end = exon_rec.end
                                    strand = exon_rec.strand
                                    exon_id = exon_rec.attributes.get("ID", "-")[0].split("-")[2]
                                    interval = pybedtools.Interval(chrom, start, end, name="exon_{}".format(exon_id), score=".", strand=strand)
                                    intervals.append(interval)
                        bed_obj = pybedtools.BedTool(intervals).sort()
                        bed_obj.saveas(tranx_bed)
                        mt.write(f"{symbol}\t{tranx_id}\t{tranx_bed}\n")
                    else:
                        logger.warning("This mRNA record does not contain an attribute called gene: {}\n".format(rec))


def identify_constitutive_exons(chr, start, end,
                                meta_table="/paedyl01/disk1/yangyxt/public_data/gene_annotation/all_exon_beds/GCF_000001405.25_GRCh37.p13_genomic.exon.meta.tsv", 
                                gene_symbol = "",
                                tranx_id = "",
                                skip_long_check = False):
    # Require the start to be 0-indexed and end to be not inclusive
    import pybedtools
    meta_df = pd.read_table(meta_table, low_memory=False, header=None, names=["symbol", "tranx_id", "bed_path"])

    if gene_symbol in meta_df.loc[:, "symbol"].values:
        tranx_beds = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "bed_path"].values
        tranx_ids = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "tranx_id"].values
        gene_symbols = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "symbol"].values
    else:
        if tranx_id in meta_df.loc[:, "tranx_id"].values:
            gene_symbol = meta_df.loc[meta_df.loc[:, "tranx_id"] == tranx_id, "symbol"].values[0]
            tranx_beds = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "bed_path"].values
            tranx_ids = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "tranx_id"].values
            gene_symbols = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "symbol"].values
        elif tranx_id.split(".")[0] in meta_df.loc[:, "tranx_id"].str.split(".", expand=True).iloc[:, 0].values:
            gene_symbol = meta_df.loc[meta_df.loc[:, "tranx_id"].str.split(".", expand=True).iloc[:, 0] == tranx_id.split(".")[0], "symbol"].values[0]
            tranx_beds = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "bed_path"].values
            tranx_ids = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "tranx_id"].values
            gene_symbols = meta_df.loc[meta_df.loc[:, "symbol"] == gene_symbol, "symbol"].values
        else:
            tranx_beds = meta_df.loc[:, "bed_path"].values
            tranx_ids = meta_df.loc[:, "tranx_id"].values
            gene_symbols = meta_df.loc[:, "symbol"].values
            logger.warning("Cannot use the gene_symbol info to narrow down the range. It is strongly not recommended to compare the region with all the bed paths, its just too time expensive")
            if skip_long_check:
                return np.nan

    bed_obj = pybedtools.BedTool(f"{chr}\t{start}\t{end}", from_string=True)

    assert isinstance(bed_obj, pybedtools.bedtool.BedTool), "The input interval is not a valid pybedtools.Interval object"
    logger.debug("The prepared bedtool object for the input region looks like:\n{}\n".format(bed_obj.to_dataframe(disable_auto_names=True, names=["chr", "start", "end"]).to_string(index=False)))

    gene_symbols = list(dict.fromkeys(gene_symbols))
    contained_in_exon = {g:{} for g in gene_symbols}
    for g in gene_symbols:
        contained_in_exon[g] = {tx: False for tx in tranx_ids}
        for i, bed in enumerate(tranx_beds):
            tranx_id = tranx_ids[i]
            tranx_exon_bed = pybedtools.BedTool(bed).sort()
            assert isinstance(tranx_exon_bed, pybedtools.bedtool.BedTool), "The input bed path is not a valid pybedtools.BedTool object"
            logger.debug("The prepared bedtool object for the exon region of transcript {} looks like:\n{}\n".format(tranx_id, tranx_exon_bed.to_dataframe(disable_auto_names=True).to_string(index=False)))
            intersect_bed = bed_obj.intersect(tranx_exon_bed, wo=True)
            if len(intersect_bed) > 0:
                contained_in_exon[g][tranx_id] = True

        if all([contained_in_exon[g][tx] for tx in contained_in_exon[g]]):
            return True
        elif any([contained_in_exon[g][tx] for tx in contained_in_exon[g]]):
            logger.debug("Not all the exons in all transcripts of gene {} contain the input region {}: {}".format(g, (chr, start, end), contained_in_exon[g]))
        else:
            logger.debug("None of the exons in all transcripts of gene {} contain the input region {}: {}".format(g, (chr, start, end), contained_in_exon[g]))
            return np.nan

    return False


def cdna_to_genomic_coordinate_gffutils(cdna_coordinate, 
                                        gff3_db = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.db"):
    # "NM_006785:c.1476A>C"
    # samtools faidx <ref.fasta> chr:start-end
    
    import gffutils
    db = gffutils.FeatureDB(gff3_db)

    # Parse the input cDNA coordinate and extract the transcript ID
    match = re.match(r'^(?P<transcript_id>[A-Z\_0-9\.]+):(?P<coord>[cdgr]\.?.+)$', cdna_coordinate)
    if not match:
        raise ValueError("Invalid cDNA coordinate format: {}".format(cdna_coordinate))
    
    transcript_id = match.group('transcript_id')
    coord = match.group('coord')
    parsed_coord = parse_hgvs_variant(coord)
    print(parsed_coord)
    
    if parsed_coord:
        ref, pos, alt = parsed_coord["ref"], parsed_coord["start"], parsed_coord["alt"]
    else:
        raise ValueError("Cannot correclty parse the input variant format: {}".format(cdna_coordinate))
    
    logger.info("The reference alleles are {}, start cDNA position is {} and alternative alleles are {}, the transcript ID is {}".format(ref, pos, alt, transcript_id))

    # Find the CDS
    if re.search(r"^NM_", transcript_id):
        transcripts = find_cds_features_by_id_pattern(db, transcript_id)
    else:
        transcripts = find_cds_features_by_gene_symbol(db, transcript_id)
    print(transcripts)
    # Find the CDS regions separated by introns
    cdss = [feature for feature in transcripts if feature.featuretype == 'CDS' and re.search(f'^.*{transcript_id}[\.0-9]*$', feature.attributes.get('Parent', ["None"])[0])]
    # print(cdss)
    if cdss[0].strand == "+":
        cdss.sort(key=lambda x: x.start)
    else:
        cdss.sort(key=lambda x: x.start, reverse=True)
        
    logger.info(f"Totally {len(cdss)} CDS regions has been identified for subsequent searching")
    if len(cdss) == 0:
        raise ValueError(f"Cannot locate the CDS records by matching Parent attribute with value {transcript_id}")
    
    # Convert cDNA position to genomic position
    cdna_pos = int(pos)
    genomic_pos = None
    for cds in cdss:
        cds_length = cds.end - cds.start + 1
        logger.info("Iterating to exon ID {} (starting at {} and end at {}, in strand {}, length {})".format(cds.attributes["Parent"][0],
                                                                                                            cds.start,
                                                                                                            cds.end,
                                                                                                            cds.strand,
                                                                                                            cds_length)) 
        if cdna_pos > cds_length:
            cdna_pos -= cds_length
        else:
            logger.info("Found the exon ID {} (starting at {} and end at {}, in strand {}) overlapping with the cDNA pos {}.".format(cds.attributes["Parent"][0],
                                                                                                                                         cds.start,
                                                                                                                                         cds.end,
                                                                                                                                         cds.strand,
                                                                                                                                         cdna_pos))
            if cds.strand == '+':
                genomic_pos = cds.start + cdna_pos - 1 
            else:
                genomic_pos = cds.end - cdna_pos + 1
            break

    if genomic_pos is None:
        raise ValueError("cDNA position {} is out of range".format(pos))

    chr_ind = re.search(r"^NC_0+([1-9]+[0-9]*)\.[0-9]+$", cdss[0].chrom).group(1)
    sex_chr_map = {"23": "X", "24": "Y", "12920": "M"}
    chr_ind = chr_ind if chr_ind not in ["23", "24", "12920"] else sex_chr_map[chr_ind]
    
    chrom = "chr" + chr_ind 
    genomic_coordinate = f"{chrom}:{genomic_pos}:{ref}>{alt}"
    logger.info("The found genomic coordinate of the input variant is: {}".format(genomic_coordinate))
    return genomic_coordinate




def parse_hgvs_variant(variant):
    patterns = {
        'snv': re.compile(r'^(?P<prefix>[cgnr]\.)(?P<start>\d+)(?P<ref>[ATCGatcg])>(?P<alt>[ATCGatcg])$'),
        'del': re.compile(r'^(?P<prefix>[cgnr]\.)(?P<start>\d+)(_(?P<end>\d+))?del(?P<ref>[ATCGatcg]*)$'),
        'ins': re.compile(r'^(?P<prefix>[cgnr]\.)(?P<start>\d+)_?(?P<end>\d+)?ins(?P<alt>[ATCGatcg]*)$'),
        'dup': re.compile(r'^(?P<prefix>[cgnr]\.)(?P<start>\d+)(_(?P<end>\d+))?dup(?P<ref>[ATCGatcg]*)$'),
        'inv': re.compile(r'^(?P<prefix>[cgnr]\.)(?P<start>\d+)_?(?P<end>\d+)?inv(?P<ref>[ATCGatcg]*)$')
    }

    variant_type = None
    result = None

    for var_type, pattern in patterns.items():
        match = pattern.match(variant)
        if match:
            variant_type = var_type
            result = match.groupdict()
            break
    # print(result)
    if result:
        oligomer_type = {
            'c.': 'coding',
            'g.': 'genomic',
            'r.': 'RNA',
            'n.': 'non-coding'
        }.get(result['prefix'], 'unknown')

        if not result.get('end', None):
            result['end'] = result['start']
        # print(variant_type)
        if variant_type not in {'snv'} and (not result.get('ref', None) or not result.get('alt', None)):
            # result['alt'] = '<{}>'.format(variant_type.upper())
            result['ref'] = result['alt'][0]

        result['oligomer_type'] = oligomer_type
        return result

    return None
    
    
@log
def pickout_outliers(df, bins=None, hist_jpg=None, target_col="", dists=[], common=False, tail="right", cutoff=0.05):
    '''
    Input dataframe and target column label
    The column should only store continuous values
    return a boolean pandas series with value being the index telling which rows contain an outlier value at target column
    '''
    
    from fitter import Fitter, get_common_distributions, get_distributions
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import operator
    import math
    import importlib
    
    try:
        original_value = pd.to_numeric(df[target_col], errors="coerce", downcast="float")
        raw_data = df[target_col].map({"nan":np.nan, ".":np.nan, "-":np.nan, "NaN":np.nan}).astype(float).fillna(original_value)
    except TypeError as te:
        logger.warning("Column {} cant be converted to all float type. Check the values:\n{}".format(target_col, df[target_col].drop_duplicates()))
        raw_data = pd.to_numeric(df[target_col], errors="coerce", downcast="float")
        logger.warning("After coercing no-convertable values to NA, current data array contain {} valid values and {} NA values.".format(len(raw_data) - raw_data.isna().sum(), 
                                                                                                                                         raw_data.isna().sum()))
    
    # First check whether there are some extreme values taking a large fraction of the observations.
    if raw_data.min() == raw_data.median():
        filter_bools = raw_data > raw_data.min()
        logger.warning("This continous metric of {} has its minimum value {} taking up a large fraction of the observations, min value occurences:{}, non-NA observations:{}. Filter out the observation of minimal values before drawing out the distribution.".format(target_col, raw_data.min(), raw_data[raw_data == raw_data.min()].count(), raw_data.count()))
    elif raw_data.max() == raw_data.median():
        filter_bools = raw_data < raw_data.max()
        logger.warning("This continous metric of {} has its max value {} taking up a large fraction of the observations, max value occurences:{}, non-NA observations:{}. Filter out the observation of maximal values before drawing out the distribution.".format(target_col, raw_data.max(), raw_data[raw_data == raw_data.max()].count(), raw_data.count()))
    else:
        filter_bools = np.array([True for i in range(0, raw_data.size)])
    
    data = raw_data[filter_bools]
    logger.warning("After filtration, current data array contain {} valid values and {} NA values.".format(data.size - data.isna().sum(), data.isna().sum()))
    
    if not bins:
        # If bin number has not been specified, use Sturge's Rule to decide it.
        bins = 1 + math.log10(data.count()) * 3.322
    logger.info("Input continous value of {} are going to be sliced into {} bins for histogram and model fitting.".format(target_col, bins))
    
    if dists == ["poisson"]:
        logger.exception("The poisson distribution only has one parameter of mean value. Does not need this function to determine.")
    else:
        if len(dists) == 0:      
            if common:
                dists = get_common_distributions()
            else:
                dists = get_distributions()
            dists.remove("uniform")

        f = Fitter(data, distributions=dists, timeout=30, bins=bins)
        f.fit()
        dist_errs = [ (dist, f._bic[dist] + f._fitted_errors[dist] + f._aic[dist]) for dist,v in f.fitted_pdf.items() ]
        dist_errs.sort(key=operator.itemgetter(1))
        dist = dist_errs[0][0] # Find out the best fitted model
        top_3_fits = [t[0] for t in dist_errs[:3]]
        logger.info("This is the best fitted distribution {} and its parameters:{}. The top 3 models are: {}".format(dist, f.fitted_param[dist], top_3_fits))
        
        if hist_jpg and type(hist_jpg) == str:
            fig, ax = plt.subplots(1, sharey=True, tight_layout=True)
            ax.hist(data, range=(data.min(), data.max()), bins=bins, density=True)
            for name in top_3_fits:
                ax.plot(data, f.fitted_pdf[name], label=name)
            fig.set_tight_layout(True)
            fig.savefig(hist_jpg)  
        
        fit_dist = getattr(importlib.import_module("scipy.stats"), dist)
        logger.info("{}".format(fit_dist))
        fit_params = f.fitted_param[dist]
        
    if tail == "right":
        df["pvalue"] = np.array([1 - fit_dist.cdf(x, *fit_params) for x in raw_data])
    elif tail == "left":
        df["pvalue"] = np.array([ fit_dist.cdf(x, *fit_params) for x in raw_data])
    elif tail == "both":
        df["pvalue"] = np.array([2*(1 - fit_dist.cdf(x, *fit_params)) if x >= raw_data.mean() else 2*fit_dist.cdf(x, *f.fitted_param[dist]) for x in raw_data])
        
    df["pvalue_sig"] = df["pvalue"] <= cutoff
    return df[[target_col, "pvalue_sig"]].drop_duplicates().set_index(target_col).loc[:, "pvalue_sig"]


def chunk_read_tsv_by_group(file_path, groupby_cols=[], group_n=4000, delimiter="\t"):
    '''
    Serve as a tool to read big tsv files (with header line) by group or chunks of group.
    Return a generator of pandas DataFrame
    '''
    import numpy as np
    import pandas as pd
    
    with open(file_path, mode="r") as fp:
        header_list = fp.readline().strip("\n ").split(delimiter)
        if np.array([ type(x) == int for x in groupby_cols ]).all():
            groupby_cols = [x-1 for x in groupby_cols]
        else:
            groupby_cols = [ i for i,x in groupby_cols.enumerate() if str(x) in header_list ]
        
        # first_line = np.array(fp.readline().strip("\n").split(delimiter))
        row = 1
        while True:
            try:
                line = tuple(fp.readline().strip("\n ").split(delimiter))
            except StopIteration:
                break
            else:
                if len(line) <= 1:
                    logger.warning("Row {} seems to be the end of the tsv file {}".format(row, file_path))
                    break
                else:
                    row += 1
                    logger.debug("Read to row {} in table file {}".format(row, file_path))
                    groupby_value = np.array(line).astype(str)[groupby_cols]
                    logger.debug(groupby_value)
                n = 1
                lines = set([line])
                while True:
                    pos = fp.tell()
                    try:
                        next_line = tuple(fp.readline().strip("\n").split(delimiter))
                    except StopIteration:
                        return_df = pd.DataFrame(lines, columns=header_list)
                        logger.warning("Now we have a chunk_df containing {} rows and {} columns, continue appending another chunk_df.".format(return_df.shape[0], return_df.shape[1]))
                        yield return_df
                        break
                    else:
                        logger.debug("Read to row {} in table file {}".format(row, file_path))
                        if len(next_line) <= 1:
                            logger.warning("Row {} seems to be the end of the tsv file {}".format(row, file_path))
                            return_df = pd.DataFrame(lines, columns=header_list)
                            logger.warning("Now we have a chunk_df containing {} rows and {} columns, continue appending another chunk_df.".format(return_df.shape[0], return_df.shape[1]))
                            yield return_df
                            break
                        elif (np.array(next_line).astype(str)[groupby_cols] == groupby_value).all():
                            row += 1
                            logger.debug("Row {} seems to be the same group with the last row.".format(row))
                            lines.add(next_line)
                            continue
                        elif n >= group_n:
                            logger.warning("Row {} seems to be different group with the last row. And we already have {} groups in this chunk, iterate to next chunk.".format(row, group_n))
                            try:
                                return_df = pd.DataFrame(lines, columns=header_list).drop_duplicates()
                            except Exception as e:
                                error_lines = [l for l in lines if len(l) != len(header_list)]
                                logger.error("Run into error {}, check headers ({}) and number of lines {}, and number of error lines {}.\nError lines have {} columns and sample of number lines:\n{}".format(e, header_list, len(lines), len(error_lines), len(error_lines[0]), error_lines[0]))
                                raise e
                            logger.warning("Now we have a chunk_df containing {} rows and {} columns, continue appending another chunk_df.".format(return_df.shape[0], return_df.shape[1]))
                            fp.seek(pos)
                            yield return_df
                            break
                        else:
                            n += 1
                            row += 1
                            logger.info("Row {} seems to be different group with the last row. Now we have {} groups in this chunk_df".format(row, n))
                            # logger.warning("Groupby value changed from {} to {}.".format(groupby_value, next_line.astype(str)[groupby_cols]))
                            groupby_value = np.array(next_line).astype(str)[groupby_cols]
                            lines.add(next_line)
                            continue
                        
                        
def iterate_xml(xmlfile, tag=None):
    from lxml import etree as et
    '''
    Seems it can only be used to parse plain text XML file.
    '''
    if tag:
        doc = et.iterparse(xmlfile, events=('start', 'end'), tag=tag, encoding="utf-8")
    else:
        doc = et.iterparse(xmlfile, events=('start', 'end'), encoding="utf-8")
    _, root = next(doc)
    start_tag = None
    for event, element in doc:
        if event == 'start' and start_tag is None:
            start_tag = element.tag
        if event == 'end' and element.tag == start_tag:
            yield element
            start_tag = None
            root.clear()     
            
            
def read_google_sheets_to_df(sheet_id="184Ne4LoZdJida-hnXVLVoDpjd_L9I7U8rmb00E5-Fts", 
                             service_account_file = "/paedyl01/disk1/google_cloud/hkupid-352116-2ceb75fa01b4.json", 
                             header = 0, 
                             sheet_name = "", 
                             output = None):
    '''
    To use this function, you need to first install gdrive
    and configure the OUAth (service account) for your google drive account
    '''
    import pygsheets as pgs
    gc = pgs.authorize(service_file=service_account_file)
    ss = gc.open_by_key(sheet_id)
    if len(sheet_name) > 0:
        ws = ss.worksheet_by_title(sheet_name)
    else:
        ws = ss.sheet1
    # Extract all values, convert to 2d numpy array, convert to dataframe
    np_array = np.array(ws.get_all_values())
    if type(header) == int or float:
        df = pd.DataFrame(np_array[(header + 1):], columns=np_array[header]).replace("", np.nan).dropna(how="all").drop(columns=[""])
    elif type(header) == list:
        df = pd.DataFrame(np_array, columns=header).replace("", np.nan).dropna(how="all").drop(columns=[""])
        
    if output:
        df.to_csv(output, index=False, sep="\t")
    else:
        return df
    
    
def md5(fpath):
    import hashlib
    if type(fpath) == bytes:
        logger.info("Directly input bytes content instead of file path")
        print(hashlib.md5(fpath).hexdigest())
        return hashlib.md5(fpath).hexdigest()
    else:
        hash_md5 = hashlib.md5()
        if not os.path.exists(fpath):
            return str(uuid.uuid4())
        
        with open(fpath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
                
        print(hash_md5.hexdigest())
        return hash_md5.hexdigest()
    
    
def get_IUIS_genes():
    from prepare_updated_IUIS_genelist import UPDATED_IUIS
    updated_pid_genes=UPDATED_IUIS().gene_list
    print(",".join(updated_pid_genes))
    return updated_pid_genes
    

def match_PID_genes(gene_lst="", updated_pid_genes=[]):
    if len(updated_pid_genes) == 0: 
        from prepare_updated_IUIS_genelist import UPDATED_IUIS
        updated_pid_genes=UPDATED_IUIS().gene_list
    elif type(updated_pid_genes) == str:
        updated_pid_genes = updated_pid_genes.split(",")

    if type(gene_lst) == str:
        if os.path.exists(gene_lst):
            with open(gene_lst, "r") as gl:
                genes = [ l.strip("\n") for l in gl.readlines() ]
            gl.close()
        else:
            # Delimiter should be comma if input a list
            genes = gene_lst.split(",")
    elif type(gene_lst) == list:
        genes = [ g for g in gene_lst if not na_value(g) ]
        assert len(genes) > 0
    else:
        logger.error("Input gene_list is {}, and it looks like {}. It is not what is expected for this function's input. Quit with error".format(type(gene_lst), gene_lst))
        raise TypeError
        
    included_pid_genes = [ g for g in genes if g in updated_pid_genes ]
    if len(included_pid_genes) > 0:
        print("True")
        return True
    else:
        print("False")
        return False


def create_sqlite_db(dbpath=":memory:", force=False):
    from pathlib import Path
    if os.path.exists(dbpath):
        if force:
            os.remove(dbpath)
            Path(dbpath).touch()
    else:
        Path(dbpath).touch()
            
    conn = None
    try:
        conn = sqlite3.connect(dbpath, timeout=360)
        logger.debug("Establish connection with database at {} with sqlite3 version {}".format(dbpath,sqlite3.sqlite_version))
    except Exception as e:
        logger.info("Establish connection with database at {} with sqlite3 version {}".format(dbpath,sqlite3.sqlite_version))
        logger.error(e)
        raise e
    finally:
        return conn
    

def create_sql_tab(tab_name, db_conn=None, pk = "primary_key", keys=[""], **kwargs):
    tab_create_sql = " CREATE TABLE IF NOT EXISTS \"{}\" (".format(tab_name) + \
                     ", ".join([ "\"{}\" {}".format(k,v) for k,v in kwargs.items() ]) + \
                     ", CONSTRAINT {} PRIMARY KEY (".format(pk) + \
                     ", ".join(keys) + "));"
                     
    logger.info("Composed this SQL query clause: \n{}".format(tab_create_sql))
    if db_conn is not None:
        cursor = db_conn.cursor()
        try:
            cursor.execute(tab_create_sql)
        except Exception as e:
            logger.error(e)
            raise e
        db_conn.commit()
        cursor.close()
        return db_conn
    else:
        return tab_create_sql
    
    
def compose_mod_sql_query(tabname="", keys=[], **kwargs):
    '''
    kwargs should be pairs that k is the name of var and v is the value of var
    '''
    mod_sql_query = "INSERT INTO \"{}\" ( ".format(tabname) + \
                    ", ".join([ k for k,v in kwargs.items() ]) + \
                    ") VALUES (" + \
                    ", ".join([ "\"{}\"".format(v) if type(v) == str and v != "NULL" else "{}".format(v) if not na_value(v) else "NULL" for k,v in kwargs.items() ]) + \
                    ") ON CONFLICT ({}) DO UPDATE SET ".format(", ".join(keys)) + \
                    ", ".join([ "{}=excluded.{}".format(k, k) for k,v in kwargs.items() ]) + \
                    ";"
    logger.info("Composed this SQL query clause: \n{}".format(mod_sql_query))
    return mod_sql_query


def output_sql_to_tab(db_path="", tab_name="", output_tab = "", header=True):
    from sqlalchemy import create_engine
    engine = create_engine("sqlite+pysqlite:///{}".format(db_path), echo=True)
    df = pd.read_sql_table(table_name=tab_name, con=engine)
    df.to_csv(output_tab, header=header, index=False, sep="\t")
    return df
    

def multi_thread_write_to_one_tab(total_tab_path="", 
                                  tab_columns="", 
                                  key_columns="", 
                                  dtypes="", 
                                  input_tab=""):
    '''
    Remember that the total_tab_path is never the file being updated by multiple threads simutaneously. Instead, the corresponding sql database is.
    '''
    # First read in the input table
    input_df = pd.read_table(input_tab, nrows=2, low_memory=False, na_values=['nan', 'NAN', 'NaN', '.', '-', '_', ',', ';'])
    # Preprocess input args
    db_path = ".".join(total_tab_path.split(".")[:-1]) + ".sql"
    if len(tab_columns) == 0:
        tab_columns = input_df.columns.tolist()
    else:
        tab_columns = tab_columns.split(",")
    keys = key_columns.split(",")
    tab_name = os.path.basename(total_tab_path).split(".")[0]
    tmp_path = total_tab_path + "." + str(uuid.uuid4())
    
    # Create the db file, return connection
    conn = create_sqlite_db(db_path)
    
    # Create the table in db (if not exists), return connection
    if len(dtypes) == 0:
        df_dtypes = input_df.dtypes.tolist()
        dtype_map = { np.dtype('O'): "TEXT",
                      np.dtype('int64'): "INTEGER",
                      np.dtype('float64'): "REAL" }
        dtypes = [ dtype_map[x] for x in df_dtypes ]
        logger.info("The guessed Datatypes are {}".format(dtypes))
    else:
        dtypes = dtypes.split(",")
        
    assert len(dtypes) == len(tab_columns)
    dtype_map = { tab_columns[i]: dtypes[i] for i in range(len(tab_columns)) }
    conn = create_sql_tab(tab_name=tab_name,
                          db_conn=conn,
                          keys=keys,
                          **dtype_map)
    
    # Create insert clause
    with open(input_tab, "r") as it:
        first_line = it.readline().strip("\n")
    if set(first_line.split("\t")) == set(tab_columns):
        insert_df = pd.read_table(input_tab)
    else:
        insert_df = pd.read_table(input_tab, names=tab_columns)
    cursor = conn.cursor()
    for i in range(0, len(insert_df)):
        value_map = { c: "NULL" if na_value(insert_df.iat[i, insert_df.columns.get_loc(c)]) else insert_df.iat[i, insert_df.columns.get_loc(c)] for c in tab_columns }
        sql_clause = compose_mod_sql_query(tabname = tab_name,
                                        keys = keys,
                                        **value_map)
        try:
            cursor.execute(sql_clause)
        except Exception as e:
            logger.error(e)
            raise e
        conn.commit()
    
    cursor.close()
    
    # Output the result to a table
    pd.read_sql_query("SELECT * FROM \"{}\"".format(tab_name), conn).replace([None], np.nan).to_csv(tmp_path, sep="\t", index=False)
    result = subprocess.run("bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh display_table {tmp} && mv {tmp} {p} && ls -lh {p}".format(tmp = tmp_path, 
                                                                                                                                                     p = total_tab_path), 
                            shell=True, 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.STDOUT, 
                            encoding="utf-8")
    logger.info(result.stdout)
        

def convert_input_value(v):
    if type(v) == str:
        if re.search(r"^[Tt][Rr][Uu][Ee]$", v):
            return True
        elif re.search(r"^[Ff][Aa][Ll][Ss][Ee]$", v):
            return False
        elif re.search(r"^[0-9]+$", v):
            return int(v)
        elif re.search(r"^[0-9]*\.[0-9]+$", v):
            return float(v)
        elif v == "None":
            return None
        elif re.search(r"^[Nn][Aa][Nn]$", v):
            return np.nan
        else:
            return v
    else:
        return v
        
        
def create_index_sql(table_name, index_cols, index_name = "idx_of_cols", sql_path = "", db_conn = None):
    if db_conn is None:
        if len(sql_path) > 0:
            db_conn = create_sqlite_db(sql_path)
        else:
            logger.error("The input connection is None and the input DB path is none.")
            raise ValueError
        
    create_index_clause = '''
                          CREATE UNIQUE INDEX IF NOT EXISTS "{idn}" ON {tabn}({cl});
                          '''.format(idn = index_name,
                                     tabn = table_name, 
                                     cl = ", ".join(index_cols))
                          
    check_index_clause = 'PRAGMA index_list("{}");'.format(table_name)
    
    if db_conn is not None:
        cur = db_conn.cursor()
        idx_names = [ idx[1] for idx in cur.execute(check_index_clause).fetchall() ]
        if index_name in idx_names:
            return db_conn
        else:
            cur = cur.execute(create_index_clause)
            db_conn.commit()
            cur.close()
            return db_conn
            
    else:
        logger.error("The input connection is None and the input DB path is none.")
        raise ValueError


def cal_heterozygosity_rate(table_name = "gnomAD_total_variants",
                            gnomAD_sql_template="/paedyl01/disk1/yangyxt/public_data/gnomAD/total_vcf_db/gnomAD_total_variants.v3.1.2.hg19.{}.sql",
                            contigs=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                    "chr20", "chr21", "chr22", "chrX", "chrY"], 
                            output_bed = ""):
    import multiprocessing as mp
    sql_paths = [ gnomAD_sql_template.format(contig) for contig in contigs ]
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(cal_heterozygosity_per_table, sql_paths)
        result_df = pd.concat(results)



def cal_heterozygosity_per_table(sql_path_name, table_name = "gnomAD_total_variants"):
    import pandas as pd
    import sqlite3
    con = sqlite3.connect(sql_path_name)
    df = pd.read_sql_query("SELECT * FROM {}".format(table_name), con)
    if df.loc[:, "CHROM"].values[0] == "chrX":
        df["heterozygosity"] = (df["AC_Controls_XX"] - 2 * df["nhomalt_Controls_XX"])/df["AN_Controls_XX"]
    elif df.loc[:, "CHROM"].values[0] == "chrY":
        df["heterozygosity"] = 0.0
    else:
        df["heterozygosity"] = (df["AC_Controls"] - 2 * df["AN_Controls"])/df["AN_Controls"]
    return df

    
def check_var_presence( input_variant,
                        table_name = "gnomAD_total_variants",
                        gnomAD_sql_template="/paedyl01/disk1/yangyxt/public_data/gnomAD/total_vcf_db/gnomAD_total_variants.v3.1.2.hg19.{}.sql",
                        contigs=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                 "chr20", "chr21", "chr22", "chrX", "chrY"],
                        target_cols = ["AC_Controls", "nhomalt_Controls",
                                       "AC_Controls_XX", "nhomalt_Controls_XX",
                                       "AC_Controls_XY", "nhomalt_Controls_XY"],
                        match_cols = {"CHROM": "Chr", "POS": "Start", "REF": "Ref", "ALT": "Alt"}):
    import vcfpy
    if type(input_variant) == vcfpy.Record:
        contig = input_variant.CHROM
        pos = input_variant.POS
        ref = input_variant.REF
        alt = input_variant.ALT[0].serialize()
    elif type(input_variant) == pd.core.series.Series:
        col_names = input_variant.index.tolist()
        if (all([v in col_names for k,v in match_cols.items()])):
            contig = input_variant[match_cols["CHROM"]]
            pos = input_variant[match_cols["POS"]]
            ref = input_variant[match_cols["REF"]]
            alt = input_variant[match_cols["ALT"]]
        else:
            logger.warning("The input pandas row does not match with Chr,Ref,Alt,Start. Let's take a look at the input type: {}, and content: {}".format(type(input_variant), input_variant))
            match_cols = { ["CHROM", "POS", "REF", "ALT"][i]:col_names[[0,1,3,4]][i] for i in range(0,4) }
            contig = input_variant[match_cols["CHROM"]]
            pos = input_variant[match_cols["POS"]]
            ref = input_variant[match_cols["REF"]]
            alt = input_variant[match_cols["ALT"]]
    else:
        # logger.debug("The input variant record is an instance of {}, show the content:\n{}".format(type(input_variant), input_variant))
        contig = input_variant[match_cols["CHROM"]]
        pos = input_variant[match_cols["POS"]]
        ref = input_variant[match_cols["REF"]]
        alt = input_variant[match_cols["ALT"]]
            
    if contig in contigs:
        gnomAD_sql_path = gnomAD_sql_template.format(contig)
        conn = create_index_sql(table_name=table_name,
                                index_cols=["CHROM", "POS", "REF", "ALT"],
                                sql_path = gnomAD_sql_path)
        if conn is not None:
            cursor = conn.cursor()
            query_clause = '''
                           SELECT {tc} 
                           FROM "{tn}"
                           WHERE "CHROM" = "{cg}" AND "POS" = {ps} AND "REF" = "{rf}" AND "ALT" = "{at}";  
                           '''.format(tc = ", ".join([ '"' + c + '"' for c in target_cols ]),
                                      tn = table_name,
                                      cg = contig,
                                      ps = pos,
                                      rf = ref,
                                      at = alt)
            logger.debug("Before executing this clause, check the SQL phrase: \n{}\n".format(query_clause))
            try:
                results = cursor.execute(query_clause).fetchall()
            except sqlite3.OperationalError as se:
                logger.error("This SQL clause caused error: {}\n".format(query_clause))
                raise se
            if len(results) == 0:
                if type(input_variant) == pd.core.series.Series:
                    for tc in target_cols: input_variant[tc] = np.nan
                    return input_variant
                else: 
                    return np.nan
            elif len(results) > 1:
                logger.warning("This variant input matches with multiple records in gnomAD database: {}\n The matched records are: \n{}\n".format(input_variant, results))
                if type(input_variant) == pd.core.series.Series:
                    for i in range(0, len(target_cols)): input_variant[target_cols[i]] = results[0][i]
                    return input_variant
                else:
                    return results[0]
            else:
                if type(input_variant) == pd.core.series.Series:
                    for i in range(0, len(target_cols)): input_variant[target_cols[i]] = results[0][i]
                    return input_variant
                else:
                    return results[0]
        else:
            if type(input_variant) == pd.core.series.Series:
                for tc in target_cols: input_variant[tc] = np.nan
                return input_variant
            else: 
                return np.nan
    else:
        if type(input_variant) == pd.core.series.Series:
            for tc in target_cols: input_variant[tc] = np.nan
            return input_variant
        else: 
            return np.nan
        

def check_column_empty(path="", labels=""):
    import pandas as pd
    import numpy as np
    try:
        df = pd.read_table(path, low_memory=False, usecols=labels.split(","), na_values=["", ".", "-"])
    except ValueError as ve:
        logger.info("Some of the input columns not found in the table: {}, running into error {}".format(labels, ve))
        print("NOT_FOUND")
    if len(df) == 0:
        print("ALL_EMPTY")
        return
    # print(df.to_string(), file=sys.stderr)
    final_bool = True
    for col in df.columns.tolist():
        if df[col].apply(na_value).all():
            final_bool = final_bool and True
        else:
            # print(df.loc[df[col].apply(na_value), :].to_string(), file=sys.stderr)
            final_bool = final_bool and False
    if final_bool:
        print("ALL_EMPTY")
    else:
        print("NOT_ALL_EMPTY")
    

def get_latest_ncbi_anno(output, update=True):
    from updated_ncbi_annotation import UPDATED_NCBI_ANNO
    UPDATED_NCBI_ANNO(parallel=True, output = output, update=update).extract_cds()
    

def combine_tables_to_excel(excel_path="", 
                            table_paths="",
                            sheet_labels=""):
    sheet_labels=sheet_labels.split(",")
    table_paths=table_paths.split(",")
    
    if len(sheet_labels) < len(table_paths):
        sheet_labels = [ os.path.basename(table_paths[i]) for i in range(0, len(table_paths)) ]
    
    with pd.ExcelWriter(excel_path) as writer:
        for i in range(0, len(table_paths)):
            df = pd.read_table(table_paths[i], low_memory=False, encoding="utf-8")
            sheet_name = sheet_labels[i]
            if len(sheet_name) > 30: sheet_name = sheet_name[:29]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
        writer.close()
        
        
def fetch_sample_metadata(sample_ID,
                          family_ID="",
                          seq_type="",
                          total_ped="/paedyl01/disk1/yangyxt/ngs_total_ped.ped"):
    total_ped_df = pd.read_table(total_ped, low_memory=False, encoding="utf-8")
    if len(seq_type) > 0:
        ped_df = total_ped_df.loc[total_ped_df["Seq_type"] == seq_type.lower(), :]
    else:
        ped_df = total_ped_df
        
    assert type(sample_ID) == str and len(sample_ID) > 0
    
    sample_df = ped_df.loc[ped_df["IndividualID"] == sample_ID, :].drop_duplicates()
    
    if len(sample_df) == 0:
        logger.error("Cannot locate the sample ID {} in total pedigree table {}".format(sample_ID,
                                                                                        total_ped))
        raise ValueError
    elif len(sample_df) > 1:
        logger.error("Locate multiple samples with the same sample ID {} in total pedigree table {}. Lets take a look at the located pedigree records now:\n{}".format(sample_ID,
                                                                                                                    total_ped,
                                                                                                                    sample_df.to_string(index=False)))
        # Now we need to use other index to limit the matching records
        if len(seq_type) > 0:
            sample_df = sample_df.loc[sample_df["Seq_type"] == seq_type, :]
        if len(sample_df) > 1:
            if len(family_ID) > 0 and (sample_df["#FamilyID"] == family_ID).sum() > 0:
                sample_df = sample_df.loc[sample_df["#FamilyID"] == family_ID, :]
        
        if len(sample_df) == 0:
            logger.error("Cannot locate the sample ID {}, (with familyID {} and seqtype {} ) in total pedigree table {}".format(sample_ID,
                                                                                                                                        family_ID,
                                                                                                                                        seq_type, 
                                                                                                                                        total_ped))
            return np.nan
        elif len(sample_df) > 1:
            logger.error("Locate multiple samples with the same sample ID {} (with familyID {} and seqtype {}) in total pedigree table {}. Lets take a look at the located pedigree records now:\n{}".format(sample_ID,
                                                                                                                                                                                                                    family_ID,
                                                                                                                                                                                                                    seq_type, 
                                                                                                                                                                                                                    total_ped,
                                                                                                                                                                                                                    sample_df.to_string(index=False)))
            return np.nan
    
    assert len(sample_df) == 1
    return sample_df.squeeze().to_dict()


def fetch_latest_file(file_list):
    if type(file_list) == str:
        file_list = file_list.split(",")
    
    try:
        assert type(file_list) == list
    except AssertionError as ae:
        logger.warning("Please input the comma delimited file paths or a python list object")
        raise ae

    valid_file_list = []
    for f in file_list:
        if not os.path.exists(f):
            logger.debug("File {} not existed.".format(f))
        else:
            valid_file_list.append(f)
    
    # Now trying to get the latest file
    if len(valid_file_list) == 0:
        logger.warning("It seems that the input file list {} does not have any one of them existed.".format(",".join(file_list)))
        return np.nan
    
    latest_file = max(valid_file_list, key=os.path.getmtime)
    print(latest_file)
    return latest_file
        
            

def fetch_latest_anno_path(sample_ID,
                           family_ID="",
                           seq_type="",
                           total_ped="/paedyl01/disk1/yangyxt/ngs_total_ped.ped",
                           probe_col = "Probe",
                           seq_type_col = "Seq_type",
                           sample_batch_col = "Sample_Batch"):
    sample_meta = fetch_sample_metadata(sample_ID, family_ID, seq_type, total_ped)
    logger.debug("The fetched sample metadata is {}".format(sample_meta))
    probe = sample_meta[probe_col]
    seqtype = sample_meta[seq_type_col]
    family_ID = sample_meta["#FamilyID"]
    sample_batch = sample_meta[sample_batch_col]
    
    # List out the possible paths
    candidate1="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{famID}/{famID}.hg19_multianno.lowfreq.filtered.reformat.domain.allvars.rmfit.hetcomp.clnvar.roh.podel.func.tsv".format(seq_type = seqtype,
                                                                                                                                                                                                       sample_batch = sample_batch,
                                                                                                                                                                                                       famID = family_ID)
    candidate2="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{famID}/{famID}.hg19_multianno.lowfreq.filtered.reformat.domain.allvars.rmfit.hetcomp.clnvar.roh.podel.func.txt".format(seq_type = seqtype,
                                                                                                                                                                                                       sample_batch = sample_batch,
                                                                                                                                                                                                       famID = family_ID)
    candidate3="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{proband}.{famID}.hg19_multianno.reformat.mq.filtered.lowfreq.exonsplic.podel.final.txt".format(seq_type = seqtype,
                                                                                                                                                                               sample_batch = sample_batch,
                                                                                                                                                                               famID = family_ID, 
                                                                                                                                                                               proband = sample_ID)
    candidate4="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{proband}.{famID}.hg19_multianno.reformat.mq.filtered.lowfreq.exonsplic.podel.final.tsv".format(seq_type = seqtype,
                                                                                                                                                                               sample_batch = sample_batch,
                                                                                                                                                                               famID = family_ID, 
                                                                                                                                                                               proband = sample_ID)
    candidate5="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{proband}.{famID}.hg19_multianno.reformat.mq.filtered.lowfreq.podel.final.txt".format(seq_type = seqtype,
                                                                                                                                                                     sample_batch = sample_batch,
                                                                                                                                                                     famID = family_ID, 
                                                                                                                                                                     proband = sample_ID)
    candidate6="/paedyl01/disk1/yangyxt/{seq_type}/{sample_batch}/annotations/{proband}.{famID}.hg19_multianno.reformat.mq.filtered.lowfreq.podel.final.tsv".format(seq_type = seqtype,
                                                                                                                                                                     sample_batch = sample_batch,
                                                                                                                                                                     famID = family_ID, 
                                                                                                                                                                     proband = sample_ID)
    
    candidates = [candidate1, candidate2, candidate3, candidate4, candidate5, candidate6]
    updated_anno_path = fetch_latest_file(candidates)
    
    if na_value(updated_anno_path):
        logger.warning("Sample {} seems to not have a valid bioinformatic annotation result judging from the file path candidates: {}".format(sample_ID,
                                                                                                                                              ",".join(candidates)))
        return np.nan
    else:
        logger.info("Sample {} does have a valid bioinformatic annotation result path: {}".format(sample_ID, updated_anno_path))
    
    return updated_anno_path



def locate_variant_in_sample_bioinfo_result(varID="",
                                            varID_col="uniq_ID",
                                            rank_score_col="ACMG_quant_score",
                                            varcoord={"Chr":"", "Start":0, "End":0, "Ref":"N", "Alt":"N"}, 
                                            sample_ID="",
                                            seq_type="",
                                            family_ID=""):
    # First fetch the latest sample annotation path
    anno_path = fetch_latest_anno_path(sample_ID, family_ID=family_ID, seq_type=seq_type)
    
    if anno_path in [np.nan]:
        logger.warning(" Cannot fetch a valid annotation file for sample {} in family {} in sequencing sample type {}".format(sample_ID,
                                                                                                                             family_ID, 
                                                                                                                             seq_type))
        return "Cannot locate the func.tsv file"
    
    excel_path = anno_path[:-4] + ".selected.xlsx"
    if not os.path.exists(excel_path):
        logger.error("Cannot locate the updated excel file for the sample {} in family {} in sequencing sample type {}".format(sample_ID,
                                                                                                                             family_ID, 
                                                                                                                             seq_type))
        return "Cannot locate the excel file"

    vardf_dict = pd.read_excel(excel_path, sheet_name=None)
    
    # Set up the initiate value for topper ranked variants
    topper_vars = []
    for sheet_name, anno_df in vardf_dict.items():
        if rank_score_col not in anno_df.columns.tolist():
            logger.warning("Seems that this table {} in {} does not contain the ACMG rank score column {}. Let's take a look at it: \n{}\n".format(sheet_name,
                                                                                                                                                   excel_path,
                                                                                                                                                   rank_score_col, 
                                                                                                                                                   anno_df[:5].to_string(index=False)))
            continue
        # Check how to match with the variant, either using variant ID or the variant coordinates combo
        if len(varID) > 0:
            rank_score = anno_df.loc[anno_df[varID_col] == varID, rank_score_col]
        else:
            coord_bools = np.array([True for i in range(0, len(anno_df))])
            for field, value in varcoord.items():
                if value == 0:
                    continue
                if field in ["Start", "End"]:
                    field_bools = anno_df[field].astype(int) == int(value)
                else:
                    field_bools = anno_df[field].astype(str) == str(value)
                coord_bools = coord_bools & field_bools
            rank_score = anno_df.loc[coord_bools, rank_score_col]
        
        # Considering the extreme case that the variant cannot be matched with the provided variant ID or the coordinates
        if len(rank_score) == 0:
            if len(varID) > 0:
                logger.debug("Cannot match the variant () in the table {} in excel file {}".format(varID, sheet_name, excel_path))
            else:
                logger.debug("Cannot match the variant () in the table {} in excel file {}".format(varcoord, sheet_name, excel_path))
        else:
            topper_vars = anno_df.loc[anno_df[rank_score_col] > rank_score.tolist()[0], varID_col]
            break

    if type(topper_vars) == list and len(topper_vars) == 0:
        logger.warning("Cannot match the variant () in excel file {}".format(varID if len(varID) == 0 else varcoord, excel_path))
        return "Cannot locate the variant in excel file"
    
    # Proceed with matched variant records
    if len(topper_vars) == 0:
        rank = 1
    else:
        rank = len(topper_vars.drop_duplicates().tolist()) + 1
    
    return rank


def get_worksheet_object(file_id,
                        service_account_file = "/paedyl01/disk1/google_cloud/hkupid-352116-2ceb75fa01b4.json",
                        save_file = None,
                        sheet_name = ""):
    import pygsheets as pgs
    gc = pgs.authorize(service_file=service_account_file)
    ss = gc.open_by_key(file_id)
    if len(sheet_name) > 0:
        ws = ss.worksheet_by_title(sheet_name)
    else:
        ws = ss.sheet1
        
    if save_file:
        np_array = np.array(ws.get_all_values())
        np.savetxt(save_file, np_array, fmt="%s", delimiter="\t")
    return ws


def get_updated_worksheet_as_df(file_id, save_file=None, **kwargs):
    ws = get_worksheet_object(file_id, **kwargs)
    np_array = np.array(ws.get_all_values())
    df = pd.DataFrame(np_array[1:], columns=np_array[0]).replace("", np.nan).dropna(how="all")
    if save_file:
        df.to_csv(save_file, sep="\t", index=False)
    return df
    

def get_latest_search_request(sheet_id = "12cIhykbPOR0ynWyAJNcTc0RowfMvB1MQqa4h0Ia6K3w",
                              value_cols = ["Email Address", 
                                            "Who's asking for it ?",
                                            "What type of data you would like to search?",
                                            "Please input the gene symbol you would like to search",
                                            "Do you need to search variants in the entire callset?", 
                                            "Is there a specific CDS change you would like to search?",
                                            "Is there a specific amino acid change you would like to search?",
                                            "Chromosome",
                                            "Variant Start position", 
                                            "Variant End position",
                                            "Variant position reference allele",
                                            "Variant position alternative allele"],
                              **kwargs):
    df = get_updated_worksheet_as_df(sheet_id, **kwargs)
    logger.info("The fetched updated search request sheet looks like: \n{}\n".format(df[:5].to_string(index=False)))
    # Generate a uniq ID for each query
    query_fields = df.loc[:, "Email Address"]
    for col in value_cols:
        query_fields = query_fields + "," + df[col].astype(str)
    df = df.assign(query_ID = query_fields)
    return df
    

def track_search_request_progress(sheet_id = "12cIhykbPOR0ynWyAJNcTc0RowfMvB1MQqa4h0Ia6K3w",
                                  sheet_name = "Query_progress_track"):
    from datetime import datetime
    # First fetch the progress sheet
    progress_df = get_updated_worksheet_as_df(sheet_id, sheet_name = sheet_name)
    # Then fefch the latest search request
    df = get_latest_search_request(sheet_id)
    
    # Extract all the query IDs and find out whether the searching results are ready or not.
    qids = df.loc[:, "query_ID"].drop_duplicates().tolist()
    p_qids = progress_df.loc[:, "query_ID"].drop_duplicates().tolist()
    
    new_queries = [q for q in qids if q not in p_qids]
    if len(new_queries) == 0:
        logger.info("No new queries at this moment: {}".format(datetime.now()))
    
    unfinished_queries = progress_df.loc[progress_df.loc[:, "result_URL"].isna(), "query_ID"].drop_duplicates().tolist()
    if len(unfinished_queries) == 0:
        logger.info("No unfinished queries at this moment: {}".format(datetime.now()))
    
    total_need_work_queries = new_queries + unfinished_queries
    if len(total_need_work_queries) == 0:
        logger.info("No new queries or unfinished queries at this moment: {}".format(datetime.now()))
        exit(0)
    
    return total_need_work_queries, progress_df.loc[np.logical_not(progress_df.loc[:, "result_URL"].isna()), ["query_ID", "result_URL"]]


def search_per_request(query_ID,
                       sheet_id = "12cIhykbPOR0ynWyAJNcTc0RowfMvB1MQqa4h0Ia6K3w",
                       sheet_name = "Form Responses 1",
                       value_cols =["Email Address", 
                                    "Who's asking for it ?",
                                    "What type of data you would like to search?",
                                    "Please input the gene symbol you would like to search",
                                    "Do you need to search variants in the entire callset?", 
                                    "Is there a specific CDS change you would like to search?",
                                    "Is there a specific amino acid change you would like to search?",
                                    "Chromosome",
                                    "Variant Start position", 
                                    "Variant End position",
                                    "Variant position reference allele",
                                    "Variant position alternative allele"]):
    df = get_latest_search_request(sheet_id, sheet_name = sheet_name)
    requests = df.loc[df.loc[:, "query_ID"] == query_ID, value_cols]
    
    raw_callset = True if re.search(r"[Yy][Ee][Ss]", str(requests.loc[:, "Do you need to search variants in the entire callset?"].values[0])) else False
    callset_arg = " -r True" if raw_callset else ""
    
    seq_types = requests.loc[:, "What type of data you would like to search?"].values[0].split(",")
    email = requests.loc[:, "Email Address"].values[0]
    name = requests.loc[:, "Who's asking for it ?"].values[0]
    gene = requests.loc[:, "Please input the gene symbol you would like to search"].values[0]
    
    coding_change = requests.loc[:, "Is there a specific CDS change you would like to search?"].values[0]
    cds_arg = " -c {}".format(coding_change) if not na_value(coding_change) else ""
    
    aa_change = requests.loc[:, "Is there a specific amino acid change you would like to search?"].values[0]
    aa_arg = " -a {}".format(aa_change) if not na_value(aa_change) else ""
    
    chromosome = requests.loc[:, "Chromosome"].values[0]
    start = requests.loc[:, "Variant Start position"].values[0]
    end = requests.loc[:, "Variant End position"].values[0]
    ref_allele = requests.loc[:, "Variant position reference allele"].values[0]
    alt_allele = requests.loc[:, "Variant position alternative allele"].values[0]
    if all([na_value(x) for x in [chromosome, start, end, ref_allele, alt_allele]]):
        logger.info("This query {} does not try to search a variant by coordinates".format(query_ID))
        if all([na_value(x) for x in [coding_change, aa_change]]):
            logger.info("This query {} does not try to search a specific variant by specifying CDS change and AA change".format(query_ID))
        elif na_value(coding_change):
            logger.info("This query {} does try to search a specific variant by specifying AA change {}".format(query_ID, aa_change))
        else:
            logger.info("This query {} does try to search a specific variant by specifying CDS change {}".format(query_ID, coding_change))
        
    cmd = "bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh search_genes_in_all_pat -g {} -e {} -n {}".format(gene, email, name) + callset_arg + cds_arg + aa_arg
    results = subprocess.run(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, encoding="utf-8")
    logger.info("The searching process goes like this: \n{}\n*************************This is the END of the log file***************************\n".format(results.stdout))
    url_lines = [ l for l in results.stdout.split("\n") if re.search(r"^The view URL of the uploaded file is", l) ]
    if len(url_lines) > 0:
        view_URL = url_lines[-1].strip().split(" ")[-1]
    else:
        view_URL = np.nan
        logger.warning(f"Cant find a matching record for variant, {chromosome}:{start}-{end}:{ref_allele}->{alt_allele}")
        return np.nan
    logger.info("The fetched view URL is: {}".format(view_URL))
    if len(view_URL) > 40:
        # Now the file is ready, we can share it with the requester
        share_google_drive_file(file_id = view_URL.strip().split("/")[-2], email_address = email)
        return view_URL
    else:
        return np.nan


# Please give me a python implementation of a function that is designed for these purposes:
# 1. The function needs to have access to a google drive file that authorized by service account file service_account_file = "/paedyl01/disk1/google_cloud/hkupid-352116-2ceb75fa01b4.json",
# 2. The function needs to share that google drive file with a specified email address
# 3. The function uses google api to share the file with the email address
# Below is the function.
def share_google_drive_file(file_id = "1LZbX9XrQ1Q2JHv0x8yW3tY0wqY3TtjZ8",
                            service_account_json = "/paedyl01/disk1/google_cloud/hkupid-352116-2ceb75fa01b4.json",
                            email_address = ""):
    import os
    import json
    import googleapiclient
    from google.oauth2.service_account import Credentials
    from googleapiclient.discovery import build

    # Set the path to the service account JSON file
    json_path = service_account_json

    # Set the ID of the file you want to share

    # Set the email address of the user you want to share the file with
    email = email_address

    # Load the service account JSON file and authorize the credentials
    with open(json_path) as json_file:
        json_data = json.load(json_file)
    creds = Credentials.from_service_account_info(json_data)

    # Create a Drive API client
    service = build('drive', 'v3', credentials=creds)

    # Add the user to the file's list of permissions
    permission = {
        'type': 'user',
        'role': 'writer',
        'emailAddress': email
    }
    
    try:
        res = service.permissions().create(fileId=file_id, body=permission).execute()
    except googleapiclient.errors.HttpError as ge:
        logger.error("This file runs into a problem trying to share it with others: {}".format(file_id))
    else:
        logger.info("Succesfully share the file {} with user {}".format(file_id, email_address))

    # Print the result
    # print(res)


    
def respond_to_search_requests(sheet_id = "12cIhykbPOR0ynWyAJNcTc0RowfMvB1MQqa4h0Ia6K3w",
                               progress_sheet_name = "Query_progress_track",
                               threads = 6):
    from datetime import datetime
    import multiprocessing as mp
    work_queries, ready_query_df = track_search_request_progress(sheet_id)
    if len(work_queries) == 0:
        logger.info("No new queries or unfinished queries at this moment: {}".format(datetime.now()))
        return
    
    # Else start to process the work_queries
    df = get_latest_search_request(sheet_id)
    query_lines = df.loc[df["query_ID"].isin(work_queries), :]
    logger.info("{}\n".format(query_lines.to_string(index=False)))
    
    # Need to fetch the ID URL mapping after the search process
    pool = mp.Pool(threads)
    query_lines["result_URL"] = pool.map(search_per_request, work_queries)
    logger.info("{}\n".format(query_lines.to_string(index=False)))
    id_url_map = query_lines.loc[:, ["query_ID","result_URL"]].drop_duplicates().dropna(subset=["result_URL"])
    
    logger.info("The fetched results sheet looks like this: \n{}".format(id_url_map.to_string(index=False)))
    
    # Need to update the progress sheet
    progress_df = pd.concat([ready_query_df, id_url_map], axis = 0, ignore_index=True)
    logger.info("{}\n".format(progress_df.to_string(index=False)))
    progress_ws = get_worksheet_object(sheet_id, sheet_name = progress_sheet_name)
    progress_ws.set_dataframe(progress_df, start="A1", nan="", copy_index=False, copy_head=True)
    # logger.info("{}\n".format(progress_ws.to_string(index=False)))

    
        
    
def update_the_responses_clean(**kwargs):
    from fetch_updated_form_response import UPDATED_RESPONSES
    UPDATED_RESPONSES(**kwargs).update_online_clean_sheet()
    


def download_folder_contents(folder_id, output_directory, service_account_json, export_format=None):
    from google.oauth2 import service_account
    from googleapiclient.discovery import build
    def list_files_in_folder(service, folder_id):
        files = []
        query = f"'{folder_id}' in parents"
        page_token = None
        while True:
            response = service.files().list(q=query,
                                            spaces='drive',
                                            fields='nextPageToken, files(id, name, mimeType)',
                                            pageToken=page_token).execute()
            files.extend(response.get('files', []))
            page_token = response.get('nextPageToken', None)
            if page_token is None:
                break
        return files

    # Load service account JSON and create credentials
    with open(service_account_json, 'r') as f:
        sa_info = json.load(f)
    credentials = service_account.Credentials.from_service_account_info(sa_info, scopes=['https://www.googleapis.com/auth/drive'])

    # Build the API client
    service = build('drive', 'v3', credentials=credentials)

    # Create the output directory if it does not exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # List all files and folders in the given folder ID
    contents = list_files_in_folder(service, folder_id)

    for item in contents:
        item_id = item['id']
        item_name = item['name']
        item_mime_type = item['mimeType']

        # Download files and recursively process subfolders
        if item_mime_type == 'application/vnd.google-apps.folder':
            new_output_directory = os.path.join(output_directory, item_name)
            download_folder_contents(item_id, new_output_directory, service_account_json, export_format)
        else:
            output_file = os.path.join(output_directory, item_name)
            download_google_drive_file(item_id, output_file, service_account_json, export_format)



def download_google_drive_file(file_id, output_file, service_account_json, export_format=None):
    from google.oauth2 import service_account
    from googleapiclient.discovery import build
    from googleapiclient.errors import HttpError
    from googleapiclient.http import MediaIoBaseDownload

    try:
        # Load service account JSON and create credentials
        with open(service_account_json, 'r') as f:
            sa_info = json.load(f)
        credentials = service_account.Credentials.from_service_account_info(sa_info, scopes=['https://www.googleapis.com/auth/drive'])

        # Build the API client
        service = build('drive', 'v3', credentials=credentials)

        # Get file metadata
        file_metadata = service.files().get(fileId=file_id).execute()
        file_size = file_metadata.get('size', None)
        file_name = file_metadata.get('name', None)
        mime_type = file_metadata.get('mimeType', None)

        logger.info(f"Downloading file '{file_name}' with ID '{file_id}'...")

        if "google-apps" in mime_type:
            # Set the export MIME type based on the provided export_format
            if mime_type == "application/vnd.google-apps.document":
                if export_format == 'pdf':
                    export_mime_type = 'application/pdf'
                elif export_format == 'docx':
                    export_mime_type = 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
                elif export_format == 'txt':
                    export_mime_type = 'text/plain'
                else:
                    raise ValueError(f"Invalid export format '{export_format}' for Google Docs")
            elif mime_type == "application/vnd.google-apps.spreadsheet":
                if export_format == 'xlsx':
                    export_mime_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
                elif export_format == 'csv':
                    export_mime_type = 'text/csv'
                else:
                    raise ValueError(f"Invalid export format '{export_format}' for Google Sheets")
            elif mime_type == "application/vnd.google-apps.presentation":
                if export_format == 'pptx':
                    export_mime_type = 'application/vnd.openxmlformats-officedocument.presentationml.presentation'
                elif export_format == 'pdf':
                    export_mime_type = 'application/pdf'
                else:
                    raise ValueError(f"Invalid export format '{export_format}' for Google Slides")

            request = service.files().export_media(fileId=file_id, mimeType=export_mime_type)
            response = request.execute()
            with open(output_file, 'wb') as f:
                f.write(response)

        else:
            # If file size is larger than 10MB, use media download
            if file_size and int(file_size) > 10 * 1024 * 1024:
                request = service.files().get_media(fileId=file_id)
                with open(output_file, 'wb') as f:
                    downloader = MediaIoBaseDownload(f, request)
                    done = False
                    while done is False:
                        status, done = downloader.next_chunk()
                        logger.info(f"Downloaded {int(status.progress() * 100)}.")
            else:
                # Download the file
                request = service.files().get_media(fileId=file_id)
                response = request.execute()
                with open(output_file, 'wb') as f:
                    f.write(response)

        logger.info(f"File '{file_name}' has been downloaded to '{output_file}'.")

    except HttpError as error:
        logger.info(f"An error occurred: {error}")
        return None


def query_the_internal_AF(**kwargs):
    from filter_var_by_inhouse_AF import main_process
    main_process(**kwargs)


def find_latest_modified_file(files, delete_old_files=False):
    from datetime import datetime
    files = files.split(",")
    
    latest_file = files[0]
    latest_mod_time = os.path.getmtime(latest_file)
    
    for file in files:
        mod_time = os.path.getmtime(file)
        if mod_time > latest_mod_time:
            if delete_old_files:
                os.remove(latest_file)
            latest_file = file
            latest_mod_time = mod_time
        elif delete_old_files:
            os.remove(file)
            
    logger.info("Latest modified file: " + latest_file)
    return latest_file


def draw_chord_diagram(SD_table_path, cytoband_file_path):
    import pandas as pd
    import numpy as np
    import holoviews as hv
    from holoviews import opts, dim
    hv.extension('bokeh')
 
    # Load cytoband file
    cytoband_df = pd.read_table(cytoband_file_path, names=['chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'])
    cytoband_df['chrom'] = cytoband_df['chrom'].map(lambda x: x.lstrip('chr')) # remove 'chr' prefix
    
    # Load SD_table
    SD_table_df = pd.read_table(SD_table_path, header=None, sep='\s+')
    
    def find_cyto_band(chr, start, end, cyto_data):
        mask = (cyto_data['chrom'] == chr) & (cyto_data['chromStart'] <= start) & (cyto_data['chromEnd'] >= end)
        bands = cyto_data.loc[mask, 'name']
        if len(bands) > 0:
            return bands.iloc[0]
        else:
            return None
    
    # create bins
    bins = []
    for index, row in SD_table_df.iterrows():
        bins.append(find_cyto_band(row[0], row[1], row[2], cytoband_df))
        bins.append(find_cyto_band(row[3], row[4], row[5], cytoband_df))
    
    # remove None values
    bins = [b for b in bins if b is not None]
    
    # create connectivity data frame
    data = []
    for index, row in SD_table_df.iterrows():
        band1 = find_cyto_band(row[0], row[1], row[2], cytoband_df)
        band2 = find_cyto_band(row[3], row[4], row[5], cytoband_df)
        if band1 is not None and band2 is not None:
            data.append([band1, band2, row[12]]) # use size for edge thickness, you can change this to gap rate or others
    
    df = pd.DataFrame(data, columns=['source', 'target', 'value'])
    
    chord = hv.Chord((df, pd.DataFrame(bins, columns=['name'])))
    chord.opts(
        opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), 
                   labels='name', node_color=dim('index').str()))
    
    return chord


def convert_gff_anno_to_bed(input_gff, output_bed):
    # The updated gff annotation download link is: 
    # https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
    if re.search(r"\.db$", input_gff):
        import gffutils
        db = gffutils.FeatureDB(input_gff)
    elif re.search(r"\.gff[\.gz]*$", input_gff):
        db = create_gffutils_db(input_gff)
    # Extract mRNA transcript IDs of all the CDS regions
    all_feature_ids = [ rec.attributes.get("ID", "")[0] for rec in db.all_features() if re.search(r"pseudogene|exon", rec.featuretype) and re.search(r"^NC_", rec.chrom) ]
    all_feature_ids = set(list(dict.fromkeys([x for x in all_feature_ids if len(x) > 0 ])))
    logger.info(f"How many IDs there are for the selected features: {len(all_feature_ids)}")
    # use the transcript IDs to extract the exon regions (exon includes UTR and CDS)
    all_target_features = [rec for rec in db.all_features() if re.search(r"^NC_", rec.chrom) and rec.attributes.get("ID","none-none-none")[0] in all_feature_ids ]
    with open(output_bed, "w") as ob:
        for rec in all_target_features:
            try:
                chr_ind = re.search(r"^NC_0+([1-9]+[0-9]*)\.[0-9]+$", rec.chrom).group(1)
            except AttributeError:
                logger.warning("This record cannot match the chromosome format: {}\n".format(rec))
                raise ValueError
            sex_chr_map = {"23": "X", "24": "Y", "12920": "M"}
            chr_ind = chr_ind if chr_ind not in ["23", "24", "12920"] else sex_chr_map[chr_ind]
            chrom = "chr" + chr_ind
            start = rec.start
            end = rec.end
            strand = rec.strand
            feature_id = rec.attributes.get("ID","none-none-none")[0].split("-")[1]
            ob.write(f"{chrom}\t{start}\t{end}\t{feature_id}\t{strand}\n")


def round_up(input_value):
    import math
    rounded = math.ceil(float(input_value))
    print(rounded)
    return rounded
    

def modify_sequence_per_contig(input_tuple):
    chrom, seq, regions = input_tuple
    seq = str(seq.seq)
    modified_seq = ""
    logger.info("This is the sequence of chromosome {} with length {}".format(chrom, len(seq)))
    if regions:
        # Sort the regions by start position
        regions = sorted(regions, key=lambda x: x[0])
        for i in range(0, len(regions)):
            region = regions[i]
            start, end = region
            last_region = regions[i-1] if i > 0 else []
            next_region = regions[i+1] if i < len(regions) - 1 else []
            if len(regions) == 1:
                logger.info("This is the only region of chromosome {} with start {} and end {}".format(chrom, start, end))
                modified_seq = seq[:start] + seq[start:end].upper() + seq[end:]
            elif i == 0:
                # The first interval (most upstream)
                logger.info("This is the first region of chromosome {} with start {} and end {}".format(chrom, start, end))
                modified_seq = seq[:start] + seq[start:end].upper()
            elif i == len(regions) - 1:
                # The last interval
                logger.info("This is the last region of chromosome {} with start {} and end {}".format(chrom, start, end))
                last_start, last_end = last_region
                modified_seq = modified_seq + seq[last_end:start] + seq[start:end].upper() + seq[end:]
            else:
                # Middle intervals (most downstream)
                # logger.info("This is the middle region of chromosome {} with start {} and end {}".format(chrom, start, end))
                last_start, last_end = last_region
                modified_seq = modified_seq + seq[last_end:start] + seq[start:end].upper()
        return (chrom, modified_seq, True)
    else:
        modified_seq = seq
        return (chrom, modified_seq, False)


def lift_softmask_target_region(unwrap_fasta="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.nowrap.fasta",
                                target_region="/paedyl01/disk1/yangyxt/public_data/gene_annotation/all_coding_region.merged.pad150.bed",
                                output_fasta="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.coding_nomask.fasta",
                                threads=24):
    import pybedtools as pb
    sorted_region_df = pb.BedTool(target_region).sort().merge().to_dataframe(disable_auto_names = True,
                                                                          names=["chr", "start", "end"],
                                                                          usecols = [0,1,2])
    target_regions = list(zip(sorted_region_df["chr"], sorted_region_df["start"].astype(int), sorted_region_df["end"].astype(int)))
    logger.info("How many intervals are there in the target region: {}".format(len(target_regions)))
    all_contigs = list(dict.fromkeys([l[0] for l in target_regions]))
    target_region_by_contig = {c:[] for c in all_contigs}
    for l in target_regions:
        target_region_by_contig[l[0]].append(tuple([int(l[1]), int(l[2])]))
    logger.info("How many intervals are there in the target region of chr1: {}".format(len(target_region_by_contig["chr1"])))

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import multiprocessing as mp

    pool = mp.Pool(threads)
    sequences = SeqIO.to_dict(SeqIO.parse(unwrap_fasta, "fasta"))
    logger.info("What is the type of sequences: {}".format(type(sequences)))
    logger.info("What is the type of sequences of chr1: {}".format(type(sequences["chr1"])))
    input_tuples = [(chrom, seq, tuple(target_region_by_contig.get(chrom, tuple([])))) for chrom, seq in sequences.items()]
    results = pool.imap_unordered(modify_sequence_per_contig, input_tuples)
    modified_sequences = { c: "" for c in all_contigs}
    for chrom, seq, modified in results:
        assert str(seq).upper() == str(sequences[chrom].seq).upper()
        if modified:
            assert seq != str(sequences[chrom].seq)
        else:
            logger.warning("Seems like the sequence for chromosome {} is not modified because its not in the target region bed file {}".format(chrom, target_region))
        seq_record = SeqRecord(Seq(seq), id=chrom, description="")
        modified_sequences[chrom] = seq_record

    pool.close()

    with open(output_fasta, "w") as op:
        SeqIO.write(modified_sequences.values(), op, "fasta")
    


def calculate_inferred_coverage(bam_file, 
                              output_tsv = "",
                              return_df = False,
                              min_mapq = 10,
                              filter_tags = [],
                              filter_logic = "and",
                              genome_file = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.contigsize.genome",
                              target_region = "",
                              target_tag = "",
                              logger=logger):
    import pysam
    import pybedtools

    if type(filter_tags) == str:
        filter_tags = filter_tags.split(",")

    if target_region and not target_tag:
        target_tag = os.path.basename(target_region).split(".")[0]

    if not target_region:
        target_tag = None
    
    if len(filter_tags) > 0:
        name_filter_tag = "_".join([filter_logic] + filter_tags)
    else:
        name_filter_tag = ""

    if len(output_tsv) == 0:
        if name_filter_tag:
            output_tsv = bam_file.replace(".bam", f".infercov.minMQ_{min_mapq}.{name_filter_tag}.tsv")
        else:
            output_tsv = bam_file.replace(".bam", f".infercov.minMQ_{min_mapq}.tsv")
        if target_tag:
            output_tsv = output_tsv.replace(".tsv", f".{target_tag}.tsv")

    output_bed = output_tsv.replace(".tsv", ".bed")
    tmp_output_bed = output_bed.replace(".bed", f".{uuid.uuid4()}.bed")

    skip = False
    if os.path.exists(output_bed):
        if os.path.getmtime(output_bed) > os.path.getmtime(bam_file):
            logger.info("The output bed file {} is newer than the input bam file {}, so we will skip the coverage calculation step".format(output_bed, bam_file))
            skip = True
            # Sometimes the adjacent interval might fetch the same pair twice. To prevent double counting in the depth calculation. We need to drop the duplicates

    if not skip:
        if target_region:
            bed_obj = pybedtools.BedTool(target_region).sort().merge()
        else:
            bed_obj = None

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            with open(tmp_output_bed, "w") as ob:
                bam_iters = []
                if bed_obj:
                    for interval in bed_obj:
                        bam_iter = bam.fetch(interval.chrom, interval.start, interval.end)
                        for read in bam_iter:
                            if filter_tags:
                                if filter_logic == "or":
                                    filter_pass = any([read.has_tag(tag) for tag in filter_tags])
                                elif filter_logic == "and":
                                    filter_pass = all([read.has_tag(tag) for tag in filter_tags])
                                elif filter_logic == "not":
                                    filter_pass = not any([read.has_tag(tag) for tag in filter_tags])
                                else:
                                    logger.warning("The filter logic {} is not supported, please use either 'or' or 'and'. Now just disable the filters".format(filter_logic))
                                    filter_pass = True
                            else:
                                filter_pass = True

                            if (not read.is_duplicate) and \
                                (not read.is_unmapped) and \
                                (read.mapping_quality >= min_mapq) and \
                                (not read.is_secondary) and \
                                (not read.is_supplementary) and \
                                (not read.is_qcfail) and \
                                (not read.mate_is_unmapped) and \
                                read.is_proper_pair and filter_pass:
                                # Ensure we process each pair only once
                                if read.is_read1:
                                    chrom = read.reference_name
                                    start = read.reference_start
                                    end = read.reference_end
                                    read_query_name = read.query_name

                                    if read.next_reference_name == chrom:
                                        ob.write(f"{chrom}\t{min(start, read.next_reference_start)}\t{max(end, read.next_reference_start + 150)}\t{read_query_name}\n")
                                    else:
                                        ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")
                            elif (not read.is_duplicate) and \
                                (not read.is_unmapped) and \
                                (read.mapping_quality >= min_mapq) and \
                                (not read.is_secondary) and \
                                (not read.is_supplementary) and \
                                (not read.is_qcfail) and \
                                (not read.is_proper_pair) and \
                                read.is_paired and filter_pass:
                                    chrom = read.reference_name
                                    start = read.reference_start
                                    end = read.reference_end
                                    read_query_name = read.query_name
                                    ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")
                            elif (not read.is_duplicate) and \
                                (not read.is_unmapped) and \
                                (read.mapping_quality >= min_mapq) and \
                                (not read.is_secondary) and \
                                (not read.is_supplementary) and \
                                (not read.is_qcfail) and \
                                (not read.is_paired) and filter_pass:
                                    chrom = read.reference_name
                                    start = read.reference_start
                                    end = read.reference_end
                                    read_query_name = read.query_name
                                    ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")
                else:
                    bam_iter = bam.fetch()        
                    for read in bam_iter:
                        if filter_tags:
                            if filter_logic == "or":
                                filter_pass = any([read.has_tag(tag) for tag in filter_tags])
                            elif filter_logic == "and":
                                filter_pass = all([read.has_tag(tag) for tag in filter_tags])
                            elif filter_logic == "not":
                                filter_pass = not any([read.has_tag(tag) for tag in filter_tags])
                            else:
                                logger.warning("The filter logic {} is not supported, please use either 'or' or 'and'. Now just disable the filters".format(filter_logic))
                                filter_pass = True
                        else:
                            filter_pass = True

                        if (not read.is_duplicate) and \
                            (not read.is_unmapped) and \
                            (read.mapping_quality >= min_mapq) and \
                            (not read.is_secondary) and \
                            (not read.is_supplementary) and \
                            (not read.is_qcfail) and \
                            (not read.mate_is_unmapped) and \
                            read.is_proper_pair and filter_pass:
                            # Ensure we process each pair only once
                            if read.is_read1:
                                chrom = read.reference_name
                                start = read.reference_start
                                end = read.reference_end
                                read_query_name = read.query_name

                                if read.next_reference_name == chrom:
                                    ob.write(f"{chrom}\t{min(start, read.next_reference_start)}\t{max(end, read.next_reference_start + 150)}\t{read_query_name}\n")
                                else:
                                    ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")
                        elif (not read.is_duplicate) and \
                            (not read.is_unmapped) and \
                            (read.mapping_quality >= min_mapq) and \
                            (not read.is_secondary) and \
                            (not read.is_supplementary) and \
                            (not read.is_qcfail) and \
                            (not read.is_proper_pair) and \
                            read.is_paired and filter_pass:
                                chrom = read.reference_name
                                start = read.reference_start
                                end = read.reference_end
                                read_query_name = read.query_name
                                ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")
                        elif (not read.is_duplicate) and \
                            (not read.is_unmapped) and \
                            (read.mapping_quality >= min_mapq) and \
                            (not read.is_secondary) and \
                            (not read.is_supplementary) and \
                            (not read.is_qcfail) and \
                            (not read.is_paired) and filter_pass:
                                chrom = read.reference_name
                                start = read.reference_start
                                end = read.reference_end
                                read_query_name = read.query_name
                                ob.write(f"{chrom}\t{start}\t{end}\t{read_query_name}\n")

        # Sometimes the adjacent interval might fetch the same pair twice. To prevent double counting in the depth calculation. We need to drop the duplicates
        bed_df = pd.read_table(tmp_output_bed, header=None, names=['chrom', 'start', 'end', 'read_query_name'], low_memory=False).dropna()
        bed_df.loc[:, "start"] = bed_df.loc[:, "start"].astype(int)
        bed_df.loc[:, "end"] = bed_df.loc[:, "end"].astype(int)
        logger.info("The bed file looks like this: \n{}".format(bed_df[:5].to_string(index=False)))
        bed_df.to_csv(tmp_output_bed, sep="\t", header=False, index=False)
        executeCmd(f"sort -k1,1 -k2,2n {tmp_output_bed} > {tmp_output_bed}.tmp && mv {tmp_output_bed}.tmp {output_bed} && rm {tmp_output_bed}")
    
    tmp_output_tsv = output_tsv.replace(".tsv", f".{uuid.uuid4()}.tsv")

    cmd = f"bedtools genomecov -bg -i {output_bed} -g {genome_file} | \
            sort -k 1,1 -k 2,2n - | \
            mawk -F '\\t' '{{ for (i=$2+1;i<=$3;i++) {{printf \"%s\\t%i\\t%i\\n\", $1, i, $4;}} }}' - > {tmp_output_tsv} && \
            bgzip -f -c {tmp_output_tsv} > {tmp_output_tsv}.gz && \
            tabix -s 1 -b 2 -e 2 {tmp_output_tsv}.gz"

    skip = False
    if os.path.exists(output_tsv):
        if os.path.getmtime(output_tsv) > os.path.getmtime(output_bed):
            skip = True
            if os.path.exists(output_tsv + ".gz"):
                if os.path.getmtime(output_tsv + ".gz") < os.path.getmtime(output_tsv):
                    cmd = f"bgzip -f -c {output_tsv} > {output_tsv}.gz && tabix -s 1 -b 2 -e 2 {output_tsv}.gz"
                    executeCmd(cmd)
                elif os.path.exists(output_tsv + ".gz.tbi"):
                    if os.path.getmtime(output_tsv + ".gz.tbi") < os.path.getmtime(output_tsv + ".gz"):
                        cmd = f"tabix -s 1 -b 2 -e 2 {output_tsv}.gz"
                        executeCmd(cmd)
                else:
                    cmd = f"tabix -s 1 -b 2 -e 2 {output_tsv}.gz"
                    executeCmd(cmd)
            else:
                cmd = f"bgzip -f -c {output_tsv} > {output_tsv}.gz && tabix -s 1 -b 2 -e 2 {output_tsv}.gz"
                executeCmd(cmd)
    if skip:
        logger.info("The output tsv file {} is newer than the input bed file {}, so we will skip the coverage calculation step".format(output_tsv, output_bed))
    else:
        executeCmd(cmd)
        executeCmd(f"mv {tmp_output_tsv} {output_tsv} && \
                     mv {tmp_output_tsv}.gz {output_tsv}.gz && \
                     mv {tmp_output_tsv}.gz.tbi {output_tsv}.gz.tbi")
    
    coverage_df = pd.read_table(output_tsv, header=None, names=['chrom', 'pos', 'depth'], low_memory=False)
    logger.info("The coverage table {} looks like: \n{}".format(output_tsv, coverage_df[:5].to_string(index=False)))

    if return_df:
        return coverage_df
    else:
        return output_tsv + ".gz"



def retrieve_infer_depth_bam(bam_file, chrom, pos, min_mapq = 10, logger=logger, **kwargs):
    # Retrieve the depth info at single-base resolution
    depth_tab_gz = calculate_inferred_coverage(bam_file, return_df = False, min_mapq = min_mapq, **kwargs)
    import pysam
    with pysam.TabixFile(depth_tab_gz) as tabix_file:
        records = list(tabix_file.fetch(chrom, int(pos) - 1, int(pos)))
        if len(records) > 0:
            return int(records[0].split('\t')[2])
        else:
            return 0


def print_element_info(element, depth=0, elements_info={}):
    """
    Recursive function to print element information with indentation.
    """
    # Print the current element with indentation
    indent = "  " * depth  # 4 spaces for each level of depth
    attributes = ', '.join(elements_info[element]['attributes'])
    print(f"{indent}Element: {element} (Attributes: {attributes})")

    # Recursively print child elements
    for child in elements_info[element]['children']:
        print_element_info(child, depth + 1, elements_info)


def convert_all_sets_to_lists(d):
    """
    Recursively convert all sets in the dictionary to lists.
    """
    if isinstance(d, set):
        return list(d)
    elif isinstance(d, dict):
        return {k: convert_all_sets_to_lists(v) for k, v in d.items()}
    elif isinstance(d, list):
        return d
    return d


def defaultdict_to_dict(d):
    from collections import defaultdict
    for key, value in d.items():
        if isinstance(value, defaultdict):
            d[key] = defaultdict_to_dict(value)
        elif isinstance(value, set):
            d[key] = list(value)
        elif isinstance(value, dict):
            d[key] = convert_all_sets_to_lists(value)
    return dict(d)


def dict_to_defaultdict(d):
    from collections import defaultdict
    """
    Recursively convert a dictionary to a defaultdict.
    Also converts lists back to sets for 'attributes' and 'children'.
    """
    if isinstance(d, dict):
        new_dict = defaultdict(lambda: {"attributes": set(), "children": set(), "parent": set()})
        for key, value in d.items():
            if key in ['attributes', 'children', "parent"]:
                # Convert lists back to sets for 'attributes' and 'children'
                new_dict[key] = set(value)
            else:
                # Recursively convert other dictionaries
                new_dict[key] = dict_to_defaultdict(value)
        return new_dict
    return d  # Return the item as is if it's not a dictionary


def summary_xml_element_map(xml_file):
    from lxml import etree
    from collections import defaultdict
    
    json_file = xml_file.replace(".xml", ".json")

    if not os.path.exists(xml_file):
        raise ValueError(f"The input XML file '{xml_file}' does not exist")
    
    if os.path.exists(json_file):
        if os.path.getmtime(json_file) > os.path.getmtime(xml_file):
            with open(json_file, 'r') as jf:
                elements_info_loaded = json.load(jf)
                elements_info = dict_to_defaultdict(elements_info_loaded)
    else:
        # Initialize data structures
        elements_info = defaultdict(lambda: {"attributes": set(), "children": set(), "parent": set()})

        # Parse the XML file
        for event, elem in etree.iterparse(xml_file, events=('start', 'end'), recover=True):
            # Update element information
            if event == 'start':
                # Add attributes of the element
                elements_info[elem.tag]["attributes"].update(elem.attrib.keys())

                # Add parent-child relationship
                parent = elem.getparent()
                if parent is not None:
                    elements_info[elem.tag]["parent"].add(parent.tag)
                    elements_info[parent.tag]["children"].add(elem.tag)

            # Memory management
            if event == 'end':
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]

        # Convert elements_info to a standard dictionary
        elements_info_dict = defaultdict_to_dict(elements_info)
        # Serialize to JSON and write to a file
        
        with open(json_file, 'w') as json_file:
            json.dump(elements_info_dict, json_file, indent=4)

    # Display the collected information
    # Find root elements (elements with no parent)
    root_elements = [elem for elem, info in elements_info.items() if len(info['parent']) == 0]

    print("The top-level element is: {}".format(root_elements))
    assert all([len(elements_info[e]["parent"]) != 0 for e in elements_info.keys() if e not in root_elements])
    print(f"The total structure are:\n\n")

    # Perform depth-first traversal from each root element
    for root in root_elements:
        print_element_info(root, elements_info=elements_info)


def get_HC_patho_clinvar_tab(output_table):
    from prepare_updated_clinvar_recs import UPDATED_CLINVAR
    import pandas as pd

    clinvar_cls = UPDATED_CLINVAR()
    df = clinvar_cls.convert_var_archives_to_tab()
    logger.info("The clinvar table looks like this: \n{}".format(df[:5].to_string(index=False)))
    df.to_csv(output_table, sep="\t", index=False)
    logger.info("The clinvar table is saved to: {}".format(output_table))



def find_initiation_codon_usage(output_prefix = "/paedyl01/disk1/yangyxt/public_data/gene_annotation/ensembl_geneids",
                                ensembl_tranx_cds_map = "/paedyl01/disk1/yangyxt/public_data/ensembl_resources/manual/GRCh37_tranx_cds_len.tsv"):
    from biomart import BiomartServer
    server = BiomartServer("http://www.ensembl.org/biomart")
    dataset = server.datasets['hsapiens_gene_ensembl']
    
    # Determine which transcripts contain exons
    tranx_cds_map = pd.read_table(ensembl_tranx_cds_map).set_index("TranxID").loc[:, "CDS_len"].to_dict()

    response = dataset.search({
        'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'transcript_start'],
        'filters': {'biotype': 'protein_coding'},
    })

    gene_transcripts = {}
    same_initiation_codon = set()
    different_initiation_codon = set()

    for row in response.iter_lines():
        row = row.decode('utf-8')
        ensembl_geneid, transcript_id, transcript_start = row.split("\t")
        
        if ensembl_geneid not in gene_transcripts:
            gene_transcripts[ensembl_geneid] = {"coding_transcripts":[], "noncoding_transcripts":[], "start_codon":[]}
            
        if float(tranx_cds_map.get(transcript_id, 0)) > 0:
            logger.info("Transcript {} has at least one CDS region".format(transcript_id))
            gene_transcripts[ensembl_geneid]["coding_transcripts"].append(transcript_id)
        else:
            logger.info("Transcript {} does not have CDS region".format(transcript_id))
            gene_transcripts[ensembl_geneid]["noncoding_transcripts"].append(transcript_id)
            continue
        
        if transcript_start not in gene_transcripts[ensembl_geneid]["start_codon"]:
            gene_transcripts[ensembl_geneid]["start_codon"] = []
        
        gene_transcripts[ensembl_geneid]["start_codon"].append(transcript_start)

    for ensembl_geneid, maps in gene_transcripts.items():
        if len(set(maps["start_codon"])) == 1:
            same_initiation_codon.add(ensembl_geneid)
        else:
            different_initiation_codon.add(ensembl_geneid)
            
    with open(output_prefix + ".same_init_codon.txt", "w") as sf:
        for g in same_initiation_codon: sf.write(g + "\n")
    logger.info("The output file contain a list of ensembl gene IDs that all transcripts are using the same init codon: {}.".format(output_prefix + ".same_init_codon.txt"))    
        
        
    with open(output_prefix + ".diff_init_codon.txt", "w") as df:
        for g in different_initiation_codon: df.write(g + "\n")
    logger.info("The output file contain a list of ensembl gene IDs that all transcripts are using at least two different init codons: {}.".format(output_prefix + ".diff_init_codon.txt"))

    return same_initiation_codon, different_initiation_codon



def add_header_line(header, field, ID, number='', type_str='', description=""):
    import vcfpy
    from collections import OrderedDict
    '''
    Add header line if original header does not contain the required meta line
    '''
    logger.info("Before adding header line of ID {} in field {}, now the header has {}-ID pairs like this: \n\t{}".format(ID, field, field, "\n\t".join([ "\n\t".join([ str(k) + ":" + str(v) for k, v in d.items() ]) for key,d in header._indices.items() if key == field ])))
    if header.has_header_line(field, ID):
        logger.warning("This field {} in {} is already in the header. We want to replace the existing headerlines with this one.".format(ID, field))
        try:
            assert field in ["INFO", "FORMAT", "contig", "ALT", "FILTER"]
            original_line = [line for line in header.lines if line.key == field and line.id == ID][0]
            new_header = vcfpy.Header([line.copy() for line in header.lines if line.key != field or line.id != ID], header.samples.copy())
        except AttributeError:
            original_line = [line for line in header.lines if line.key == field and line._value == ID][0]
            new_header = vcfpy.Header([line.copy() for line in header.lines if line.key != field or line._value != ID], header.samples.copy())
        header = new_header
        description = original_line.description if len(description) == 0 else description
        try:
            type_str = original_line.type if len(type_str) == 0 else type_str
        except AttributeError:
            type_str = "String"
        number = "1" if len(number) == 0 else number
    else:
        number = "1" if len(number) == 0 else number
        type_str = "String" if len(type_str) == 0 else type_str
        description = "ID {} for field {}".format(ID, field) if len(description) == 0 else description
    
    logger.info("The new header line for this field {} in {} will have number {}, type {}, description: {}.".format(ID, field, number, type_str, description))

    ##Add alt allele to header
    if field == "ALT":
        DEL = vcfpy.AltAlleleHeaderLine('ALT','<ID=' + ID + ',Description="' + description + '">',{'ID':ID, 'Description':description})
        header.add_line(DEL)
    ##Add contig line
    elif field == "contig":
        mapping=OrderedDict({'ID':ID,'length':description})
        header.add_contig_line(mapping)
    #Add filter line
    elif field == "FILTER":
        mapping=OrderedDict({'ID':ID,'Description':description})
        header.add_filter_line(mapping)
    #Add format line
    elif field == "FORMAT":
        mapping=OrderedDict({'ID':ID,'Number':number,'Type':type_str,'Description':description})
        header.add_format_line(mapping)
    ##Add info line
    elif field == "INFO":
        mapping=OrderedDict({'ID':ID, 'Number':number, 'Type':type_str, 'Description':description})
        header.add_info_line(mapping)
        
    logger.info("After adding header line, now the header has {}-ID pairs like this: \n\t{}".format(field, "\n\t".join([ "\n\t".join([ str(k) + ":" + str(v) for k, v in d.items() ]) for key,d in header._indices.items() if key == field ])))
    return header



class HOMOSEQ_REGION:
    def __init__(self, chrom, start, end, strand, ref_fasta=""):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.ref_fasta = ref_fasta
        self.rela_start = 0
        self.rela_end = end - start
        self.data = tuple([self.chrom, self.start, self.end, self.strand])

    def __hash__(self):
        return hash(self.data)

    def __repr__(self):
        return ", ".join(str(x) for x in [self.chrom, self.start, self.end, self.strand])

    def __eq__(self, other):
        return self.chrom == other[0] and self.start == other[1] and self.end == other[2] and self.strand == other[3]

    def __iter__(self):
        for item in [self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand]:
            yield item

    def __getitem__(self, index):
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self.data))
            return [self.data[i] for i in range(start, stop, step)]
        else:
            return self.data[index]

    def __len__(self):
        return len(self.data)

    def fix_coord(self):
        return tuple([self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand])

        # Now that we should only locate one sample record
# def read_google_sheets_to_df(sheet_id="184Ne4LoZdJida-hnXVLVoDpjd_L9I7U8rmb00E5-Fts", header=0):
#     '''
#     To use this function, you need to first install gdrive
#     and configure the OUAth for your google drive account
#     '''
#     import subprocess
#     import re
#     import sys
#     import os
#     cmd = "gdrive info {}".format(sheet_id)
#     result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
#     name = [l for l in result.stdout.split("\n") if re.search(r"Name:", l)][0].strip("\n ").split(": ")[-1]
#     cmd = "gdrive export -f --mime text/tab-separated-values {}".format(sheet_id)
#     print(subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").stdout, file=sys.stderr)
#     return pd.read_table(os.path.join(os.getcwd(), name + ".tsv"), low_memory=False, header=header)

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--function", type=str, help="The function name", required=True)
    parser.add_argument("-a", "--arguments", type=str, help="The function's input arguments, delimited by semi-colon ;", required=False, default=None)
    parser.add_argument("-k", "--key_arguments", type=str, help="Keyword arguments for the function, delimited by semi-colon ;", required=False, default=None)
    
    args = parser.parse_args()
    try:
        fargs = [ convert_input_value(a) for a in args.arguments.split(";") ] if type(args.arguments) == str else []
        fkwargs = { t.split("=")[0]: convert_input_value(t.split("=")[1]) for t in args.key_arguments.split(";") } if type(args.key_arguments) == str else {}
        logger.info("Running function: {}, input args are {}, input kwargs are {}".format(args.function, fargs, fkwargs))
    except Exception as e:
        logger.error("Input argument does not meet the expected format, encounter Parsing error {}, Let's check the input:\n-f {}, -a {}, -k {}".format(
            e,
            args.function,
            args.arguments,
            args.key_arguments
        ))
        raise e
    globals()[args.function](*fargs, **fkwargs)