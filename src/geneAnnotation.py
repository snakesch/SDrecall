#!/usr/bin/env python3
import pandas as pd
from pandarallel import pandarallel as pa
import argparse as ap
import logging

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s',
                    datefmt='%a %b-%m %I:%M:%S%P', level = "INFO")
logger = logging.getLogger("root")

pa.initialize(progress_bar=False, verbose=0)

def get_coor(row):
    '''Takes a row with all exon/intron coordinates and returns a dataframe with individual exon/intron details'''

    logger = logging.getLogger("root")

    if not "exonStarts" in row.index or not "exonEnds" in row.index:
        raise KeyError("exonStarts or exonEnds not found in fetched contents. ")

    if not "cdsStart" in row.index or not "cdsEnd" in row.index:
        raise KeyError("cdsStart or cdsEnd not found in fetched contents. ")

    exonStarts, exonEnds = list(map(int, row["exonStarts"].split(",")[:-1])), list(map(int, row["exonEnds"].split(",")[:-1]))
    intronEnds = list(map(lambda x: int(x) - 1, exonStarts[1:])) # 1 unit offset at the beginning
    intronStarts = list(map(lambda x: int(x) + 1, exonEnds[:-1])) # 1 unit offset at the end

    for start, end in zip(intronStarts, intronEnds):
        if start == end: # Case: Contiguous exons with no intron between
            intronStarts.remove(start)
            intronEnds.remove(end)

    isCoding = (int(row['cdsStart']) < int(row['cdsEnd']))

    if row['strand'] == '+':
        utr_5_starts = row['txStart']
        utr_5_ends = row['cdsStart'] - 1 if isCoding else row['txStart']
        utr_3_starts, utr_3_ends = row['cdsEnd'] + 1, row['txEnd']
    elif row["strand"] == "-":
        # Reverse coordinates
        exonStarts, exonEnds = exonStarts[::-1], exonEnds[::-1]
        intronStarts, intronEnds = intronStarts[::-1], intronEnds[::-1]
        utr_5_starts, utr_5_ends = row['cdsEnd'] + 1, row['txEnd']
        utr_3_starts = row['txStart']
        utr_3_ends = row['cdsStart'] - 1 if isCoding else row['txStart']
    else:
        raise ValueError(f"Unrecognized strand symbol {row['strand']} . ")

    intronStarts.insert(0, utr_5_starts)
    intronStarts.append(utr_3_starts)
    intronEnds.insert(0, utr_5_ends)
    intronEnds.append(utr_3_ends)

    logger.debug(f"Intron start: {intronStarts}; \nIntron ends: {intronEnds}")

    ### Create a dataframe of all exon and intron coordinates
    exon_dict = [{'Chr': row['chrom'], 'Start': start, 'End': end,
                 'Symbol': row['name2'], 'Feature': 'exon_' + str(idx+1),
                 'Strand': row['strand'], 'Interval_length': end-start,
                 'NCBI_ID': row['name']}
                 for idx, (start, end) in enumerate(zip(exonStarts, exonEnds))]
    ### Need to resolve 5-UTR and 3-UTR for introns
    intron_dict = []
    complexExonCount = 0
    for idx, (start, end) in enumerate(zip(intronStarts, intronEnds)):
        __dict = {'Chr': row['chrom'], 'Start': start, 'End': end,
                     'Symbol': row['name2'], 'Feature': 'intron_' + str(idx - complexExonCount),
                     'Strand': row['strand'], 'Interval_length': end-start,
                     'NCBI_ID': row['name']}
        # Note: Need to fix start/end positions for UTRs
        if idx == 0:
            __dict['Feature'] = '5-UTR'
            if row['strand'] == '-' and end - start == -1:
                __dict['End'] = start
                __dict['Interval_length'] = 0
        elif idx == len(intronStarts) - 1:
            __dict['Feature'] = '3-UTR'
            if row['strand'] == '+' and end - start == -1:
                __dict['End'] = start
                __dict['Interval_length'] = 0
        if end < start:
            complexExonCount += 1
            continue
        intron_dict.append(__dict)

    combined_dict = exon_dict + intron_dict

    return pd.DataFrame(combined_dict)

def merge_ncbi_id(df):
    
    dfs = []
    for _, grp_df in df.groupby(["Chr", 'Start', 'End', 'Symbol', 'Feature', 'Strand', 'Interval_length'], sort=False):
        if grp_df.shape[0] == 1:
            dfs.extend(grp_df.to_dict('records'))
        else:
            ncbi_list = grp_df['NCBI_ID'].tolist()
            ori = grp_df.iloc[0, :]
            modified = {'Chr': ori['Chr'], 
                        'Start': ori['Start'], 
                        'End': ori['End'], 
                        'Symbol': ori['Symbol'], 
                        'Feature': ori['Feature'], 
                        'Strand': ori['Strand'], 
                        'Interval_length': ori['Interval_length'], 
                        'NCBI_ID': ';'.join(ncbi_list)}
            dfs.append(modified)

    return pd.DataFrame.from_records(dfs)

def update_hgnc(df, url = "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt"):

    logger = logging.getLogger("root")

    def searchSymbol(row, ref):

        symbol = row["name2"]
        calls = ref[ (ref["symbol"] == symbol) | (ref["alias_symbol"].str.contains(symbol + '|', regex=False))
                    | (ref["alias_symbol"].str.endswith(symbol)) | (ref["prev_symbol"].str.contains(symbol + '|', regex=False))
                    | (ref["prev_symbol"].str.endswith(symbol))]

        calls = calls[calls.apply(lambda row: symbol in row['alias_symbol'].split("|") + row['prev_symbol'].split("|"), axis=1)]
        try:
            ### No match found
            if calls.shape[0] == 0:
                return symbol

            ### Return symbol if there is a match
            if calls[calls["symbol"] == symbol].shape[0] >= 1:
                return symbol
            else:
                out_symbols = list(set(calls["symbol"].tolist()))
                if len(out_symbols) > 1:
                    logger.warning(f"Multiple HGNC symbols for {symbol} : {out_symbols}")
                    return symbol
                else:
                    return out_symbols[0]
        except:
            raise ValueError(f"Unrecognized symbol - {symbol}. Matched entries: {calls}")

    # Download data from EMBL-EBI
    logger.info(f"Fetching HGNC symbols from : {url}")
    hgnc_raw = pd.read_csv(url, sep="\t", low_memory=False, usecols=["symbol", "alias_symbol",
                                                                      "prev_symbol"])
    hgnc_raw = hgnc_raw.fillna("NONE")

    df["name2"] = df.parallel_apply(lambda row: searchSymbol(row, hgnc_raw), axis = 1)

    return df

def write_regions(df, outd):
    '''Writes dataframe of exons only, introns only and all'''

    import os
    logger = logging.getLogger("root")


    exon_df = df[df['Feature'].str.startswith('exon_')]
    intron_df = df[~df['Feature'].str.startswith('exon_')]

    exon_df.to_csv(os.path.join(outd, 'refgene_latest_anno_exon.bed'), sep = '\t', index=False )
    intron_df.to_csv(os.path.join(outd, 'refgene_latest_anno_intron.bed'), sep = '\t', index=False )
    df.to_csv(os.path.join(outd, 'refgene_latest_anno.bed'), sep = '\t', index=False )
    logger.info(f"Written all genomic coordinates to {outd} . ")

    return

def report_coor(outd, genePanel_fp=None):

    ### Fetch exon coordinates from RefGene
    refGene_url = "http://api.genome.ucsc.edu/getData/track?genome=hg19;track=refGene"
    logger.info(f"Fetching exon coordinates from {refGene_url}")
    dict_lists = pd.read_json(refGene_url, orient='records')['refGene'].tolist()
    df_list = [ pd.DataFrame.from_records(dict_list) for dict_list in dict_lists ]
    refGene_anno_df = pd.concat(df_list, ignore_index=True)

    ### Update HGNC symbols first
    updated_refGene_anno_df = update_hgnc(refGene_anno_df)
    if genePanel_fp:
        # Sanitize gene list
        with open(genePanel_fp, "r") as ifs:
            genes = list(set([gene.strip().split(u"\xa0")[0] for gene in ifs.readlines() if gene and gene.isupper() ]))
            logger.info(f"Loaded {len(genes)} genes from gene panel. ")
        updated_refGene_anno_df = updated_refGene_anno_df[updated_refGene_anno_df["name2"].isin(genes)]

    ### Compute genomic coordinates
    logger.info("Computing genomic coordinates ... ")
    df = pd.concat(refGene_anno_df.parallel_apply(get_coor, axis=1).to_numpy(), ignore_index=True)
    logger.debug("Sorting table by genomic coordinates ... ")
    df = df.sort_values(['Chr', 'Start', 'End'], ascending = True).reset_index(drop=True)
    logger.debug(f"Extracted {df.shape[0]} genomic regions . ")
    # Merge NCBI ID (usually seen in introns)
    logger.info("Merging NCBI IDs ... ")
    df = merge_ncbi_id(df)

    write_regions(df, outd)

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, help="Output directory", required=True)
    parser.add_argument("-v", "--verbose", type=str, default="INFO", help="Verbosity level")
    parser.add_argument("-l", "--target_genes", type=str, help="Target gene list (one gene per row)", default=None)
    args = parser.parse_args()

    report_coor(args.outdir, genePanel_fp=args.target_genes)


