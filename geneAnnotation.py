#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse as ap
import logging
import sys
import re
import gc

def update_hgnc(table, d=";", label=None, url="ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt"):
    # a stands for the absolute path to the annotation table, if not inputted, we'll use a default URL
    # b stands for the absolute path to the to_be_modified table, or it can be a file containing a list of genes with one gene per row. b can also be a DataFrame Object
    # c stands for the column label in the to_be_modified table where we compare gene symbols. If b is a gene list file, do not input c, let it be None
    # d stands for the delimiter used in gene symbol field

    try:
        hgnc_anno_table = pd.read_csv(url, sep='\t', low_memory=False)
    except:
        logging.critical("Failed to fetch data from default URL.")
        sys.exit(-1)
    if label:
        if isinstance(table, pd.DataFrame):
            to_be_modified_table = table
        elif isinstance(table, str):
            to_be_modified_table = pd.read_csv(table, sep='\t', low_memory=False)
    else:
        to_be_modified_table = pd.read_csv(table, sep='\t', low_memory=False)
        gene_symbol_regex = re.compile("^[A-Z]{1}[A-Z0-9]{2,}$|^C[0-9]+orf[0-9]+$")
        label = [ label for label in to_be_modified_table.columns.tolist() if re.search(gene_symbol_regex, str(to_be_modified_table.at[1,label]))][0]

    new_table = pd.DataFrame()
    new_table['original_symbols'] = list(dict.fromkeys(to_be_modified_table[label]))
    new_table['HGNC_symbol'] = new_table.apply(convert_per_row, args=(hgnc_anno_table, 'original_symbols', d), axis=1)
    new_table.set_index('original_symbols', inplace=True)
    to_be_modified_table = to_be_modified_table.merge(new_table, left_on=label, right_index=True)
    new_table = None
    del new_table
    gc.collect()

    logging.info("Now we have a new column recording the updated HGNC_symbol in: \n" + str(to_be_modified_table[:4].to_string(index=False)))
    original_column_index = to_be_modified_table.columns.get_loc(label)
    to_be_modified_table.drop(columns=label, inplace=True)
    tobe_mov = to_be_modified_table.pop('HGNC_symbol')
    to_be_modified_table.insert(loc=original_column_index, column=label, value=tobe_mov)
    if label == "Later_defined_header":
        to_be_modified_table.to_csv(table, sep='\t', index=False, header=False)
    elif not isinstance(table, pd.DataFrame):
        to_be_modified_table.to_csv(table, sep='\t', index=False)
    else:
        return to_be_modified_table

def get_intron_coor(row):
    """
    This function generates two columns: intron start coordinates and intron end coordinates.
    """
    exon_starts = list(map(int, row['exonStarts'][:-1].split(',')))
    exon_ends = list(map(int, row['exonEnds'][:-1].split(',')))
    # Determine whether the record represents a non-coding gene. (cdsStart == cdsEnd)
    cds_available = int(row['cdsStart']) < int(row['cdsEnd'])
    __pre, __cur = 0, 1
    consec = False
    intron_starts = []
    intron_ends = []
    while __cur < len(exon_starts):
        if exon_starts[__cur] > exon_ends[__pre] + 1:
            intron_starts.append(exon_ends[__pre] + 1)
            intron_ends.append(exon_starts[__cur] - 1)
            if consec:
                consec = False
        elif consec:
            pass
        else:
            consec = True
        __pre += 1
        __cur += 1
    if row['strand'] == '+':
        if cds_available:
            utr_5_starts = row['txStart']
            utr_5_ends = row['cdsStart'] - 1
            utr_3_starts = row['cdsEnd'] + 1
            utr_3_ends = row['txEnd']
        else:
            utr_5_starts = row['txStart']
            utr_5_ends = row['txStart']
            utr_3_starts = row['cdsEnd'] + 1
            utr_3_ends = row['txEnd']
    elif row['strand'] == '-':
        intron_starts.reverse()
        intron_ends.reverse()
        if cds_available:
            utr_5_starts = row['cdsEnd'] + 1
            utr_5_ends = row['txEnd']
            utr_3_starts = row['txStart']
            utr_3_ends = row['cdsStart'] - 1
        else:
            utr_5_starts = row['cdsEnd'] + 1
            utr_5_ends = row['txEnd']
            utr_3_starts = row['txStart']
            utr_3_ends = row['txStart']
    # Insert utr coordinates into the intron coor list
    intron_starts.insert(0, utr_5_starts)
    intron_starts.append(utr_3_starts)
    intron_ends.insert(0, utr_5_ends)
    intron_ends.append(utr_3_ends)
    logging.info("For this row {ind}, the returned intron starts are: \n{sta}\nthe returned intron ends are: \n{ends}\n\n".format(ind=row.name, sta=",".join(map(str,intron_starts)), ends=",".join(map(str,intron_ends))))
    return ",".join(map(str,intron_starts)), ",".join(map(str,intron_ends))


# Define a function that generate 8-column df unit based on intron or exon coor pairs on each row.
def generate_8_column_df_unit(row, *column_pairs):
    start_column = [ x for x in column_pairs if re.search("Starts$", x) ][0]
    end_column = [ x for x in column_pairs if re.search("Ends$", x) ][0]

    # Split coordinates into lists. Check whether coordinates string end with comma ",".
    if re.search(",$",row[start_column]):
        start_list = row[start_column][:-1].split(',')
    else:
        start_list = row[start_column].split(',')
    if re.search(",$",row[end_column]):
        end_list = row[end_column][:-1].split(',')
    else:
        end_list = row[end_column].split(',')

    # Pair-up start_end coordinates in dict
    dict_list = []
    if re.search("^exon", start_column):
        if row['strand'] == '-':
            start_list.reverse()
            end_list.reverse()
        for i in range(0, len(start_list)):
            df_dict = {"Chr": row['chrom'], "Start": start_list[i], "End": end_list[i], "Symbol": row['name2'], "Feature":"exon_" + str(i+1), "Strand": row['strand'], "Interval_length":int(end_list[i]) - int(start_list[i]), "NCBI_ID": row['name']}
            dict_list.append(df_dict)
        # convert the dict_list object into dataframe
        df = pd.DataFrame(dict_list)
    else:
        for i in range(0, len(start_list)):
            if i == 0:
                df_dict = {"Chr": row['chrom'], "Start": start_list[i], "End": end_list[i], "Symbol": row['name2'], "Feature":"5-UTR", "Strand": row['strand'], "Interval_length":int(end_list[i]) - int(start_list[i]), "NCBI_ID": row['name']}
            elif i == len(start_list) - 1:
                df_dict = {"Chr": row['chrom'], "Start": start_list[i], "End": end_list[i], "Symbol": row['name2'], "Feature":"3-UTR", "Strand": row['strand'], "Interval_length":int(end_list[i]) - int(start_list[i]), "NCBI_ID": row['name']}
            else:
                df_dict = {"Chr": row['chrom'], "Start": start_list[i], "End": end_list[i], "Symbol": row['name2'], "Feature":"intron_" + str(i), "Strand": row['strand'], "Interval_length":int(end_list[i]) - int(start_list[i]), "NCBI_ID": row['name']}
            dict_list.append(df_dict)
        # convert the dict_list object into dataframe
        df = pd.DataFrame(dict_list).drop_duplicates()

    return df

# Define a function that checks the error UTR coordinates.
def modify_utr_coor(row):
    if int(row['Interval_length']) == -1:
        logging.warning("This row {} start coor is bigger than end coor".format(row.name))
        if row['Strand'] == '+' and row['Feature'] == '3-UTR':
            return str(int(row['Start'])-2), row['End'], 0
        elif row['Strand'] == '-' and row['Feature'] == '5-UTR':
            return str(int(row['Start'])-2), row['End'], 0
        else:
            return row['Start'], str(int(row['End'])+2), 0
    else:
        return row['Start'], row['End'], row['Interval_length']

def merge_ncbi_id(row, df, dicts):
    """
    This function merges rows that only differ by NCBI ID (gene name). Input df should have been sorted by Chr, Start, End, Symbol, Feature and Strand.
    """
    if len(dicts) >= 1:
        dic = dicts[-1]
        checked = ((dic['Chr'] == row['Chr']) & (dic['Start'] == row['Start']) & (dic['End'] == row['End']) & (dic['Symbol'] == row['Symbol']) & (dic['Feature'] == row['Feature']) & (dic['Strand'] == row['Strand']))

        if not checked:
            df_section = df.iloc[row.name:row.name+100, :]
            bool_array = (df_section['Chr'] == row['Chr']) & (df_section['Start'] == row['Start']) & (df_section['End'] == row['End']) & (df_section['Symbol'] == row['Symbol'])
            unit_df = df_section.loc[bool_array, :]
            row_dict = {'Chr': row['Chr'], 'Start':row['Start'], 'End':row['End'], 'Symbol':row['Symbol'], 'Feature':row['Feature'], 'Strand':row['Strand'], 'Interval_length':row['Interval_length'], 'NCBI_ID': ';'.join(unit_df['NCBI_ID'].tolist())}
            logging.info("The merged NCBI ID is: " + str(row_dict['NCBI_ID']))
            dicts.append(row_dict)
            if len(dicts) >= 2:
                dicts = dicts[1:]
            return row_dict
        else:
            logging.info("We found that this row has already been recorded for NCBI_ID merging, short circuiting works. Skip this row {}".format(row.name))
            return
    else:
        df_section = df.iloc[row.name:row.name+100, :]
        bool_array = (df_section['Chr'] == row['Chr']) & (df_section['Start'] == row['Start']) & (df_section['End'] == row['End']) & (df_section['Symbol'] == row['Symbol'])
        unit_df = df_section.loc[bool_array, :]
        row_dict = {'Chr': row['Chr'], 'Start':row['Start'], 'End':row['End'], 'Symbol':row['Symbol'], 'Feature':row['Feature'], 'Strand':row['Strand'], 'Interval_length':row['Interval_length'], 'NCBI_ID': ';'.join(unit_df['NCBI_ID'].tolist())}
        logging.info("The merged NCBI ID is: " + str(row_dict['NCBI_ID']))
        dicts.append(row_dict)
        if len(dicts) >= 2:
            dicts = dicts[1:]
        return row_dict

def convert_per_row(row, table, label, delimiter):
    # Select the cell where it matches with the Gene symbol in to_be_modified_table
    selected_record = table.loc[table['symbol'].isin(str(row[label]).split(delimiter)), 'symbol']
    # Convert the series into list
    # If we found one match
    if len(selected_record.tolist()) == 1:
        return selected_record.tolist()[0]
    # If we found more than one match
    elif len(selected_record.tolist()) > 1:
        logging.info("This record (index {}) contains multiple updated HGNC symbols {}".format(row.name, str(row[label])))
        return delimiter.join(selected_record.tolist())
    # If we haven't found one match, it means the symbol in to_be_modified table is an alias.
    else:
        logging.info("This record (index {}) does not contain an updated HGNC symbol {}".format(row.name, str(row[label])))
        # Prepare a series of booleans based on checking if there are intersection
        # bool_ser_alias = table['alias_symbol'].apply(intersection, args=str(row[label]).split(delimiter))
        bool_ser_alias = table['alias_symbol'].str.strip('"').str.split("|", expand=True).isin(str(row[label]).split(delimiter)).any(axis=1)
        # bool_ser_prev = table['prev_symbol'].apply(intersection, args=str(row[label]).split(delimiter))
        bool_ser_prev = table['prev_symbol'].str.strip('"').str.split("|", expand=True).isin(str(row[label]).split(delimiter)).any(axis=1)
        alias_or_prev_record_list = table.loc[bool_ser_alias | bool_ser_prev, 'symbol'].tolist()
        if len(alias_or_prev_record_list) == 1:
            return alias_or_prev_record_list[0]
        elif len(alias_or_prev_record_list) > 1:
            # We may found the name match with one prev symbol for one gene and one alias symbol for another gene
            # We need to determine the true gene for the current symbol
            search_symbols = str(row[label]).split(delimiter)
            if bool_ser_alias.any(): alias_symbols = table.loc[bool_ser_alias, 'alias_symbol'].str.split("|", expand=True)
            if bool_ser_prev.any(): prev_symbols = table.loc[bool_ser_prev, 'prev_symbol'].str.split("|", expand=True)
            if 'alias_symbols' in locals() and 'prev_symbols' in locals():
                symbol_dict = {alias_symbols.shape[1]: 'alias_symbol', prev_symbols.shape[1]:'prev_symbol'}
                if symbol_dict[min(alias_symbols.shape[1], prev_symbols.shape[1])] == "alias_symbol":
                    return delimiter.join(table.loc[bool_ser_alias,'symbol'].tolist())
                elif symbol_dict[min(alias_symbols.shape[1], prev_symbols.shape[1])] == "prev_symbol":
                    return delimiter.join(table.loc[bool_ser_prev, 'symbol'].tolist())
            else:
                return delimiter.join(alias_or_prev_record_list)
        else:
            return str(row[label])

def main_convert(output_merged, output_exon, output_intron, target_genes=None):

    dict_lists = pd.read_json('http://api.genome.ucsc.edu/getData/track?genome=hg19;track=refGene', orient='records')['refGene'].tolist()
    # Notice that if we do not need to select chrom, instead we need to download the whole genome, the downloaded json['refGene'].tolist() would be a structure like:
    # [[chr1_dict_list],[chr2_dict_list],...,[chrY_dict_list]]
    df_list = []
    for dict_list in dict_lists:
        refGene_anno_df = pd.DataFrame.from_records(dict_list)
        df_list.append(refGene_anno_df)

    refGene_anno_df = pd.concat(df_list, ignore_index=True)
    try:
        refGene_anno_df['intronStarts'], refGene_anno_df['intronEnds'] = zip(*refGene_anno_df.apply(get_intron_coor, axis=1))
        logging.info(str(refGene_anno_df.head()))
    except KeyError:
        logging.warning("Some column headers are not right, print the raw table headers: " + ",".join(map(str,refGene_anno_df.columns.tolist())))
        logging.warning("Let's also see an example of the raw dict list:" + str([len(dict_list) for dict_list in dict_lists]))
        return None

    exon_df = pd.concat(refGene_anno_df.apply(generate_8_column_df_unit, args=('exonStarts', 'exonEnds'), axis=1).to_numpy(), ignore_index=True)
    intron_df = pd.concat(refGene_anno_df.apply(generate_8_column_df_unit, args=('intronStarts', 'intronEnds'), axis=1).to_numpy(), ignore_index=True)

    # Concat two dfs and sort by coordinates
    total_df = pd.concat([exon_df, intron_df], ignore_index=True)
    total_df = total_df.astype({"Chr":'str', 'Start':'int', 'End':'int', 'Symbol':'str', 'Feature':'str', 'Strand':'str', 'Interval_length':'int', 'NCBI_ID':'str'}, errors='raise')
    total_df.sort_values(by=['Chr', 'Start', 'End'], ascending=True, inplace=True)

    # Correct errors of UTR region
    total_df['Start'], total_df['End'], total_df['Interval_length'] = zip(*total_df.apply(modify_utr_coor, axis=1))

    # Condense df by merging NCBI IDs
    dict_list = []
    total_df.sort_values(by=['Chr', 'Start', 'End', 'Symbol', 'Feature', 'Strand'], ascending=True, inplace=True)
    logging.info(str(total_df.head()))
    df_array = total_df.apply(merge_ncbi_id, args=(total_df, dict_list), axis=1).dropna(how='all').to_numpy()

    # convert list of dicts to dict and split it to two dfs
    merged_df = pd.DataFrame.from_records(df_array)
    merged_df.dropna(how='all', inplace=True)
    merged_df = merged_df.replace(to_replace = "", value = ".")
    logging.info(str(merged_df.head()))
    merged_exon_df = merged_df.loc[merged_df['Feature'].str.startswith("exon"), :]
    merged_intron_df = merged_df.loc[np.logical_not(merged_df['Feature'].str.startswith("exon")), :]

    merged_df = update_hgnc(table=merged_df, label="Symbol", d=";")
    merged_df.to_csv(output_merged, sep='\t', index=False)

    merged_exon_df = update_hgnc(table=merged_exon_df, label="Symbol", d=";")
    merged_exon_df.to_csv(output_exon, sep='\t', index=False)

    merged_intron_df = update_hgnc(table=merged_intron_df, label="Symbol", d=";")
    merged_intron_df.to_csv(output_intron, sep='\t', index=False)

    if target_genes:
        # Use the function will return a file with latest HGNC gene symbol
        # target_genes = update_hgnc(table=target_genes, d=";", label="Symbol")
        with open(target_genes, 'r') as tg:
            genes = tg.readlines()
        genes = set(genes)
        selected_df = merged_df.loc[merged_df['Symbol'].isin(genes), :]
        logging.info(str(selected_df.head()))

        # Isolate the exon part and intron part respectively.
        selected_exon_df = selected_df.loc[selected_df['Feature'].str.startswith("exon"), :]
        selected_intron_df = selected_df.loc[np.logical_not(selected_df['Feature'].str.startswith("exon")), :]

        # Output the selected tables:
        target_gene_group_name = target_genes.split('/')[-1].split('.')[0]
        selected_df.to_csv(output_merged[:-4]+"." + target_gene_group_name + ".txt", sep='\t', index=False)
        selected_exon_df.to_csv(output_exon[:-4]+"." + target_gene_group_name + ".txt", sep='\t', index=False)
        selected_intron_df.to_csv(output_intron[:-4]+"." + target_gene_group_name + ".txt", sep='\t', index=False)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-om", "--output_merged", type=str, help="Output path of merged annotation table (exon & intron)",
                        required=True)
    parser.add_argument("-oe", "--output_exon", type=str, help="Output path of exon annotation table",
                        required=True)
    parser.add_argument("-oi", "--output_intron", type=str, help="Output path of intron annotation table",
                        required=True)
    parser.add_argument("-v", "--verbose", type=str, default="INFO", help="Verbosity level")
    parser.add_argument("-tg", "--target_genes", type=str, help="Target gene list (one gene per row)",
                        required=False, default=None)

    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    main_convert(args.output_merged, args.output_exon, args.output_intron, args.target_genes)
