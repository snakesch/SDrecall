def annotate_from_ref(INPATH: str, REF_REGIONS: str, OUTPATH: str, logLevel="INFO"):
    
    import subprocess
    import logging
    import sys
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = logLevel.upper())

    if subprocess.run("which bedtools", shell=True, capture_output=True).returncode != 0:
        raise ModuleNotFoundError("BEDTools not detected.")
    code = subprocess.run("tail -n+2 {} | bedtools intersect -a {} -b - -wo > {}".format(REF_REGIONS, INPATH, OUTPATH), shell=True).returncode
    if code != 0:
        logging.error("Error in bedtools intersect.")
        sys.exit(-1)
       
    return

def update_bed_file(bed_df, bed_path, logLevel="INFO"):
    
    import pandas as pd
    import subprocess
    import logging
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = logLevel.upper())
    
    bed_df.drop_duplicates().to_csv(bed_path, sep="\t", index=False, header=False)
    cmd = "cat {tbed} | sort -V -k1,3 | cut -f 1-3 | bedtools merge -i - > {tbed}.tmp && mv {tbed}.tmp {tbed}".format(tbed=bed_path)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        logging.error("Error in creating " + bed_path + ".")
        sys.exit(-1)

    return

def merge_bed_intervals(bed_path, logLevel="INFO"):
    
    import subprocess
    import logging
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = logLevel.upper())
    
    cmd = "bedtools sort -i {bed} | bedtools merge -i - > {bed}.tmp && mv {bed}.tmp {bed}".format(bed=bed_path)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        logging.error("Error in merging BED intervals.")
        sys.exit(-1)
        
    return   

def one_way(two_way_df, col_1 = [0, 1, 2, 16], col_2 = [3, 4, 5, 16]):
    '''
    We assume gene name is on the 17th column (0-based). Generate one-way map.
    '''
    
    import pandas as pd
    
    df_1 = two_way_df.iloc[:, col_1]
    df_1.columns = ["chr", "start", "end", "gene"]
    df_2 = two_way_df.iloc[:, col_2]
    df_2.columns = ["chr", "start", "end", "gene"]

    return pd.concat([df_1, df_2], axis=0).drop_duplicates().reset_index(drop=True)

def condense(df, outpath="", gene_col=17):
    
    import pandas as pd
    
    cols = df.columns.tolist()
    by_pairs = df.groupby(by=cols[:6])
    condensed_df = by_pairs.apply(condense_per_group, cols, gene_col)
    if len(outpath) > 0:
        condensed_df.to_csv(outpath, sep='\t', mode="a", header=False, index=False)
    else:
        return condensed_df
    
def condense_per_group(group_df, cols, gene_col):
    
    import pandas as pd
    
    condense_col = cols[gene_col-1]
    group_df[condense_col] = ",".join(list(dict.fromkeys(group_df[condense_col].tolist())))
    return group_df.drop_duplicates()
   
def generate_bed_per_gene(group_df, gene_col, outpath):
    
    import pandas as pd
    import os
    
    if group_df.shape[0] == 0:
        return
    gene = group_df.iloc[0, gene_col].replace(",", "_")
    
    # Homologous regions
    total_bed = one_way(group_df, col_1=[0, 1, 2, gene_col], col_2=[3, 4, 5, gene_col])
    total_bed.iloc[:,[1,2]] = total_bed.iloc[:,[1,2]].astype(int)
    update_bed_file(total_bed, os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"))
    merge_bed_intervals(os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"))
    
    # Principal components
    selected_bed = group_df.iloc[:, [0,1,2,gene_col]]
    selected_bed.iloc[:,[1,2]] = selected_bed.iloc[:,[1,2]].astype(int)
    update_bed_file(selected_bed, os.path.join(outpath, "principal_components", gene + ".bed"))
    merge_bed_intervals(os.path.join(outpath, "principal_components", gene + ".bed"))
    
    out = pd.read_csv(os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"), sep="\t", header=None)

    return out

