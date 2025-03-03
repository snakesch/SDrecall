import os
import pandas as pd
import pybedtools as pb

from src.log import logger
from src.const import SDrecallPaths


def stat_all_RG_region_size(sdrecall_paths: SDrecallPaths):
    # Test this function, only takes 3 min to finish
    # use glob module to list out all the bed files under wkd named as "PC*_related_homo_regions.raw.bed"
    beds = sdrecall_paths.all_homo_regions_bed_paths()
    logger.info(f"The beds are {beds}")
    
    beds = sorted(beds, key = lambda f: pb.BedTool(f).sort().total_coverage(), reverse=True)

    records = []
    for bed in beds:
        bedf = pd.read_table(bed, header=None)
        fc_bedf = bedf.loc[bedf.iloc[:, -1].str.contains(r"^FC:"), :].drop_duplicates().sort_values(by=bedf.columns.tolist()[-1])
        sub_ids = [i for i in range(fc_bedf.shape[0])]
        rg_label = os.path.basename(bed).split("_")[0]
        records.append((rg_label,sub_ids))

    df = pd.DataFrame(records, columns=["rg_label", "subgroup_id"])
    df_exploded = df.explode("subgroup_id")
    return df_exploded