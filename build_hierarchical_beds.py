#!/usr/bin/env python3

import logging

logger = logging.getLogger("root")

def establish_beds_per_PC_cluster(label, pc_sd_pair, base_folder, avg_frag_size, std_frag_size):
    '''
    pc_sd_pair: { "PC0": {"PCs": {0: region0, 1: region1, 2: region2}, 
                    "SD_counterparts": {0: [], 1: [], 2: []}} }
    '''
    from itertools import chain
    
    import pybedtools as pb
    
    from utils import init_dtree
    
    # Initialize output directory tree
    path_catalog = init_dtree(base_folder, label)
    logger.debug(f"The PC bed file is {path_catalog['PC_bed']}, the counterparts bed file is {path_catalog['Counterparts_bed']}, the total bed file is {path_catalog['All_region_bed']}")
    
    # Build PC BED objects
    current_pc = pc_sd_pair[label]
    current_pc_list = list(map(lambda x: x.to_tuple(), current_pc["PCs"].values()))
    current_pc_bed_obj = pb.BedTool(current_pc_list)
    current_pc_bed_obj.saveas(path_catalog['PC_bed'].replace(".bed", ".raw.bed"))
    current_pc_bed_obj = current_pc_bed_obj.sort().merge()
    current_pc_bed_obj.saveas(path_catalog['PC_bed'])

    # Build SD counterparts BED object and subtract any overlapping PC
    resolved_counterpart_list = list(chain(*current_pc["SD_counterparts"].values())) 
    current_counterpart_list = list(map(lambda x: x.to_tuple(), resolved_counterpart_list))
    current_counterpart_bed_obj = pb.BedTool(current_counterpart_list)
    
    non_overlapping_counterparts = current_counterpart_bed_obj.subtract(current_pc_bed_obj).sort().merge()
    non_overlapping_counterparts.saveas(path_catalog['Counterparts_bed'])
    
    # Combine PCs and SD counterparts
    current_pc_bed_obj.cat(non_overlapping_counterparts).saveas(path_catalog['All_region_bed'])
                              
    return

def get_mask_genome_per_PC_cluster(label, base_folder, reference_genome, build, avg_frag_size, std_frag_size):
    
    import pybedtools as pb
    import os
    from utils import init_dtree
    
    import pyfaidx as pf
    
    path_catalog = init_dtree(base_folder, label)
    bed_to_mask = path_catalog['PC_bed']
    
    factor = int(avg_frag_size + 2*std_frag_size)
    regions_to_mask = pb.BedTool(bed_to_mask).merge(d=2*factor).slop(b=factor+1000, genome=build).sort()
    
    # Get output paths
    outd, fn = os.path.split(bed_to_mask)
    
    # Reference genome 
    genome = pf.Fasta(reference_genome, indexname=reference_genome + ".fai")
    
    # For each PC, build masked genomes from masked sequences
    masked_genome = open( os.path.join( outd, label + ".masked.fasta"), "w")
    for iv in regions_to_mask.to_dataframe().itertuples():
        iv_seq = genome[iv.chrom][iv.start:iv.end].seq
        mask_seq = "N" * 1000
        masked_iv_seq = mask_seq + iv_seq[1000:-1000] + mask_seq
        assert masked_iv_seq != iv_seq, f"Sequences from {iv.chrom}:{iv.start}-{iv.end} are already flanked by Ns on both ends. "
        target_seq = pf.Sequence(name = iv.chrom, seq = masked_iv_seq, start = iv.start, end = iv.end)
        masked_genome.write(repr(target_seq))
        masked_genome.write("\n")
    masked_genome.close()
    
    # Create index file for masked genomes
    pf.Fasta(os.path.join( outd, label + ".masked.fasta"))
    
    return