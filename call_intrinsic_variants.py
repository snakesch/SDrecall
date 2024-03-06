#!/usr/bin/env python3

import logging
logger = logging.getLogger("root")

def getSubseq(all_homo_bed, read_fq_path, mate_fq_path, reference_genome, avg_frag_size) -> None:
        
        import pyfaidx as pf
        import pybedtools as pb
        import pandas as pd
        
        # Here the input bed file should be the bed file of all homologous reegions
        if avg_frag_size <= 100:
            raise ValueError(f"Read length is too small (<100). ")
        avg_frag_size = avg_frag_size + 1 if avg_frag_size % 2 == 1 else avg_frag_size
        
        # Build static windows for realignment
        bed_obj = pb.BedTool(all_homo_bed)
        window_bed_obj = bed_obj.window_maker(b=bed_obj, w=avg_frag_size)
        truncated_ends = [ idx for idx, iv in enumerate(window_bed_obj) if iv.length < 100]
        window_bed_df = window_bed_obj.to_dataframe()
        for truncated_end in truncated_ends:
            if window_bed_df.iloc[truncated_end - 1, 0] == window_bed_df.iloc[truncated_end, 0]:
                window_bed_df.iloc[truncated_end - 1, 2] = window_bed_df.iloc[truncated_end, 2]
        window_bed_obj = pb.BedTool.from_dataframe(window_bed_df.drop(truncated_ends))
        assert bed_obj.total_coverage() == window_bed_obj.total_coverage(), f"Error in getSubseq: Inconsistent coverages (all_homo_bed: {bed_obj.total_coverage()}; windowed_bed: {window_bed_obj.total_coverage()})"
        
        genome = pf.Fasta(reference_genome, indexname=reference_genome + ".fai", sequence_always_upper=True)
        
        read_fs, mate_fs = open(read_fq_path, "w"), open(mate_fq_path, "w")
        for iv in window_bed_obj:
            cur_seq = genome.get_seq(iv.chrom, iv.start + 1, iv.end) # offset by 1
            read, mate = cur_seq[:150], cur_seq[-150:].complement[::-1]
            read_fs.write(repr(read))
            read_fs.write("\n+\n" + "F" * 150 + "\n")
            mate_fs.write(repr(mate))
            mate_fs.write("\n+\n" + "F" * 150 + "\n")                
        read_fs.close()
        mate_fs.close()
        
        return read_fq_path, mate_fq_path