class Genome:
    """
    A class implementation for genome masking.
    
    Usage:
    rg = Genome(path)
    rg.mask(bedf) # bedf needs not to be complemented
    
    """
       
    def __init__(self, path):
        self.faiIndex = self._fai
        self.dictIndex = self._dictIndex
        self.path = path
        self.total_bed = self._total_bed
        
    def _fai(self):
        
        import os
        from .utils import executeCmd
        
        if not os.path.exists(os.path.join(self.path + ".fai")):
            executeCmd(f"samtools faidx {self.path} ")
        return os.path.join(self.path + ".fai")
    
    def _dictIndex(self):
        
        import os
        from .utils import executeCmd
        
        dict_path = ".".join(self.path.split(".")[:-1]) + ".dict"
        if not os.path.exists(dict_path):
            executeCmd(f"gatk CreateSequenceDictionary -R {self.path} -O {dict_path}")
        
        return dict_path

    def _total_bed(self):
        
        import os
        
        total_bed_path = ".".join(self.path.split(".")[:-1]) + ".bed"
        return total_bed_path if os.path.exists(total_bed_path) else ""        
    
    def getContigGenome(self):
        
        import pandas as pd
        
        contig_genome_fn = ".".join(self.path.split(".")[:-1]) + ".contigsize.genome"
        df = pd.read_csv(self.faiIndex(), sep="\t", header=None).iloc[:, [0,1]]
        df.to_csv(contig_genome_fn, sep="\t", index=False, header=False)
        
    def mask(self, bedf: str):
        
        import os
        import pandas as pd
        from .utils import executeCmd
        
        if not self.total_bed:
            raise FileNotFoundError(f"Please provide a BED file of the genome.")
        
        if not os.path.exists(bedf):
            raise FileNotFoundError(f"Invalid BED file : {bedf} ")
        
        _base_fn = os.path.basename(bedf)
        if "_related_homo_regions" in _base_fn:
            region = _base_fn.split("_r")[0]
        else:
            region = _base_fn.split(".")[-1]

        # Get target contigs
        contig_out = bedf + ".contig"
        df = pd.read_csv(bedf, sep="\t", header=None).iloc[:, 0]
        df.drop_duplicates().to_csv(contig_out, index=False, header=False)
        
        # Create complementary BED file
        masked_genome_tmp = ".".join(self.path.split(".")[:-1]) + "_" + region + "_masked.fasta.tmp"
        masked_genome = ".".join(self.path.split(".")[:-1]) + "_" + region + "_masked.fasta"
        cbed = complement_bed(bedf, self)
        cmd = f"seqtk subseq -l 60 {self.path} {contig_out} > {masked_genome_tmp} "
        executeCmd(cmd)
        
        cmd = f"seqtk seq -l 60 -M {cbed} -n N {masked_genome_tmp} > {masked_genome} "
        executeCmd(cmd)
        
        os.remove(masked_genome_tmp)
        os.remove(contig_out)
        if os.path.exists(self.path + ".seqkit.fai"):
            os.remove(self.path + ".seqkit.fai")

        return masked_genome

def complement_bed(bedf: str, rg: Genome) -> str:
    """
    This function takes the path of a BED file and makes complement of it.

    Return value is the path of complement BED.
    """

    import os
    from .utils import executeCmd

    contig_genome = ".".join(rg.path.split(".")[:-1]) + ".contigsize.genome"
    if not os.path.exists(contig_genome):
        rg.getContigGenome()

    _merged_out = bedf[:-3] + "merged.bed.tmp"
    executeCmd(f"bedtools slop -b 100 -g {contig_genome} -i {bedf} | bedtools merge -d 150 -i - > {_merged_out} ") 
    
    merged_out = bedf[:-3] + "merged.bed"
    cmd = f"bedtools subtract -a {rg.total_bed()} -b {_merged_out} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {merged_out} "
    executeCmd(cmd)
    
    os.remove(_merged_out)

    return merged_out