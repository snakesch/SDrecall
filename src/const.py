import os
import uuid
import re
from typing import Dict, List, Optional, Tuple

shell_utils = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'shell_utils.sh')

class SDrecallPaths:
    """
    Class to manage file paths and naming conventions for SDrecall.
    
    This centralizes all path handling to make maintenance easier and 
    ensure consistency across different modules.
    """
    
    # Class-level variable to hold the singleton instance
    _instance = None
    
    @classmethod
    def get_instance(cls):
        """Return the singleton instance, or raise error if not initialized"""
        if cls._instance is None:
            raise RuntimeError("SDrecallPaths not initialized. Call initialize() first.")
        return cls._instance
    
    @classmethod
    def initialize(cls, 
                  ref_genome: str, 
                  input_bam: str, 
                  reference_sd_map: str, 
                  output_dir: str,
                  target_bed: str = "",
                  sample_id: str = None, 
                  target_tag: str = None):
        """
        Initialize the singleton instance with the provided parameters.
        
        Args:
            ref_genome: Path to reference genome
            input_bam: Path to input BAM file
            reference_sd_map: Path to reference SD map
            output_dir: Parent directory where all results will be stored
            target_bed: Optional path to target BED file
            sample_id: Optional sample identifier (extracted from BAM if not provided)
            target_tag: Optional target region tag (extracted from target_bed if not provided)
        
        Returns:
            The initialized SDrecallPaths instance
        """
        cls._instance = cls(
            ref_genome=ref_genome,
            input_bam=input_bam,
            reference_sd_map=reference_sd_map,
            output_dir=output_dir,
            target_bed=target_bed,
            sample_id=sample_id,
            target_tag=target_tag
        )
        return cls._instance
    
    def __init__(self, 
                ref_genome: str, 
                input_bam: str, 
                reference_sd_map: str,
                output_dir: str,
                target_bed: str = "",
                sample_id: str = None, 
                target_tag: str = None):
        """
        Initialize the path manager with key inputs.
        
        Args:
            ref_genome: Path to reference genome
            input_bam: Path to input BAM file
            reference_sd_map: Path to reference SD map
            output_dir: Parent directory where all results will be stored
            target_bed: Optional path to target BED file
            sample_id: Optional sample identifier (extracted from BAM if not provided)
            target_tag: Optional target region tag (extracted from target_bed if not provided)
        """
        # Store input file paths
        self.ref_genome = os.path.abspath(ref_genome)
        self.input_bam = os.path.abspath(input_bam)
        self.reference_sd_map = os.path.abspath(reference_sd_map)
        self.target_bed = os.path.abspath(target_bed) if target_bed else ""
        self.output_dir = os.path.abspath(output_dir)
        
        # Extract assembly version primarily from reference_sd_map path
        self.assembly = self._extract_assembly_version(reference_sd_map, ref_genome)
        
        # Extract sample ID
        self.sample_id = sample_id or self._extract_sample_id(input_bam)
        
        # Extract target tag
        self.target_tag = target_tag or self._extract_target_tag(target_bed)
        
        # Create base directory name based on sample_id, assembly, and target_tag
        dir_parts = [self.sample_id, self.assembly]
        if self.target_tag:
            dir_parts.append(self.target_tag)
        
        # Generate base directory name under output_dir
        base_name = "_".join(dir_parts)
        self.work_dir = os.path.join(self.output_dir, base_name)
        
        # Make sure the work_dir exists
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Create standard directory structure
        self.dirs = self._initialize_directories()
        
        # Store the basename for file naming
        self.basename = os.path.basename(self.input_bam).replace(".bam", "")
    
    def _extract_assembly_version(self, reference_sd_map: str, ref_genome: str) -> str:
        """Extract assembly version from reference SD map path and reference genome"""
        # Try to detect from reference_sd_map
        if "hg19" in reference_sd_map or "GRCh37" in reference_sd_map:
            return "hg19"
        elif "hg38" in reference_sd_map or "GRCh38" in reference_sd_map:
            return "hg38"
        
        # As fallback, try ref_genome
        if "hg19" in ref_genome or "GRCh37" in ref_genome:
            return "hg19"
        elif "hg38" in ref_genome or "GRCh38" in ref_genome:
            return "hg38"
        
        # Default to extracting from filename
        filename = os.path.basename(ref_genome)
        if "hg19" in filename or "GRCh37" in filename:
            return "hg19"
        elif "hg38" in filename or "GRCh38" in filename:
            return "hg38"
        else:
            raise ValueError(f"Could not determine assembly version from {reference_sd_map} or {ref_genome}")
    
    def _extract_sample_id(self, input_bam: str) -> str:
        """Extract sample ID from BAM filename"""
        bam_name = os.path.basename(input_bam)
        # Just use filename without extension
        return bam_name.split('.')[0]
    
    def _extract_target_tag(self, target_bed: str) -> str:
        """Extract target tag from target BED filename"""
        if not target_bed:
            return "whole_genome"
            
        bed_name = os.path.basename(target_bed)
        # Get first part of filename delimited by dot
        return bed_name.split('.')[0]
    
    def _initialize_directories(self) -> Dict[str, str]:
        """Set up the standard directory structure"""
        dirs = {
            "root": self.work_dir,
            "recall_results": os.path.join(self.work_dir, "recall_results")
        }
        
        # Create directories if they don't exist
        for dir_path in dirs.values():
            os.makedirs(dir_path, exist_ok=True)
            
        return dirs
    
    # File path getter methods
    def get_multi_align_bed_path(self) -> str:
        """Get path for multi-aligned regions BED file"""
        return os.path.join(self.dirs["intermediate"], f"{self.basename}.{self.target_tag}.multialign.bed")
    
    def get_raw_sd_binary_map_path(self) -> str:
        """Get path for raw SD binary map"""
        return os.path.join(self.dirs["intermediate"], "raw_SD_binary_map.tsv")
    
    def get_filtered_sd_binary_map_path(self) -> str:
        """Get path for filtered SD binary map"""
        return os.path.join(self.dirs["intermediate"], "filtered_SD_binary_map.tsv")
    
    def get_multiplex_graph_path(self) -> str:
        """Get path for multiplex graph"""
        return os.path.join(self.dirs["intermediate"], f"{self.basename}-multiplexed_homologous_sequences.graphml")
    
    def get_annotated_graph_path(self) -> str:
        """Get path for annotated graph"""
        return self.get_multiplex_graph_path().replace(".graphml", ".trim.annoPC.graphml")
    
    def get_directed_graph_path(self, chrom: str) -> str:
        """Get path for directed graph for a specific chromosome"""
        return self.get_multiplex_graph_path().replace(".graphml", f".directed.overlap.{chrom}.graphml")
    
    def get_target_overlapping_query_sd_bed_path(self) -> str:
        """Get path for target overlapping query SD BED"""
        return os.path.join(self.dirs["intermediate"], f"{self.basename}target_overlapping_query_SD.bed")
    
    def get_rg_dir(self, rg_label: str) -> str:
        """Get directory for a specific PC cluster"""
        rg_dir = os.path.join(self.dirs["rg_clusters"], rg_label)
        os.makedirs(rg_dir, exist_ok=True)
        return rg_dir
    
    def get_rg_counterparts_dir(self, rg_label: str) -> str:
        """Get directory for PC counterparts"""
        counterparts_dir = os.path.join(self.get_rg_dir(rg_label), f"{rg_label}_counterparts")
        os.makedirs(counterparts_dir, exist_ok=True)
        return counterparts_dir
    
    def get_rg_all_dir(self, rg_label: str) -> str:
        """Get directory for all PC-related regions"""
        all_dir = os.path.join(self.get_rg_dir(rg_label), f"{rg_label}_all")
        os.makedirs(all_dir, exist_ok=True)
        return all_dir
    
    def get_rg_main_dir(self, rg_label: str) -> str:
        """Get main directory for PC"""
        main_dir = os.path.join(self.get_rg_dir(rg_label), rg_label)
        os.makedirs(main_dir, exist_ok=True)
        return main_dir
    
    def get_rg_bed_path(self, rg_label: str) -> str:
        """Get path for PC bed file"""
        return os.path.join(self.get_rg_main_dir(rg_label), f"{rg_label}.bed")
    
    def get_counterparts_bed_path(self, rg_label: str) -> str:
        """Get path for counterparts bed file"""
        return os.path.join(self.get_rg_counterparts_dir(rg_label), f"{rg_label}_counterparts_regions.bed")
    
    def get_all_homo_regions_bed_path(self, rg_label: str) -> str:
        """Get path for all homologous regions bed file"""
        return os.path.join(self.get_rg_all_dir(rg_label), f"{rg_label}_related_homo_regions.bed")
    
    def get_all_homo_regions_fastq_path(self, rg_label: str) -> str:
        """Get path for all homologous regions fastq file"""
        return os.path.join(self.get_rg_all_dir(rg_label), f"{rg_label}_related_homo_regions.raw.fastq")
    
    def get_masked_genome_path(self, rg_label: str) -> str:
        """Get path for masked genome"""
        return os.path.join(self.get_rg_main_dir(rg_label), f"{rg_label}.masked.fasta")
    
    def get_masked_genome_fai_path(self, rg_label: str) -> str:
        """Get path for masked genome FAI index"""
        return f"{self.get_masked_genome_path(rg_label)}.fai"
    
    def get_masked_genome_contigsize_path(self, rg_label: str) -> str:
        """Get path for masked genome contig size file"""
        return os.path.join(self.get_rg_main_dir(rg_label), f"{rg_label}.masked.contigsize.genome")
    
    def get_minimap_index_path(self, rg_label: str) -> str:
        """Get path for minimap2 index"""
        return os.path.join(self.get_rg_main_dir(rg_label), f"{rg_label}.masked.mmi")
    
    def get_intrinsic_sam_path(self, rg_label: str) -> str:
        """Get path for intrinsic alignment SAM"""
        return os.path.join(self.get_rg_main_dir(rg_label), f"{rg_label}.raw.sam")