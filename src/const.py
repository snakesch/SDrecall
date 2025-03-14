import os
import re
import shutil
import sys
from typing import Dict, List, Set

from src.insert_size import get_insert_size_distribution

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
                  target_tag: str = None,
                  clean_dirs: bool = True):
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
            clean_dirs: Whether to clean existing content in work directories
        
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
            target_tag=target_tag,
            clean_dirs=clean_dirs
        )
        return cls._instance
    
    def __init__(self, 
                ref_genome: str, 
                input_bam: str, 
                reference_sd_map: str,
                output_dir: str,
                target_bed: str = "",
                sample_id: str = None, 
                target_tag: str = None,
                clean_dirs: bool = True):
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
            clean_dirs: Whether to clean existing content in work directories
        """
        # Store input file paths
        self.ref_genome = os.path.abspath(ref_genome)
        self.input_bam = os.path.abspath(input_bam)
        self.reference_sd_map = os.path.abspath(reference_sd_map)
        self.target_bed = os.path.abspath(target_bed) if target_bed else ""
        self.output_dir = os.path.abspath(output_dir)
        self.repo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        
        # Extract assembly version primarily from reference_sd_map path
        self.assembly = self._extract_assembly_version(reference_sd_map, ref_genome)
        self.avg_frag_size, self.median_frag_size, self.frag_size_std = get_insert_size_distribution(self.input_bam)
        
        # Extract sample ID
        self.sample_id = sample_id or self._extract_sample_id(input_bam)
        
        # Extract target tag
        self.target_tag = target_tag or self._extract_target_tag(target_bed)
        
        # Create base directory name based on sample_id, assembly, and target_tag
        dir_parts = [self.sample_id, self.assembly]
        if self.target_tag:
            dir_parts.append(self.target_tag)
        
        # Generate base directory name under output_dir
        base_name = "_".join(dir_parts) + "_SDrecall"
        self.basename = base_name
        self.work_dir = os.path.join(self.output_dir, base_name)
        
        # Make sure the work_dir exists
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Track realign groups
        self.realign_groups = set()
        
        # Create standard directory structure
        self._initialize_directories(clean_dirs)
        
    
    def _extract_assembly_version(self, reference_sd_map: str, ref_genome: str) -> str:
        """Extract assembly version from reference SD map path and reference genome"""
        sd_map_name = os.path.basename(reference_sd_map)
        # Try to detect from reference_sd_map
        if "hg19" in sd_map_name or "GRCh37" in sd_map_name:
            return "hg19"
        elif "hg38" in sd_map_name or "GRCh38" in sd_map_name:
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
            return "exome"
        if "default_target" in target_bed:
            return "exome"
            
        bed_name = os.path.basename(target_bed)
        # Get first part of filename delimited by dot
        return bed_name.split('.')[0]
    
    def _clean_directory(self, directory: str, file_only: bool = False) -> None:
        """Remove all content from a directory"""
        if os.path.exists(directory):
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isfile(item_path):
                    os.unlink(item_path)
                elif os.path.isdir(item_path):
                    if not file_only:
                        shutil.rmtree(item_path)
    
    def _initialize_directories(self, clean_dirs: bool = True) -> Dict[str, str]:
        """
        Set up the standard directory structure
        
        Args:
            clean_dirs: Whether to clean existing directory contents
        """
        dirs = {
            "root": self.work_dir,
            "recall_results": os.path.join(self.work_dir, "recall_results"),
            "realign_groups": os.path.join(self.work_dir, "realign_groups"),
            "intermediates": os.path.join(self.work_dir, "intermediates")
        }
        
        self.recall_results_dir = dirs["recall_results"]
        self.realign_groups_dir = dirs["realign_groups"]
        self.tmp_dir = dirs["intermediates"]
        self.dirs = dirs

        # Create directories if they don't exist and clean if requested
        for dir_name, dir_path in dirs.items():
            os.makedirs(dir_path, exist_ok=True)
            if clean_dirs:  # Don't clean root dir
                if dir_name != "root":
                    self._clean_directory(dir_path, file_only=True)
                else:
                    self._clean_directory(dir_path, file_only=False)
            elif dir_name == "recall_results" or dir_name == "intermediates":
                self._clean_directory(dir_path, file_only=True)

        if not clean_dirs:
            self._discover_and_register_realign_groups()
            
        return dirs
    

    def pooled_raw_bam_path(self) -> str:
        """Path for pooled raw bam file"""
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.pooled.raw.bam")
    
    def recall_raw_vcf_path(self) -> str:
        """Path for recall raw vcf file"""
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.raw.vcf.gz")
    
    def pooled_filtered_bam_path(self) -> str:
        """Path for pooled filtered bam file"""
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.pooled.clean.bam")
    
    def recall_filtered_vcf_path(self) -> str:
        """Path for recall filtered vcf file"""
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.clean.vcf.gz")


    def final_recall_vcf_path(self) -> str:
        """Path for final recall vcf file"""
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.vcf.gz")

    def realign_meta_table_path(self) -> str:
        """Path for realign meta table"""
        return os.path.join(self.work_dir, f"{self.sample_id}_realign_meta_table.tsv")
    

    def ref_genome_fai_path(self) -> str:
        """Path for reference genome FAI index"""
        if os.path.exists(f"{self.ref_genome}.fai") and os.path.getmtime(f"{self.ref_genome}.fai") > os.path.getmtime(self.ref_genome):
            return f"{self.ref_genome}.fai"
        else:
            import subprocess
            subprocess.run(f"samtools faidx {self.ref_genome}", shell=True)
            return f"{self.ref_genome}.fai"
        
    
    def register_realign_group(self, rg_label_or_index) -> None:
        """
        Register a realign group to track it
        
        Args:
            rg_label_or_index: Either an RG label ("RG0") or an index (0)
        """
        rg_label = self._normalize_rg_label(rg_label_or_index)
        self.realign_groups.add(rg_label)
        
        # Ensure the RG directory exists
        rg_dir = self.rg_dir(rg_label)
        os.makedirs(rg_dir, exist_ok=True)


    def _discover_and_register_realign_groups(self):
        """
        Discover existing realign groups from the file system and register them
        """
        if not os.path.exists(self.realign_groups_dir):
            return
            
        # Look for directories matching the RG pattern in the realign_groups directory
        for item in os.listdir(self.realign_groups_dir):
            item_path = os.path.join(self.realign_groups_dir, item)
            # If it's a directory and matches the RG pattern
            if os.path.isdir(item_path) and re.match(r'^RG\d+$', item):
                # Register this realign group
                self.register_realign_group(item)
                
        # Log what we found
        if self.realign_groups:
            print(f"Discovered {len(self.realign_groups)} existing realign groups: {', '.join(sorted(self.realign_groups))}", file=sys.stderr)

    
    def _normalize_rg_label(self, rg_label_or_index) -> str:
        """
        Convert various inputs to a standardized RG label format.
        
        Args:
            rg_label_or_index: Either an RG label ("RG0") or an index (0)
            
        Returns:
            Standardized RG label in the format "RG{index}"
            
        Raises:
            ValueError: If the input can't be converted to a valid RG label
        """
        # If it's an integer, format as "RG{index}"
        if isinstance(rg_label_or_index, int):
            return f"RG{rg_label_or_index}"
        
        # If it's a string
        if isinstance(rg_label_or_index, str):
            # If it's a numeric string, treat as index
            if rg_label_or_index.isdigit():
                return f"RG{rg_label_or_index}"
            
            # If it already matches the RG pattern, use directly
            if re.match(r'^RG\d+$', rg_label_or_index):
                return rg_label_or_index
        
        # Invalid input
        raise ValueError(f"Invalid RG label or index: {rg_label_or_index}. "
                         f"Must be an integer, numeric string, or string matching 'RG<digits>' pattern.")
    
    def get_all_realign_groups(self) -> Set[str]:
        """Get all registered realign groups"""
        return self.realign_groups
    
    def qnode_grouping_graph(self) -> str:
        """Get path for qnode grouping graph"""
        return os.path.join(self.work_dir, f"{self.basename}_qnode_grouping.graphml")
    
    # File path getter methods
    def multi_align_bed_path(self) -> str:
        """Path for multi-aligned regions BED file"""
        return os.path.join(self.work_dir, f"{self.basename}.{self.target_tag}.multialign.bed")
    
    def raw_sd_binary_map_path(self) -> str:
        """Path for raw SD binary map"""
        return os.path.join(self.work_dir, "raw_SD_binary_map.tsv")
    
    def filtered_sd_binary_map_path(self) -> str:
        """Path for filtered SD binary map"""
        return os.path.join(self.work_dir, "filtered_SD_binary_map.tsv")
    
    def multiplex_graph_path(self) -> str:
        """Get path for multiplex graph"""
        return os.path.join(self.work_dir, f"{self.basename}_multiplexed_SDs.graphml")
    
    def annotated_graph_path(self) -> str:
        """Get path for annotated graph"""
        return self.multiplex_graph_path().replace(".graphml", ".trim.annoPC.graphml")
    
    def directed_graph_path(self, chrom: str) -> str:
        """Get path for directed graph for a specific chromosome"""
        return self.multiplex_graph_path().replace(".graphml", f".directed.overlap.{chrom}.graphml")
    
    def target_overlapping_query_sd_bed_path(self) -> str:
        """Get path for target overlapping query SD BED"""
        return os.path.join(self.dirs["root"], f"{self.basename}target_overlapping_query_SD.bed")
    
    # RG (Realign Group) related paths and directories
    def rg_dir(self, rg_label_or_index) -> str:
        """Get directory for a specific RG"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        rg_dir = os.path.join(self.dirs["realign_groups"], rg_label)
        os.makedirs(rg_dir, exist_ok=True)
        return rg_dir

    def rg_raw_masked_bam_path(self, rg_label_or_index) -> str:
        """Path for RG raw masked bam file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.only_{rg_label}.raw.bam")

    def rg_realign_fastqs_path(self, rg_label_or_index) -> str:
        """Path for RG realign fastqs file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.only_{rg_label}.r1.fastq"), \
               os.path.join(self.recall_results_dir, f"{self.sample_id}.sdrecall.only_{rg_label}.r2.fastq")
    
    def rg_query_bed_path(self, rg_label_or_index) -> str:
        """Path for RG bed file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}.bed")
    
    def rg_counterparts_bed_path(self, rg_label_or_index) -> str:
        """Path for counterparts bed file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}_counterparts.bed")
    
    def all_homo_regions_bed_path(self, rg_label_or_index) -> str:
        """Path for all homologous regions bed file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}_related_homo_regions.bed")
    
    def all_homo_regions_fastq_path(self, rg_label_or_index) -> str:
        """Get path for all homologous regions fastq file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}_related_homo_regions.raw.fastq")
    
    def masked_genome_path(self, rg_label_or_index) -> str:
        """Get path for masked genome"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}.masked.fasta")
    
    def masked_genome_fai_path(self, rg_label_or_index) -> str:
        """Get path for masked genome FAI index"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return f"{self.masked_genome_path(rg_label)}.fai"
    
    def masked_genome_contigsize_path(self, rg_label_or_index) -> str:
        """Get path for masked genome contig size file"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}.masked.contigsize.genome")
    
    def minimap_index_path(self, rg_label_or_index) -> str:
        """Path for minimap2 index"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}.masked.mmi")
    
    def total_intrinsic_bam_path(self) -> str:
        """Path for total intrinsic alignment BAM"""
        return os.path.join(self.work_dir, f"total_intrinsic_alignments.bam")
    
    def raw_intrinsic_bam_path(self) -> str:
        """Path for raw intrinsic alignment BAM"""
        raw_intrin_bam = os.path.join(self.repo_dir, "data", self.assembly, "raw_intrin_align", f"assembly_intrinsic_align.WGAC.{self.assembly}.reformat.bam")
        if not os.path.exists(raw_intrin_bam):
            raise FileNotFoundError(f"Raw intrinsic alignment BAM file not found at {raw_intrin_bam}, please check the input assembly tag.")
        elif os.path.getsize(raw_intrin_bam) < 10000:
            raise FileExistsError(f"Raw intrinsic alignment BAM seems not properly downloaded at {raw_intrin_bam}, please check the git lfs status to see if the big BAM file is properly downloaded to the local file system.")
        else:
            return raw_intrin_bam
    
    def intrinsic_bam_path(self, rg_label_or_index) -> str:
        """Path for intrinsic alignment SAM"""
        rg_label = self._normalize_rg_label(rg_label_or_index)
        return os.path.join(self.rg_dir(rg_label), f"{rg_label}.intrinsic.bam")
    
    # Methods to get files across all realign groups
    def all_rg_query_bed_paths(self) -> List[str]:
        """Get paths to all RG bed files"""
        return [self.rg_query_bed_path(rg) for rg in self.realign_groups]
    
    def total_recall_SD_region_bed_path(self) -> str:
        """Get the bed file path containing all the realign target regions"""
        return os.path.join(self.realign_groups_dir, "all_target_recall_SD_regions.bed")
    
    def all_masked_genome_paths(self) -> List[str]:
        """Get paths to all masked genome files"""
        return [self.masked_genome_path(rg) for rg in self.realign_groups]
    
    def all_intrinsic_bam_paths(self) -> List[str]:
        """Get paths to all intrinsic alignment SAM files"""
        return [self.intrinsic_bam_path(rg) for rg in self.realign_groups]
    
    def all_rg_counterparts_bed_paths(self) -> List[str]:
        """Get paths to all counterparts bed files"""
        return [self.rg_counterparts_bed_path(rg) for rg in self.realign_groups]
    
    def all_homo_regions_bed_paths(self) -> List[str]:
        """Get paths to all homologous regions bed files"""
        return [self.all_homo_regions_bed_path(rg) for rg in self.realign_groups]
    
    def total_homo_regions_bed_path(self) -> str:
        """Get path to all homologous regions bed file"""
        return os.path.join(self.work_dir, "all_RG_related_homo_regions.bed")