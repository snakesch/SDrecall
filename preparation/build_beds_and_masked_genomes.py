import sys
from multiprocessing import Pool
from pybedtools import BedTool
from pybedtools import helpers

from src.utils import executeCmd, sortBed_and_merge, merge_bed_files
from src.const import shell_utils, SDrecallPaths
from src.log import error_handling_decorator, logger
from preparation.homoseq_region import HOMOSEQ_REGION
from preparation.genome import Genome
from preparation.intrinsic_alignment import getIntrinsicBam


def build_beds_and_masked_genomes(grouped_qnode_cnodes: list,
                                  sd_paralog_pairs: dict,
                                  sdrecall_paths: SDrecallPaths,
                                  nthreads=12,
                                  avg_frag_size=400,
                                  std_frag_size=140):    
    # Label SD-paralog pairs. Name disconnected qnodes as RG0, connected qnodes as PC1, PC2, ...
    # We need to restructure the grouped_qnode_cnodes to pass the information of sd-paralog pairs to the establish_beds_per_RG_cluster function
    helpers.set_tempdir(sdrecall_paths.tmp_dir)
    new_results = []
    for result in grouped_qnode_cnodes:
        new_result = {"SD_qnodes": {}, "SD_counterparts": {}}
        cluster_idx = 0
        for i in range(0, len(result["SD_qnodes"])):
            fc_node = result["SD_qnodes"][i]
            if fc_node in sd_paralog_pairs:
                new_result['SD_qnodes'][cluster_idx] = [fc_node]
                new_result['SD_counterparts'][cluster_idx] = sd_paralog_pairs[fc_node]
                cluster_idx += 1
        new_results.append(new_result)

    # Load balancing
    new_results = sorted(new_results, key = lambda x: sum(v[0][2] - v[0][1] for k,v in x["SD_qnodes"].items()) * sum(len(v) for k,v in x["SD_counterparts"].items()), reverse=True)
    labels = [ "RG" + str(n) for n in range(0, len(grouped_qnode_cnodes))]

    # Register all realign groups
    for label in labels:
        sdrecall_paths.register_realign_group(label)

    # Create arguments for parallel processing
    worker_args = []
    for i, (result, label) in enumerate(zip(new_results, labels)):
        # For each worker, we create a tuple of all needed arguments
        worker_args.append((
            result,
            label,
            (
                ("Query_bed", sdrecall_paths.rg_query_bed_path(label)),
                ("Counterparts_bed", sdrecall_paths.rg_counterparts_bed_path(label)),
                ("All_region_bed", sdrecall_paths.all_homo_regions_bed_path(label)),
                ("Masked_genome", sdrecall_paths.masked_genome_path(label)), 
                ("Intrinsic_bam", sdrecall_paths.intrinsic_bam_path(label)),
                ("Ref_genome", sdrecall_paths.ref_genome)
            ),
            avg_frag_size,
            std_frag_size,
            2  # threads per worker
        ))
    
    # Run in parallel
    pool = Pool(nthreads)
    results = pool.imap_unordered(imap_establish, worker_args)
    i = 0
    intrinsic_bams = []
    for success, result, logs in results:
        logger.debug(f"{i}th subprocess started")
        if not success:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
        else:
            intrinsic_bams.append(result)
        logger.debug(logs)
        logger.debug(f"{i}th subprocess completed")
        i+=1
    
    pool.close()
    if not all([t[0] for t in results]):
        raise RuntimeError("Error occurred during the parallel execution of establish_beds_per_RG_cluster. Exit.")

    ## total_intrinsic_alignments.bam is a merged alignment file for BILC model.
    intrinsic_bams = [ sdrecall_paths.intrinsic_bam_path(label) for label in labels ]
    total_intrinsic_bam = sdrecall_paths.total_intrinsic_bam_path()
    
    ## Create total_intrinsic_alignments.bam
    intrinsic_bam_header = total_intrinsic_bam.replace(".bam", ".bam.header")
    cmd = f"source {shell_utils} && modify_bam_sq_lines {intrinsic_bams[0]} {sdrecall_paths.ref_genome} {intrinsic_bam_header}"
    executeCmd(cmd, logger=logger)

    intrinsic_bam_list = total_intrinsic_bam.replace(".bam", ".bams.list.txt")
    with open(intrinsic_bam_list, "w") as f:
        f.write("\n".join(intrinsic_bams))

    cmd = f"samtools merge -@ {nthreads} -h {intrinsic_bam_header} -b {intrinsic_bam_list} -o - | \
            samtools sort -T {sdrecall_paths.tmp_dir} -O bam -o {total_intrinsic_bam} && \
            samtools index {total_intrinsic_bam} && \
            ls -lht {total_intrinsic_bam} || \
            echo Failed to concatenate all the filtered realigned BAM files."
    executeCmd(cmd)
               
    # Fourth, extract all PC*_related_homo_regions.bed file and concat them together.
    beds = sdrecall_paths.all_homo_regions_bed_paths()
    combined_bed = merge_bed_files(beds, tmp_dir=sdrecall_paths.tmp_dir)
    combined_bed.saveas(sdrecall_paths.total_homo_regions_bed_path())

    return sdrecall_paths

def imap_establish(tup_args):
    return establish_beds_per_RG_cluster(*tup_args)

@error_handling_decorator
def establish_beds_per_RG_cluster(cluster_dict={"SD_qnodes":{},
                                                "SD_counterparts":{}},
                                  label = "RG0",
                                  path_tups = (),
                                  avg_frag_size = 400,
                                  std_frag_size = 140,
                                  threads = 2,
                                  logger = logger):

    """
    input cluster dict now looks like this:
    {
        "SD_qnodes": {0: [], 1: [], 2: []}, 
        "SD_counterparts": {0: [], 1: [], 2: []}
    }
    """

    ## Initialize file structure
    paths = {k:v for k,v in path_tups}
    
    # Create the query region bed file
    with open(paths["Query_bed"], "w") as f:
        for idx, records in cluster_dict["SD_qnodes"].items():
            for record in records:
                if len(record) >= 3:
                    # Invoke __iter__ method instead of __getitem__
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3]]) + "\n")
                elif len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3]]) + "\n")
    sortBed_and_merge(paths["Query_bed"])
            
    ## Create the counterparts region bed file
    with open(paths["Counterparts_bed"], "w") as f:
        for idx, records in cluster_dict["SD_counterparts"].items():
            fc_node = cluster_dict["SD_qnodes"][idx][0]
            for record in records:
                assert isinstance(record, HOMOSEQ_REGION), "The record in the SD_counterparts list is not a HOMOSEQ_REGION object: {}".format(record)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)     
                ## Invoke __iter__ method instead of __getitem__
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3]]) + "\n")

    ## Remove PC regions from counterpart BEDs and force strandedness
    BedTool(paths["Counterparts_bed"]).subtract(BedTool(paths["Query_bed"]), s=True).saveas(paths["Counterparts_bed"])
    sortBed_and_merge(paths["Counterparts_bed"])
    
    ## Create the total region bed file by combining PC BED and counterpart BED
    ## Invoke __iter__ method instead of __getitem__
    with open(paths["All_region_bed"], "w") as f:
        for idx, records in cluster_dict["SD_qnodes"].items():
            for record in records:
                if len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3], f"FC:{label}_{idx}"]) + "\n")
                elif len(record) >= 3:
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3], f"FC:{label}_{idx}"]) + "\n")
        for idx, records in cluster_dict["SD_counterparts"].items():
            fc_node = cluster_dict["SD_qnodes"][idx][0]
            for record in records:
                if not isinstance(record, HOMOSEQ_REGION):
                    logger.critical("Some elements of the SD_counterparts list are not HOMOSEQ_REGION objects: {}".format(record))
                    sys.exit(1)
                fc_node_rela_start, fc_node_rela_end = record.qnode_relative_region(fc_node)
                # The rela_start and rela_end are based on the corresponding PC region
                f.write("\t".join([str(value) for value in record][:3] + [str(fc_node_rela_start), str(fc_node_rela_end), record[3], f"NFC:{label}_{idx}"]) + "\n")

    # Now all the bed files are finalized for this RG cluster
    # We can proceed to prepare the masked genome and intrinsic bam
    ref_genome = paths["Ref_genome"]
    masked_genome_path = paths["Masked_genome"]
    contig_sizes = ref_genome.replace(".fasta", ".fasta.fai")
    # Prepare masked genomes
    masked_genome = Genome(ref_genome).mask(paths["Query_bed"], 
                                            avg_frag_size = avg_frag_size, 
                                            std_frag_size = std_frag_size, 
                                            genome=contig_sizes, 
                                            logger=logger, 
                                            path=masked_genome_path)
    
    # Call intrinsic variants
    bam_path = getIntrinsicBam( rg_bed = paths["Query_bed"], 
                                all_homo_regions_bed = paths["All_region_bed"], 
                                rg_masked = masked_genome,
                                ref_genome = ref_genome,
                                intrinsic_bam = paths["Intrinsic_bam"],
                                avg_frag_size = avg_frag_size,
                                std_frag_size = std_frag_size,
                                threads = threads )
    return bam_path