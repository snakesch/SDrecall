import networkx as nx
from intervaltree import Interval, IntervalTree
from functools import partial
from collections import defaultdict
from python_utils import convert_input_value
from extract_SD_pairs_from_bam_json import main_parse_json_and_process
from pick_multialign_region import main_func_pick_region
from community import community_louvain
from HashableDict import HashableDict
import ast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

bash_utils_hub = "/paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh"

class HashableGraph:
    def __init__(self, graph):
        self.graph = graph

    def __hash__(self):
        node_data = sorted(self.graph.nodes(data=True))
        edge_data = sorted(self.graph.edges(data=True))
        return hash((tuple(node_data), tuple(edge_data)))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __getattr__(self, name):
        # Delegate attribute access to the internal graph object
        return getattr(self.graph, name)

    def to_graph(self):
        return self.graph


def calculate_bed_size(bed_file):
    if isinstance(bed_file, str):
        bed_file = BedTool(bed_file)
    total_size = 0
    for feature in bed_file.sort().merge():
        total_size += len(feature)
    return total_size


def prepare_tmp_file(tmp_dir="/paedyl01/disk1/yangyxt/test_tmp", **kwargs):
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass

    return tempfile.NamedTemporaryFile(dir = "/paedyl01/disk1/yangyxt/test_tmp", delete = False, **kwargs)
class Genome:
    """
    A class implementation for genome masking.

    Usage:
    rg = Genome(path)
    rg.mask(bedf) # bedf needs not to be complemented

    """

    def __init__(self, path):
        self.path = path
        self.faiIndex = self._fai()
        self.bwtIndex = self._bwtIndex()
        self.total_bed = self._total_bed()

    def _fai(self):

        import os

        if not os.path.exists(os.path.join(self.path + ".fai")) or (os.path.getmtime(self.path) > os.path.getmtime(self.path + ".fai")):
            executeCmd(f"samtools faidx {self.path} ")
        return os.path.join(self.path + ".fai")


    def _bwtIndex(self):

        import os

        if not os.path.exists(os.path.join(self.path + ".bwt")) or (os.path.getmtime(self.path) > os.path.getmtime(self.path + ".bwt")):
            executeCmd(f"bwa index {self.path}")
        return os.path.join(self.path + ".bwt")


    def _total_bed(self):

        import os

        total_bed_path = ".".join(self.path.split(".")[:-1]) + ".bed"
        return total_bed_path if os.path.exists(total_bed_path) else ""


    def getContigGenome(self):

        import pandas as pd

        contig_genome_fn = ".".join(self.path.split(".")[:-1]) + ".contigsize.genome"
        df = pd.read_csv(self.faiIndex(), sep="\t", header=None).iloc[:, [0,1]]
        df.to_csv(contig_genome_fn, sep="\t", index=False, header=False)


    def mask(self, bedf: str, avg_frag_size=400, std_frag_size = 130, genome="hg19"):
        '''
        This function actually reverse mask the input bed regions in the genome instead of mask the bed regions
        '''

        import os
        import math
        import pandas as pd

        if not os.path.exists(bedf):
            raise FileNotFoundError(f"Invalid BED file : {bedf} ")

        if not os.path.exists(self.total_bed):
            raise FileNotFoundError(f"Invalid BED file : {self.total_bed} ")

        _base_fn = os.path.basename(bedf)
        if "_related_homo_regions" in _base_fn:
            region = _base_fn.split("_")[0]
        else:
            region = _base_fn.split(".")[0] # We assume the format is REGION.bed

        # Get target contigs
        target_bed = BedTool(bedf).sort().merge(d=int(2*(avg_frag_size + 2*std_frag_size))).slop(b=int(avg_frag_size + 2*std_frag_size + 1000), g=genome).sort()
        # Instead of using seqtk, we might need to use biopython here
        ref_genome_seq = SeqIO.to_dict(list(SeqIO.parse(self.path, "fasta")))
        # Every interval needs to be an independent contig
        # BedTool intervals are semi-open 0-indexed intervals, which works in the same way as python string slicing
        masked_genome_contigs = []
        for interval in target_bed:
            interval_seq = ref_genome_seq[interval.chrom].seq[interval.start:interval.stop]
            mutable_seq = str(interval_seq)
            mask_seq = "".join(["N" for i in range(0,1000)])
            mutable_seq = mask_seq + mutable_seq[1000:-1000] + mask_seq
            try:
                assert str(mutable_seq) != str(interval_seq)
            except AssertionError as ae:
                logger.warning(f"Sequences from {interval.chrom}:{interval.start}-{interval.stop} in {self.path} seems already have Ns at both ends.")
            assert len(mutable_seq) == len(str(interval_seq))
            assert mutable_seq.count("N") >= 2000
            interval_seq = Seq(mutable_seq)
            interval_contig = SeqRecord(interval_seq, id=f"{interval.chrom}:{interval.start}", description="")
            masked_genome_contigs.append(interval_contig)

        # Now we need to concat the fasta sequences on all chromosomes into one genome and write it into a fasta file
        tmp_tag = str(uuid.uuid4())
        masked_genome = os.path.join(os.path.dirname(bedf), os.path.basename(".".join(self.path.split(".")[:-1]) + ".sub." + region + ".masked.fasta"))
        logger.info("The masked genome will be written to {}".format(masked_genome))
        tmp_genome = masked_genome.replace(".fasta", ".{}.fasta".format(tmp_tag))
        SeqIO.write(masked_genome_contigs, tmp_genome, "fasta")

        updated = update_plain_file_on_md5(masked_genome, tmp_genome)
        Genome(masked_genome)  # Will run building fai index and dict index and bwt index in the __init__ function

        if os.path.exists(self.path + ".seqkit.fai"):
            os.remove(self.path + ".seqkit.fai")

        return masked_genome


# Define a helper function to explode the rows
def explode_row(row, columns):
    other_value_scalar = max(*[len(row[c]) for c in columns])
    other_columns = [c for c in row.index.tolist() if c not in columns]

    if len(set([len(row[c]) for c in columns])) > 1:
        raise ValueError("This record does not have equal length of lists in these input columns: {}. Take a look at the record:\n{}\n".format(columns, row))

    new_rows = {c: [row[c]]*other_value_scalar for c in other_columns}
    new_rows.update({c: row[c] for c in columns})

    return pd.DataFrame(new_rows)


def makeWindows(df, length, logger=logger):
    """
    This function takes a BED-like dataframe and makes windows of size LENGTH.

    Regions with length of trailing block < 100 are appended to the previous block.

    """
    start = df.iloc[:, 1].tolist()
    end = df.iloc[:, 2].tolist()
    name = df.iloc[:, 3].tolist()
    score = df.iloc[:, 4].tolist()
    strand = df.iloc[:, 5].tolist()
    ret_start, ret_end, ret_name, ret_score, ret_strand = [], [], [], [], []
    for i in range(len(start)):
        _start, _end, _name, _score, _strand = [], [], [], [], []
        cur_start, cur_end = start[i], end[i]
        short = True
        while cur_end - cur_start > length:
            _start.append(cur_start)
            _end.append(cur_start + length)
            _name.append(name[i])
            _score.append(score[i])
            _strand.append(strand[i])
            cur_start += length
            short = False
        if not short and cur_end - cur_start < 100:
            _end[-1] = cur_end
        else:
            _start.append(cur_start)
            _end.append(cur_end)
            _name.append(name[i])
            _score.append(score[i])
            _strand.append(strand[i])
        ret_start.append(_start)
        ret_end.append(_end)
        ret_name.append(_name)
        ret_score.append(_score)
        ret_strand.append(_strand)

    out_df = pd.DataFrame({
        "Chr": df.iloc[0,0],
        "Start": ret_start,
        "End": ret_end,
        "Name": ret_name,
        "Score": ret_score,
        "Strand": ret_strand
    })

    try:
        return pd.concat(out_df.dropna().apply(explode_row, axis=1, args=(["Start", "End", "Name", "Score", "Strand"],)).tolist(), ignore_index=True)
    except ValueError as ve:
        logger.error("Running into ValueError, let's take a look at the input groupdf file first: \n{}\nThen look at the decomposed groupdf:\n{}\n".format(df.to_string(index=False),
                                                                                                                                                        out_df.to_string(index=False)))
        raise ve


def split_fastq(input_fastq, output_fastq1, output_fastq2, length=150):
    with open(input_fastq, "r") as input_handle, \
         open(output_fastq1, "w") as output_handle1, \
         open(output_fastq2, "w") as output_handle2:

        for record in SeqIO.parse(input_handle, "fastq"):
            # Extract first 150 bases
            record1 = record[:length]
            # Extract first 150 bases of the reverse complement
            record2 = record.reverse_complement()[:length]
            record2.id = record.id
            record2.description = record.description

            SeqIO.write(record1, output_handle1, "fastq")
            SeqIO.write(record2, output_handle2, "fastq")


def getSubseq(bedf_path: str, fastq_path: str, ref_genome: str, length = None, logger=logger) -> None:
    # Here the input bed file should be the bed file recording the entire SD set involved in one particular PC
    if length:
        if length <= 100:
            raise ValueError(f"Read length {length} is too small (<100). ")
    length = length + 1 if length % 2 == 1 else length
    bedf = pd.read_csv(bedf_path, sep="\t", header=None)

    if length:
        window_df = bedf.groupby(by=0, group_keys=False).apply(lambda x: makeWindows(x, length, logger=logger)).reset_index(level=0, drop=True)

        window_df["Length"] = window_df.iloc[:, 2] - window_df.iloc[:, 1]
        bedf[3] = bedf[2] - bedf[1]
        if window_df.Length.sum() != bedf[3].sum():
            raise RuntimeError(f"Error in preprocessing {bedf_path}. Lengths of regions {window_df.Length.sum()} do not match. {bedf[3].sum()} ")
    else:
        window_df = bedf
    window_out = bedf_path + ".preproc"
    window_df.to_csv(window_out, sep="\t", index=False, header=False)
    cmd = f"if [[ -d {fastq_path} ]]; then rm -R {fastq_path}; fi"
    executeCmd(cmd, logger=logger)

    cmd = f"seqtk subseq -s {ref_genome} {window_out} | seqtk seq -F 'F' - | mawk 'FNR %4 == 2{{ print toupper($0); }} FNR % 4 !=2 {{ print; }}' > {fastq_path} "
    executeCmd(cmd, logger=logger)
    logger.info(f"Subsequence written to {fastq_path}. ")

    split_fastq(fastq_path, fastq_path.replace(".fastq", ".1.fastq"), fastq_path.replace(".fastq", ".2.fastq"), length=int(length/2))
    os.remove(window_out)

    return (fastq_path.replace(".fastq", ".1.fastq"), fastq_path.replace(".fastq", ".2.fastq"))


def getIntrinsicVcf(pc_bed,
                    all_homo_regions_bed,
                    counter_bed,
                    pc_masked,
                    ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                    avg_frag_size = 400,
                    std_frag_size = 120,
                    threads = 2,
                    logger = logger):
    # Here fp stands for file_path
    # This function requires an input of:
    # 1. masked genome
    # 2. BED file of PC region
    # 3. BED file of all homologous regions related to the PC region
    # 4. Reference genomef
    # 5. Length stands for the read length

    if (not os.path.exists(pc_masked + ".bwt")) or \
        (os.path.getmtime(pc_masked + ".bwt") < os.path.getmtime(pc_masked)):
        cmd = f"bwa index {pc_masked} "
        executeCmd(cmd, logger=logger)
    else:
        logger.info(f"PC reverse-masked genome index detected - {pc_masked + '.bwt'}")

    length = int(avg_frag_size)
    homo_sd = os.path.dirname(all_homo_regions_bed)
    fastq_dir = homo_sd
    bam_dir = os.path.dirname(pc_bed)
    vcf_dir = os.path.dirname(pc_bed)
    fq_path = os.path.join(fastq_dir, os.path.basename(all_homo_regions_bed)[:-3] + str(length) + ".fastq")
    bam_path = os.path.join(bam_dir, os.path.basename(pc_bed)[:-3] + str(length) + ".bam")
    vcf_path = bam_path.replace(".bam", ".vcf.gz")
    pc_label = os.path.basename(pc_bed.replace(".bed",""))

    # Get Fastq files, note that these reference genome sequences are extracted to single-end sequences intead of paired end sequences.
    # First predetermine the path of the fastq file
    fq_path1, fq_path2 = getSubseq(all_homo_regions_bed, fq_path, ref_genome, length, logger=logger)
    logger.info(f"The paired fastq files to get the intrinsic VCF for {pc_bed} are {fq_path1} and {fq_path2}")

    # After getting the fastq file it should be used to map against the masked genome and perform straight forward variants calling
    if not os.path.exists(bam_path) or (os.path.getmtime(bam_path) < os.path.getmtime(pc_masked)):
        pc_masked_index = pc_masked.replace(".fasta", ".mmi")
        cmd = f"bash {bash_utils_hub} independent_minimap2_masked -f {fq_path1} -r {fq_path2} -a {pc_masked} -o {bam_path} -s all_PC -t {threads} -i {pc_label} && ls -lh {bam_path}*"
        executeCmd(cmd, logger=logger)

    os.makedirs(vcf_dir, exist_ok=True)
    tmp_vcf = vcf_path.replace(".vcf", ".tmp.vcf")

    cmd = f"export OPENBLAS_NUM_THREADS={threads} && bcftools mpileup -a FORMAT/AD,FORMAT/DP -f {ref_genome} {bam_path} | " + \
          f"bcftools call -mv -P 0 -Oz -o {tmp_vcf} && tabix -p vcf {tmp_vcf} "
    executeCmd(cmd, logger=logger)
    cmd = f"export OPENBLAS_NUM_THREADS={threads} && bcftools norm -m -both -f {ref_genome} --multi-overlaps 0 --keep-sum AD -a {tmp_vcf} | \
            bcftools norm -d exact - | \
            bcftools view -i 'ALT!=\"*\"' | \
            bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp -Oz -o {vcf_path} && tabix -f -p vcf {vcf_path} && rm {tmp_vcf} {fq_path}"
    executeCmd(cmd, logger=logger)

    logger.info(f"Writing intrinsic VCF to {vcf_path}")
    return bam_path




def executeCmd(cmd, stdout_only = False, logger = logger) -> None:

    import subprocess
    from subprocess import PIPE

    if stdout_only:
        result = subprocess.run(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    else:
        result = subprocess.run(cmd, shell=True, stderr=subprocess.STDOUT, stdout=PIPE)

    logger.info(f"Running the following shell command inside python:\n{cmd}\nAnd the output goes like this:\n{result.stdout.decode()}\n\n")
    code = result.returncode
    cmd_lst = cmd.split(" ")
    if code != 0:
        if cmd_lst[1][0] != "-":
            raise RuntimeError("Error in " + " ".join(cmd_lst))
        else:
            raise RuntimeError("Error in {}".format(cmd_lst))

    return result.stdout.decode()


def create_community_graph(graph, partitions):
    community_graph = nx.Graph()
    communities = { community: [k for k,v in partitions.items() if v == community] for community in set(partitions.values()) }
    overlap_graph_nodes = set(partitions.keys())
    # Now, we add the community nodes to the graph, along with their attributes\
    total_overlap_avgs = []
    for community_index, nodes in communities.items():
        chr_community = min(nodes, key=lambda x: x[0])[0]  # The chromosome is the smallest chromosome in lexicographical order
        start_community = min(nodes, key=lambda x: x[1])[1]  # The start is the smallest start value
        end_community = max(nodes, key=lambda x: x[2])[2]  # The end is the largest end value

        nodes_overlap_fractions = [(n[2] - n[1])/(end_community - start_community) for n in nodes]
        total_overlap_avgs.append(np.mean(nodes_overlap_fractions))

        tmp_graph = prepare_tmp_file().name
        nx.write_gpickle(graph.subgraph(nodes).copy(), tmp_graph)
        community_graph.add_node(community_index, chr = chr_community, start = start_community, end = end_community, raw_nodes_graph = tmp_graph )

    logger.info("The total overlap fraction averages of the nodes in communities are:{}".format(np.mean(total_overlap_avgs)))

    # Add edges between communities if there is an edge between their nodes in the original graph
    non_overlapping_nodes = set([])
    for u, v, data in graph.edges(data=True):
        if u not in overlap_graph_nodes:
            # This mean u is not highly overlapped with other nodes so it is not included in the overlap graph
            non_overlapping_nodes.add(u)
            continue
        if v not in overlap_graph_nodes:
            # This mean v is not highly overlapped with other nodes so it is not included in the overlap graph
            non_overlapping_nodes.add(v)
            continue
        if partitions[u] != partitions[v]:  # u and v belong to different communities
            if not community_graph.has_edge(partitions[u], partitions[v]):  # there is no edge between the communities of u and v
                community_graph.add_edge(partitions[u], partitions[v], type="segmental_duplication")

    reformatted_graph = nx.Graph()

    for node, data in community_graph.nodes(data=True):
        new_node = (data['chr'], data['start'], data['end'])
        reformatted_graph.add_node(new_node, size=int(data['end']) - int(data['start']), **data)

    for u, v, data in community_graph.edges(data=True):
        new_u = (community_graph.nodes[u]['chr'], community_graph.nodes[u]['start'], community_graph.nodes[u]['end'])
        new_v = (community_graph.nodes[v]['chr'], community_graph.nodes[v]['start'], community_graph.nodes[v]['end'])
        reformatted_graph.add_edge(new_u, new_v, **data)

    for node in non_overlapping_nodes:
        reformatted_graph.add_node(node, size = int(node[2]) - int(node[1]))

    return reformatted_graph


def pickout_bridging_nodes(dgraph, frag_size=400, std_frag_size = 120, min_out_overlap_frac = 0.8, max_in_overlap_frac = 0.3):
    bridge_node_candidates = []
    # Return a list of nodes that are bridging nodes
    for node, data in dgraph.nodes(data=True):
        size = data["size"]
        if size < frag_size * 1.5:
            continue
        out_edges = list(dgraph.out_edges(node, data=True))
        in_edges = list(dgraph.in_edges(node, data=True))
        out_edges_ws = sum([e[2]["weight"] for e in out_edges])
        in_edges_ws = sum([e[2]["weight"] for e in in_edges])
        max_in_edge_weight = max([e[2]["weight"] for e in in_edges]) if len(in_edges) > 0 else 0
        if out_edges_ws > 2*min_out_overlap_frac:
            bridge_node_candidates.append(node)

    final_bridge_nodes = []
    for node in bridge_node_candidates:
        out_edges = list(dgraph.out_edges(node, data=True))
        smaller_intervals = [e[1] for e in out_edges]
        # Convert the smaller intervals to bedtool object and perform merge and sort
        combined_small_interval_region = BedTool("\n".join(["\t".join([str(x) for x in interval]) for interval in smaller_intervals]), from_string=True).sort().merge()
        node_bed = BedTool("\t".join([str(x) for x in node]), from_string=True)

        # Calculate total overlap with the big interval
        total_overlap_bed = node_bed.intersect(combined_small_interval_region, wo=True)
        if total_overlap_bed.count() == 0:
            continue
        total_overlap = total_overlap_bed.to_dataframe( disable_auto_names=True,
                                                        names=["chr_1", "start_1", "end_1",
                                                            "chr_2", "start_2", "end_2"]).iloc[:, -1].sum()
        overlap_frac = total_overlap / (node[2] - node[1])
        dgraph.nodes[node]["overlapped_region_size"] = overlap_frac * (node[2] - node[1])
        dgraph.nodes[node]["overlapped_region_frac"] = overlap_frac
        if overlap_frac >= 0.8 or (1-overlap_frac) * (node[2] - node[1]) <= (frag_size - .675 * std_frag_size):
            final_bridge_nodes.append(node)
            dgraph.nodes[node]["umbrella"] = "True"

    return final_bridge_nodes, dgraph



def compose_multiplex_graph_per_chr(args):
    chr_id, sd_data_chr, overlap_frac_min, avg_frag_size, std_frag_size, resolution, graph_file_path = args
    sd_interval_tree = IntervalTree()
    directed_overlap_graph = nx.DiGraph()
    directed_graph_path = graph_file_path.replace(".graphml", f".directed.overlap.{chr_id}.graphml")
    assert directed_graph_path != graph_file_path, "The directed graph path is the same as the graph path, please check the input graph path"

    # Only compose overlap graph per chromosome
    for i in range(0, len(sd_data_chr)):
        row = sd_data_chr.iloc[i, :]
        directed_overlap_graph.add_node((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]), size = row["end_1"] - row["start_1"])
        directed_overlap_graph.add_node((row["chr_2"], row["start_2"], row["end_2"], row["strand2"]), size = row["end_2"] - row["start_2"])
        interval_1 = Interval(row["start_1"], row["end_1"], (row["chr_1"], row["start_1"], row["end_1"], row["strand1"]))
        interval_2 = Interval(row["start_2"], row["end_2"], (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]))
        for interval in [interval_1, interval_2]:
            if chr_id == interval.data[0]: sd_interval_tree.add(interval)
            for overlap in sd_interval_tree.overlap(interval.begin, interval.end):
                overlap_span = min([interval.end, overlap.end]) - max([interval.begin, overlap.begin])
                overlap_size = overlap.end - overlap.begin
                interval_size = interval.end - interval.begin
                weight = max(overlap_span/overlap_size, overlap_span/interval_size)
                uncovered_size = min(overlap_size, interval_size) - overlap_span
                if interval_size >= overlap_size and interval.data[0] == overlap.data[0] and overlap.data != interval.data:
                    directed_overlap_graph.add_edge(interval.data, overlap.data, weight=weight, overlap="True", uncovered_smaller_interval=uncovered_size)
                elif interval.data[0] == overlap.data[0] and overlap.data != interval.data:
                    directed_overlap_graph.add_edge(overlap.data, interval.data, weight=weight, overlap="True", uncovered_smaller_interval=uncovered_size)

    # Now we would like to save the directed graph as a graphml file and see the results in cytoscape
    nx.write_graphml(directed_overlap_graph, directed_graph_path)
    gc.collect()
    return (chr_id, directed_graph_path)


def summary_exon_overlap(groupdf):
    # For directed graph
    # This function is to summarize the exon overlap information for each node in the directed graph
    groupdf["FC_overlap"] = groupdf["symbol"].astype(str) + ":" + groupdf["tranxID"].astype(str) + ":" + groupdf["exon_No"].astype(str) + ":" + groupdf["strand"].astype(str) + ":" + groupdf["overlap_len"].astype(str) + "bases"
    summarized_overlap = ",".join(groupdf["FC_overlap"].drop_duplicates().tolist())
    groupdf["FC_overlap"] = summarized_overlap
    return groupdf.loc[:, ["chr_graph", "start_graph", "end_graph", "graph_strand", "FC_overlap"]].drop_duplicates()


def create_multiplex_graph(sd_data, graph_filepath=None,
                            overlap_frac_min = 0.8,
                            avg_frag_size=350,
                            std_frag_size = 120,
                            resolution = .1,
                            threads = 10,
                            target_bed = "",
                            base_dir = ""):

    # sd_data is the binary sd map of genomic intervals
    # Input the sd_data to networkx Graph object
    G = nx.Graph()
    sd_interval_tree = {}
    all_chrs = set(sd_data["chr_1"]).union(set(sd_data["chr_2"]))

    if os.path.exists(target_bed):
        target_chrs = pd.read_table(target_bed, header=None).iloc[:, 0].drop_duplicates().tolist()
        all_chrs = [c for c in all_chrs if c in target_chrs]

    for chr_id in all_chrs:
        sd_interval_tree[chr_id] = IntervalTree()

    # adding edges between nodes in the same line
    # The below for loop not just build up a Graph composed of intervals sharing sequence similarity with each other
    # but also builds up an interval tree for each chromosome using the same intervals
    sd_data_by_chr = {chr_id: sd_data.loc[(sd_data["chr_1"]==chr_id) | (sd_data["chr_2"]==chr_id), :] for chr_id in all_chrs}
    for i in range(0, len(sd_data)):
        row = sd_data.iloc[i, :]
        # First adding the binary pairwise SDs to the graph
        G.add_edge((row["chr_1"], row["start_1"], row["end_1"], row["strand1"]),
                   (row["chr_2"], row["start_2"], row["end_2"], row["strand2"]),
                    type = "segmental_duplication")
    logger.info(f"How many nodes in the graph: {G.number_of_nodes()} How many edges are SD edges: {G.number_of_edges()} ")

    with Pool(threads) as pool:
        # Compose overlap graph within each chromsome
        partitions_by_chr = dict(pool.imap_unordered(compose_multiplex_graph_per_chr,
                                                    ((chr_id, sd_data_by_chr[chr_id], overlap_frac_min, avg_frag_size, std_frag_size, resolution, graph_filepath) for chr_id in all_chrs)))

    directed_graphs = [read_graphml(dg_path) for chr_id, dg_path in partitions_by_chr.items()]  # k is node and v is the directed graph path
    # Merge the directed graph into one graph
    merged_directed_graph = nx.DiGraph()
    for dg in directed_graphs:
        merged_directed_graph = nx.compose(merged_directed_graph, dg)

    # Then we add edges to merged_directed_graph based on the edges in the sd_graph networks which have their type attributes set to "segmental_duplication"
    for u, v, d in G.edges(data=True):
        if d.get("type", None) == "segmental_duplication":
            # Test if u and v are in the directed graph nodes
            if u in merged_directed_graph.nodes() and v in merged_directed_graph.nodes():
                # Test if these two nodes already has an edge between them
                if not merged_directed_graph.has_edge(u, v):
                    merged_directed_graph.add_edge(u, v, **d)
                else:
                    # Just add some attributes to the existing edge instead of force cover the original edge's attributes
                    merged_directed_graph[u][v]["type"] = "segmental_duplication"
            else:
                logger.warning("These two genomic regions sharing sequences seem not contained by the overlap graph: \n{}\n{}\n".format(u, v))

    # Now we add some annotations to the directed graph regarding which transcript is overlapped with each node.
    # First convert the graph to bed
    dg_bed = graph_to_sorted_bedtool(merged_directed_graph)
    target_bed_obj = BedTool(target_bed)
    dg_exon_df = dg_bed.intersect(target_bed_obj, wo = True, F=0.3, f=0.3, e=True).to_dataframe(disable_auto_names = True,
                                                                                                names =["chr_graph", "start_graph", "end_graph", "graph_strand",
                                                                                                        "chr_feature", "start_feature", "end_feature",
                                                                                                        "symbol", "tranxID",
                                                                                                        "exon_No", "strand",
                                                                                                        "overlap_len"])
    logger.info("Take a look at the overlap dataframe here:\n{}\n".format(dg_exon_df[:5].to_string(index=False)))
    dg_node_overlap_dict = {}
    by_graph_interval = dg_exon_df.groupby(by=["chr_graph", "start_graph", "end_graph", "graph_strand"], as_index=False)
    graph_groups = [group for _, group in by_graph_interval]
    with Pool(threads) as pool:
        results = pool.imap_unordered(summary_exon_overlap, graph_groups)
        summary_df = pd.concat(results, ignore_index=True).drop_duplicates()
    logger.info(f"The summary exon overlap function applied by groupby returned an object of {type(summary_df)}, take a look at the df:\n{summary_df[:5].to_string()}\n")
    for ind, row in summary_df.iterrows():
        node = (row["chr_graph"], int(row["start_graph"]), int(row["end_graph"]), row["graph_strand"])
        dg_node_overlap_dict[node] = row[-1]

    for node, overlap_descript in dg_node_overlap_dict.items():
        if node in merged_directed_graph.nodes():
            merged_directed_graph.nodes[node]["FC_overlap"] = overlap_descript

    sd_edge_count = [e.get("type", None) for u,v,e in merged_directed_graph.edges(data=True)].count("segmental_duplication")
    logger.info(f"How many nodes are there in the total SD graph with edges implying overlaps and homology? {merged_directed_graph.number_of_nodes()} How many edges in this graph: {merged_directed_graph.number_of_edges()}. How many of them belonged to type segmental_duplication: {sd_edge_count} ")

    if graph_filepath is not None:
        nx.write_graphml(merged_directed_graph, graph_filepath)
        logger.info("The merged directed graph is now saved as {}".format(graph_filepath))

    return merged_directed_graph




def merge_community_nodes(graph, nodes_to_merge):
    """
    Merge a list of nodes into a single node in the graph.

    Args:
    graph (networkx.Graph): The graph from where nodes should be merged.
    nodes_to_merge (list): The list of nodes to be merged.

    Returns:
    The updated graph.
    """
    new_graph = graph.copy()

    # Calculate the merged node attributes based on the nodes to merge.
    chr_values = set(node[0] for node in nodes_to_merge)
    assert len(chr_values) == 1, "All nodes to merge must have the same chromosome."

    new_node = (list(chr_values)[0], min(node[1] for node in nodes_to_merge), max(node[2] for node in nodes_to_merge))
    new_graph.add_node(new_node)

    edges_to_remove = []
    edges_to_add = []

    # For each node, redirect all its edges to the new_node
    for node in nodes_to_merge:
        neighbors = list(new_graph.neighbors(node))

        for neighbor in neighbors:
            # Move the edge to the new node, while keeping edge attributes.
            edge_data = new_graph.get_edge_data(node, neighbor)
            edges_to_remove.append((node, neighbor))
            edges_to_add.append((new_node, neighbor, edge_data))

        # Now, we can safely remove the node
        new_graph.remove_node(node)

    new_graph.remove_edges_from(edges_to_remove)
    new_graph.add_edges_from(edges_to_add)

    return new_graph



def pick_from_each_group(groups):
    '''
    Each PC should only contain the sequences in the same chromosome to reduce the size of the masked genome
    '''
    # First extract all the chromosomes in the groups
    # Each group in the groups is a list of tuples. Each tuple represents a query node
    transposed = list(itertools.zip_longest(*groups))
    # remove None values from each sublist (which are fill values for exhausted groups)
    result = list(dict.fromkeys([tuple([ e for e in sublist if e not in [np.nan, None]]) for sublist in transposed]))
    return result



def calculate_node_size(node):
    return int(float(node[2]) - float(node[1]))


def convert_node_into_bed_obj(node):
    return BedTool("\t".join([str(n) for n in node]), from_string=True)


def pick_sd_counterparts(sd_counterparts, query_node, sub_graph, fragment_size = 500):
    # First pickout the potential genomic intervals bridging not similar sequences together in the graph
    # The potential bridging genomic intervals should be the ones that are significantly larger (deprecated)
    # Or we can further detect community in the SD graph and extract the nodes belonged to the same community with the query node
    # First we assign each edge a weight (weight stands for the ratio of two interval size)
    for u, v, e in sub_graph.edges(data=True):
        try:
            usize = calculate_node_size(u)
            vsize = calculate_node_size(v)
        except KeyError:
            raise ValueError("This node pair, edge combo seems erroneous: node1: {}, node2:{}, edge:{}".format(u, v, e))
        else:
            e["weight"] = min(usize, vsize)/max(usize, vsize)

    partitions = community_louvain.best_partition(sub_graph, weight="weight")
    query_node_community = partitions[query_node]

    components = set([n for n in sub_graph.nodes() if partitions[n] == query_node_community]).difference({query_node,})
    return components




def bin_on_ploidy_fold_change(list_of_nodes,
                              node_ploidy_map):
    '''
    The input list of nodes are [(node1, node2, node3), (node1, node2) ... ]
    '''

    updated_list_of_nodes = []
    for group in list_of_nodes:
        # Separate the nodes into groups based on the ploidies of the region
        low_ploidy_group = [ n for n in group if node_ploidy_map[n] <= 3 ]
        mid_ploidy_group = [ n for n in group if node_ploidy_map[n] > 3 and node_ploidy_map[n] <= 5 ]
        high_ploidy_group = [ n for n in group if node_ploidy_map[n] > 5 and node_ploidy_map[n] <= 8 ]
        ultra_high_ploidy_group = [ n for n in group if node_ploidy_map[n] > 8 and node_ploidy_map[n] <= 10 ]
        extreme_ploidy_group = [ n for n in group if node_ploidy_map[n] > 10 ]

        # Display the number of nodes in each group
        logger.info("There are {} nodes in the low ploidy group, {} nodes in the mid ploidy group, {} nodes in the high ploidy group, {} nodes in the ultra high ploidy group, {} nodes in the extreme ploidy group".format(len(low_ploidy_group),
                                                                                                                                                                                                                            len(mid_ploidy_group),
                                                                                                                                                                                                                            len(high_ploidy_group),
                                                                                                                                                                                                                            len(ultra_high_ploidy_group),
                                                                                                                                                                                                                            len(extreme_ploidy_group)))

        groups = [g for g in [low_ploidy_group, mid_ploidy_group, high_ploidy_group, ultra_high_ploidy_group, extreme_ploidy_group] if len(g) > 0]
        updated_list_of_nodes.extend(groups)
    return updated_list_of_nodes


# Deprecated. Not working well enough
"""
def cluster_on_network_distance(graph, matrix_path,
                                avg_frag_size = 400,
                                std_frag_size = 120):
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from sklearn.manifold import MDS
    import matplotlib.pyplot as plt

    # Step 1. Prepare a matrix of distance between nodes
    # Initialize a matrix
    N = len(graph.nodes())
    distance_matrix = np.full((N, N), 999.0)
    plain_dist_matrix = np.full((N, N), 999)

    # Map nodes to indices
    index_to_node_tup = { i:tup for i, tup in enumerate(graph.nodes(data=True)) }

    # Simple lengths and weighted lengths are both precomputed
    paths = dict(nx.all_pairs_shortest_path(graph))
    total_paths_no = sum([len(d) for k, d in paths.items()])
    logger.info("There are {} nodes and {} edges in the graph. Theoretically, there should be {} paths between {} nodes, actually all_pairs_shortest_path found {} paths in the found shortest paths".format(graph.number_of_nodes(),
                                                                                                                                                                                                            graph.number_of_edges(),
                                                                                                                                                                                                            N**2, N, total_paths_no))

    # Fill in direct connections
    for i in range(0, N):
        node1 = index_to_node_tup[i][0]
        for j in range(i, N):
            node2 = index_to_node_tup[j][0]
            shortest_path = paths.get(node1, {}).get(node2, None)
            if i == j:
                distance_matrix[i, j] = 0
            else:
                if not shortest_path:
                    plain_dist_matrix[i, j] = 999.0
                    plain_dist_matrix[j, i] = 999.0
                    distance_matrix[i, j] = 999.0
                    distance_matrix[j, i] = 999.0
                elif len(shortest_path) <= 2:
                    plain_dist_matrix[i, j] = 999.0
                    plain_dist_matrix[j, i] = 999.0
                    distance_matrix[i, j] = 999.0
                    distance_matrix[j, i] = 999.0
                else:
                    weighted_length = sum([graph[shortest_path[i]][shortest_path[i+1]]["distance"] for i in range(0, len(shortest_path)-1)])
                    distance_matrix[i, j] = 1/weighted_length
                    distance_matrix[j, i] = 1/weighted_length
                    plain_dist_matrix[i, j] = 1/len(shortest_path)
                    plain_dist_matrix[j, i] = 1/len(shortest_path)


    # Show me the matrix:
    logger.info("The distance matrix looks like:\n{}\n".format(distance_matrix))
    # Save the matrix to a file
    np.savetxt(matrix_path, distance_matrix, delimiter="\t")
    np.savetxt(matrix_path.replace(".txt", ".plain.txt"), plain_dist_matrix, delimiter="\t")
    logger.info("The distance matrix is saved to {}".format(matrix_path))
    logger.info("The plain distance matrix is saved to {}".format(matrix_path.replace(".txt", ".plain.txt"))

    # Step 2. Cluster the distance matrix
    bandwidth = 0.5
    logger.info(f"The bandwidth is {bandwidth}")

    mean_shift = MeanShift(bandwidth = bandwidth)
    mean_shift.fit(distance_matrix)
    cluster_labels = mean_shift.labels_

    logger.info("The MeanShift clustering found {} clusters".format(len(set(cluster_labels))))

    mds = MDS(n_components=2, dissimilarity="precomputed")
    X_2D = mds.fit_transform(distance_matrix)
    colors = [plt.cm.jet(x) for x in np.linspace(0, 1, len(np.unique(cluster_labels)))]
    # Plotting
    plt.figure(figsize=(16, 16))
    for label in np.unique(cluster_labels):
        plt.scatter(X_2D[cluster_labels == label, 0], X_2D[cluster_labels == label, 1], c=[colors[i]], label=f'Cluster {label}')

    plt.legend()
    plt.title('Clustered Points in 2D')
    plt.xlabel('MDS1')
    plt.ylabel('MDS2')

    # Save the plot to a file
    plt.savefig(matrix_path + "clustered.png")
    logger.info("The clustered plot is saved to {}".format(matrix_path + "clustered.png"))

    # Step 3. Map the cluster labels back to the nodes
    clusters = {}
    for i, label in enumerate(cluster_labels):
        if label == -1:
            clusters[str(uuid.uuid4())] = [index_to_node_tup[i]]
            continue

        if str(label) not in clusters:
            clusters[str(label)] = []
        node_tup = index_to_node_tup[i]
        clusters[str(label)].append(node_tup)

    logger.info("The MeanShift divide the FC intervals into {} clusters".format(len(clusters)))

    # Now we need to further cluster the clusters based on the similarity of the ploidy change of the nodes
    subdivided_clusters = {}
    for label, clusters in clusters.items():
        subdivided_matrix_path = matrix_path.replace(".txt", ".subdivided.{}.txt".format(label))
        #  Here clusters is just a list of node tuples: [(node, data_dict), (node, data_dict), ...]
        if len(clusters) > 1:
            subdivided_cluster = bin_on_ploidy_fold_change(clusters, label, subdivided_matrix_path)
            subdivided_clusters.update(subdivided_cluster)
        else:
            subdivided_clusters[label] = clusters

    logger.info("WE further bin the FC intervals into {} clusters based on the ploidy change similarity".format(len(subdivided_clusters)))
    return subdivided_clusters
"""


def query_connected_nodes(fc_nfc_list,
                          multiplex_graph = nx.Graph(),
                          threads = 10,
                          avg_frag_size = 500,
                          std_frag_size = 150,
                          graph_path = ""):
    # The input fc_nfc list is a list of tuples
    # Each tuple contains two items: (fc_node, nfc_nodes)
    # Where fc_node is tuple (chr, start, end)
    # Where nfc_nodes is [nfc_node1, nfc_node2, ... , nfc_nodeN]
    # Where nfc_node1 is (chr, start, end)
    '''
    Note that an FC node can also be an NFC node.
    '''
    fc_nfc_dict = { n:ns for n, ns in fc_nfc_list }
    # We also need to build a dict that records which FC node each NFC node maps to
    fc_nodes = [ n for n, ns in fc_nfc_list]
    nfc_nodes = [ nfc for n, ns in fc_nfc_list for nfc in ns ]
    # total_fc_bed = BedTool("\n".join(["\t".join([str(x) for x in fc[:3]] + [".", ".", fc[3]]) for fc in fc_nodes ]), from_string=True) # Cannot sort and merge

    # fc_internal_overlap = total_fc_bed.intersect(total_fc_bed, wao = True, s=True).to_dataframe(disable_auto_names = True,
    #                                                                                     names = ["chr_fc1", "start_fc1", "end_fc1", "name1", "score1", "strand1",
    #                                                                                              "chr_fc2", "start_fc2", "end_fc2", "name2", "score2", "strand2",
    #                                                                                              "overlap_len"] ).drop_duplicates()
    # enwrapped_fcs = []
    # for idx, row in fc_internal_overlap.iterrows():
    #     interval_1_size = int(row["end_fc1"]) - int(row["start_fc1"])
    #     interval_2_size = int(row["end_fc2"]) - int(row["start_fc2"])
    #     interval_1 = (row["chr_fc1"], int(row["start_fc1"]), int(row["end_fc1"]))
    #     interval_2 = (row["chr_fc2"], int(row["start_fc2"]), int(row["end_fc2"]))
    #     if interval_1 == interval_2:
    #         continue
    #     if row["overlap_len"] == interval_1_size:
    #         enwrapped_fcs.append(interval_1)
    #     elif row["overlap_len"] == interval_2_size:
    #         enwrapped_fcs.append(interval_2)
    # enwrapped_fcs = set(enwrapped_fcs)

    # logger.info(f"These FC overlapping SD regions are compeletely enwrapped by other SDs so leave them out: \n{enwrapped_fcs}\n\n")

    # # Filter out the enwrapped_fcs
    # fc_nfc_list = [ (n, ns) for n, ns in fc_nfc_list if n not in enwrapped_fcs ]
    # fc_nfc_dict = { n:ns for n, ns in fc_nfc_list }
    # fc_nodes = [ n for n, ns in fc_nfc_list ]

    fc_node_ploidy_change = {}
    for fc_node in fc_nodes:
        nfc_nodes = fc_nfc_dict[fc_node]
        fc_bed = convert_node_into_bed_obj(fc_node)
        nfc_bed = BedTool("\n".join(["\t".join([str(x) for x in nfc]) for nfc in nfc_nodes ]), from_string=True)
        fc_bed_size = fc_bed.sort().total_coverage()
        nfc_bed_size = nfc_bed.sort().total_coverage()
        ploidy_fold_change = (nfc_bed_size + fc_bed_size)/fc_bed_size
        fc_node_ploidy_change[fc_node] = 2* ploidy_fold_change

    # Now we just need to extract all the components in the multiplex_graph
    # The multiplex_graph is a nx.Graph object
    undirected_graph = multiplex_graph.to_undirected()
    components = list(nx.connected_components(undirected_graph))
    # Test which FC nodes are in the same component
    fc_components = [[ fc_node for fc_node in fc_nodes if fc_node in component ] for component in components ]

    iso_fc_nodes = list(dict.fromkeys([ n for c in fc_components for n in c if len(c) == 1 ]))
    logger.info(f"How many FC nodes are disconnected from other FC nodes: {len(iso_fc_nodes)}")
    disconnected_result = { "PCs": iso_fc_nodes, "SD_counterparts": list(dict.fromkeys([nfc for fc in iso_fc_nodes for nfc in fc_nfc_dict[fc]])) }

    overlap_fc_nodes = list(dict.fromkeys([ n for n in fc_nodes if n not in iso_fc_nodes ]))
    logger.info(f"How many FC nodes overlaps with homologous sequence of other FC nodes: {len(overlap_fc_nodes)}")

    same_component_fc_nodes = [ c for c in fc_components if len(c) > 1 ]
    component_delimited_fcnodes = pick_from_each_group(same_component_fc_nodes)
    logger.info(f"How many groups of FC nodes are there with components as delimiter: {len(component_delimited_fcnodes)}")

    final_fcnode_groups = bin_on_ploidy_fold_change(component_delimited_fcnodes, fc_node_ploidy_change)
    logger.info("Further divide the group by ploidy change similarity, now there are {} groups of FC nodes".format(len(final_fcnode_groups)))

    connected_results = []
    for group in final_fcnode_groups:
        sd_counterparts = [nfc_node for fcnode in group for nfc_node in fc_nfc_dict[fcnode]]
        assert isinstance(sd_counterparts[0], HOMOSEQ_REGION)
        connected_result = { "PCs":[ n for n in group ], "SD_counterparts": sd_counterparts }
        if len(group) == 0:
            logger.warning(f"There is an empty group in the component_delimited_fcnodes, please check the input same_component_fc_nodes: \n{same_component_fc_nodes}\nAnd component_delimited_fc_nodes: \n{component_delimited_fcnodes}\n\n")
            continue
        if len(sd_counterparts) == 0:
            logger.warning(f"There is an empty set of SD counterparts in the connected_result, please check the input fc_nfc_dict: \n{fc_nfc_dict}\nAnd the fc nodes group: \n{group}\n\n")
            continue
        connected_results.append(connected_result)

    logger.info(f"There are {len(connected_results) + 1} FC-NFC pairs for the preparation of realignment bed file\n\n")
    return disconnected_result, connected_results, fc_nfc_dict



def construct_folder_struc(base_folder="/paedyl01/disk1/yangyxt/indexed_genome/SD_priority_component_pairs",
                           label="",
                           logger = logger):
    parent_folder_name = label+ "_related_homo_regions"
    parent_folder_full_path = os.path.join(base_folder, parent_folder_name)

    if not os.path.exists(parent_folder_full_path):
        os.mkdir(parent_folder_full_path)

    total_bed_name = label + "_related_homo_regions.bed"
    total_bed_path = os.path.join(parent_folder_full_path, total_bed_name)

    counterparts_bed_name = label + "_counterparts_regions.bed"
    counterparts_bed_path = os.path.join(parent_folder_full_path, counterparts_bed_name)

    PC_folder_name = label
    PC_folder_full_path = os.path.join(parent_folder_full_path, PC_folder_name)

    if not os.path.exists(PC_folder_full_path):
        os.mkdir(PC_folder_full_path)

    PC_bed_name = label + ".bed"
    PC_bed_path = os.path.join(PC_folder_full_path, PC_bed_name)

    return {"base_folder_path":parent_folder_full_path,
            "PC_bed": PC_bed_path,
            "All_region_bed": total_bed_path,
            "Counterparts_bed": counterparts_bed_path}



def perform_bedtools_sort_and_merge(bed_file,
                                    output_bed_file=None,
                                    logger = logger):
    # Load your BED file
    bed = BedTool(bed_file)

    # Sort the BED file
    sorted_bed = bed.sort()

    # Merge the BED file
    merged_bed = sorted_bed.merge(s=True, c='4,5,6', o='first,first,first')

    if not output_bed_file:
        output_bed_file = bed_file

    # You can save the results to a new file
    merged_bed.saveas(output_bed_file)



def update_plain_file_on_md5(old_file, new_file, logger=logger):
    import hashlib

    if os.path.exists(old_file):
        old_md5 = hashlib.md5(open(old_file,'r').read().encode()).hexdigest()
    else:
        old_md5 = str(uuid.uuid4())

    if not os.path.exists(new_file):
        raise FileNotFoundError("The new file {} does not exist so no updates should be carried out.".format(new_file))

    if os.stat(new_file).st_size == 0:
        raise FileExistsError("The new file {} input is completely empty. Quit using it to update the original file {}".format(new_file, old_file))

    new_md5 = hashlib.md5(open(new_file,'r').read().encode()).hexdigest()

    if new_md5 == old_md5:
        logger.warning("The new file {} shares the identical content with the old one {} so no updates should be carried out. And the new file {} should be deleted".format(new_file,
                                                                                                                                                                                old_file,
                                                                                                                                                                                new_file))
        os.remove(new_file)
        return False
    else:
        import shutil
        logger.info("The new file {} is different from the old one {} so the old one will be replaced with the new one.".format(new_file, old_file))
        shutil.move(new_file, old_file)
        executeCmd(f"ls -lht {old_file}")
        return True


@error_handling_decorator
def establish_beds_per_PC_cluster(cluster_dict={"PCs":{},
                                                "SD_counterparts":{}},
                                  base_folder = "/paedyl01/disk1/yangyxt/indexed_genome/SD_priority_component_pairs",
                                  label = "PC0",
                                  ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                  length = None,
                                  avg_frag_size = 400,
                                  std_frag_size = 140,
                                  threads = 2,
                                  logger = logger):

    """
    input cluster dict now looks like this:
    {
        "PCs": {0: [], 1: [], 2: []},
        "SD_counterparts": {0: [], 1: [], 2: []}
    }
    """

    logger.info("The input cluster_dict is:\n{}\n".format(cluster_dict))
    assert isinstance(cluster_dict["SD_counterparts"][0][0], HOMOSEQ_REGION)
    paths = construct_folder_struc(base_folder=base_folder, label=label, logger=logger)
    logger.info("The PC bed file is {}, the counterparts bed file is {}, the total bed file is {}".format(paths["PC_bed"], paths["Counterparts_bed"], paths["All_region_bed"]))

    # First convert the disconnected nodes to beds, each node is a tuple consisting of three values (chr, start, end)
    tmp_id = str(uuid.uuid4())
    tmp_pc_bed = paths["PC_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_pc_bed = paths["PC_bed"].replace(".bed", ".raw.bed")
    tmp_counterparts_bed = paths["Counterparts_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_counterparts_bed = paths["Counterparts_bed"].replace(".bed", ".raw.bed")
    tmp_total_bed = paths["All_region_bed"].replace(".bed", "." + tmp_id + ".bed")
    raw_total_bed = paths["All_region_bed"].replace(".bed", ".raw.bed")

    with open(tmp_pc_bed, "w") as f:
        for idx, records in cluster_dict["PCs"].items():
            for record in records:
                if len(record) >= 3:
                    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3]]) + "\n")
                elif len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3]]) + "\n")

    subprocess.run(f"cp -f {tmp_pc_bed} {raw_pc_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_pc_bed, logger=logger)
    update_plain_file_on_md5(paths["PC_bed"], tmp_pc_bed, logger=logger)

    # Then compose the counterparts region bed file
    with open(tmp_counterparts_bed, "w") as f:
        for idx, records in cluster_dict["SD_counterparts"].items():
            for record in records:
                assert isinstance(record, HOMOSEQ_REGION), "The record in the SD_counterparts list is not a HOMOSEQ_REGION object: {}".format(record)
                if len(record) == 2:
                    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3]]) + "\n")
                elif len(record) >= 3:
                    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3]]) + "\n")
    subprocess.run(f"cp -f {tmp_counterparts_bed} {raw_counterparts_bed}", shell=True)

    # Now we need to subtract the PC from counterpart beds incase they have overlaps,  force strandness
    BedTool(tmp_counterparts_bed).subtract(BedTool(paths["PC_bed"]), s=True).saveas(tmp_counterparts_bed)
    perform_bedtools_sort_and_merge(tmp_counterparts_bed, logger=logger)

    update_plain_file_on_md5(paths["Counterparts_bed"], tmp_counterparts_bed, logger=logger)

    # Then concatenate the two bed file together to generate the total region bed file
    # Here you want to trigger the __iter__ method in HOMOSEQ_REGION instead of the __getitem__ method
    with open(tmp_total_bed, "w") as f:
        for idx, records in cluster_dict["PCs"].items():
            for record in records:
                if len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3], f"FC:{label}_{idx}"]) + "\n")
                elif len(record) >= 3:
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3], f"FC:{label}_{idx}"]) + "\n")
        for idx, records in cluster_dict["SD_counterparts"].items():
            for record in records:
                if len(record) == 2:
                    f.write("\t".join([str(value) for value in record[0]][:3] + [".", ".", record[0][3], f"NFC:{label}_{idx}"]) + "\n")
                elif len(record) >= 3:
                    f.write("\t".join([str(value) for value in record][:3] + [".", ".", record[3], f"NFC:{label}_{idx}"]) + "\n")

    subprocess.run(f"cp -f {tmp_total_bed} {raw_total_bed}", shell=True)
    perform_bedtools_sort_and_merge(tmp_total_bed, logger = logger)
    update_plain_file_on_md5(paths["All_region_bed"], tmp_total_bed, logger=logger)

    contig_sizes = ".".join(ref_genome.split(".")[:-1]) + ".contigsize.genome"
    # Alongside preparing the bed files, we can also perform masked genome preparation
    masked_genome = Genome(ref_genome).mask(paths["PC_bed"], avg_frag_size = avg_frag_size, std_frag_size=std_frag_size, genome=contig_sizes)
    executeCmd("ls -lht {}".format(paths["PC_bed"]), logger=logger)
    assert os.path.getmtime(masked_genome) > os.path.getmtime(paths["PC_bed"])

    # And we need to prepare the intrinsic VCF based on the masked genome and generated BED files
    bam_path = getIntrinsicVcf(pc_bed = paths["PC_bed"],
                    all_homo_regions_bed = paths["All_region_bed"],
                    counter_bed = paths["Counterparts_bed"],
                    pc_masked = masked_genome,
                    avg_frag_size = avg_frag_size,
                    std_frag_size = std_frag_size,
                    threads = 2,
                    logger = logger)
    return bam_path


def complement_bed(bedf: str, rg: Genome, padding_size = 400) -> str:
    """
    This function takes the path of a BED file and makes complement of it.

    Return value is the path of complement BED.
    """

    import os

    tmp_tag = str(uuid.uuid4())

    _merged_out = bedf[:-3] + "merged.bed.tmp"
    merged_out = bedf[:-3] + "merged.bed"
    tmp_output = merged_out + "." + tmp_tag

    if "merged" in bedf:
        logger.debug(f"Complementary BED files found - {bedf}. Skipping ... ")
        return bedf

    contig_genome = ".".join(rg.path.split(".")[:-1]) + ".contigsize.genome"
    if not os.path.exists(contig_genome):
        rg.getContigGenome()

    executeCmd(f"bedtools slop -b {padding_size} -g {contig_genome} -i {bedf} | bedtools merge -d 150 -i - > {_merged_out} ")

    cmd = f"bedtools subtract -a {rg.total_bed} -b {_merged_out} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {tmp_output} "
    executeCmd(cmd)

    update_plain_file_on_md5(merged_out, tmp_output)
    os.remove(_merged_out)

    return merged_out


def imap_establish(tup_args):
    return establish_beds_per_PC_cluster(*tup_args)


def imap_traverse(tup_args):
    return traverse_network_to_get_homology_counterparts(*tup_args)


def convert_nodes_into_hierachical_beds(disconnected_result = {'PCs':[], 'SD_counterparts':[]},
                                        connected_results = [],
                                        fc_nfc_dict = {},
                                        output_folder="/paedyl01/disk1/yangyxt/indexed_genome/SD_priority_component_pairs",
                                        ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                        nthreads = 12,
                                        length = None,
                                        avg_frag_size = 400,
                                        std_frag_size = 140):
    # Establish a series of labels for PC-counterparts pairs
    # For disconnected query nodes cluster, we simply name them collectively as PC0
    # For connected query nodes cluster, we name them as PC1, PC2, PC3, ...
    import shutil
    import multiprocessing as mp
    from itertools import repeat
    import resource

    # Second convert the connected nodes to beds
    all_results = [disconnected_result] + connected_results
    logger.info("Start to put beds into the rest of PC-counterparts paired beds into subfolders")
    # We kind of need to restructure the all_results to pass the information of fc-nfc pairs to the establish_beds_per_PC_cluster function
    new_results = []
    for result in all_results:
        new_result = {"PCs": {}, "SD_counterparts": {}}
        for i in range(0, len(result["PCs"])):
            fc_node = result["PCs"][i]
            new_result['PCs'][i] = [fc_node]
            new_result['SD_counterparts'][i] = fc_nfc_dict[fc_node]
        new_results.append(new_result)

    labels = [ "PC" + str(n) for n in range(0, len(all_results))]

    pool = Pool(int(round(nthreads/2)))
    results = pool.imap_unordered(imap_establish, zip(new_results,
                                                      repeat(output_folder),
                                                      labels,
                                                      repeat(ref_genome),
                                                      repeat(length),
                                                      repeat(avg_frag_size),
                                                      repeat(std_frag_size)))
    i = 0
    intrinsic_bams = []
    for success, result, logs in results:
        print(f"\n\n*********************************** {i}_subprocess_start ***************************************", file = sys.stderr)
        if not success:
            error_mes, tb_str = result
            logger.error(f"An error occurred: {error_mes}\nTraceback: {tb_str}\n")
        else:
            intrinsic_bams.append(result)
        print(logs, file = sys.stderr)
        print(f"*********************************** {i}_subprocess_end ***************************************\n\n", file = sys.stderr)
        i+=1

    pool.close()
    if not all([t[0] for t in results]):
        raise RuntimeError("Error happened during the parallel execution of establish_beds_per_PC_cluster. So stop the total python script here.")

    execute = False
    valid_intrinsic_bams = " ".join(intrinsic_bams)
    total_intrinsic_bam = os.path.join(output_folder, "total_intrinsic_alignments.bam")
    test_cmd = f"bash {bash_utils_hub} check_bam_validity {total_intrinsic_bam}"
    try:
        executeCmd(test_cmd)
    except RuntimeError:
        logger.info(f"The merged intrinsic alignment BAM file {total_intrinsic_bam} is not ready")
        execute = True
    else:
        if all([os.path.getmtime(total_intrinsic_bam) > os.path.getmtime(vb) for vb in intrinsic_bams]):
            execute = False
        else:
            execute = True

    intrinsic_bam_header = total_intrinsic_bam.replace(".bam", ".bam.header")
    cmd = f"bash {bash_utils_hub} modify_bam_sq_lines {intrinsic_bams[0]} {ref_genome} {intrinsic_bam_header}"
    executeCmd(cmd, logger=logger)

    intrinsic_bam_list = total_intrinsic_bam.replace(".bam", ".bams.list.txt")
    with open(intrinsic_bam_list, "w") as f:
        f.write("\n".join(intrinsic_bams))

    cmd = f"samtools merge -@ {nthreads} -h {intrinsic_bam_header} -b {intrinsic_bam_list} -o - | \
            samtools sort -O bam -o {total_intrinsic_bam} && \
            samtools index {total_intrinsic_bam} && \
            ls -lht {total_intrinsic_bam} || \
            echo Failed to concatenate all the filtered realigned BAM files. It wont be a fatal error but brings troubles to debugging and variant tracing."
    if execute:
        executeCmd(cmd)

    # Third remove extra PC folders that left from previous generation
    subdir_gen = os.walk(output_folder)
    # Prepare all target folder names in this time's generation, the naming syntax should follow the function called construct_folder
    if len(connected_results) > 0:
        target_folder_names = set([f"PC{x}_related_homo_regions" for x in range(0, len(connected_results)+1)])
    else:
        target_folder_names = ["PC0_related_homo_regions"]
    logger.info(f"This time, we only have these PC regions generated: {target_folder_names}")
    first_level_dirs = next(subdir_gen)[1]  # The first item should be a 3-item tuple containin: 1. Current dir full path 2. All subdir names 3. All subfile names
    for subdir in first_level_dirs:
        if subdir not in target_folder_names:
            # Remove the subdir forcely
            shutil.rmtree(os.path.join(output_folder, subdir))

    # Fourth, extract all PC*_related_homo_regions.bed file and concat them together.
    first_level_dirs = next(os.walk(output_folder))[1]
    # Build a temp file to store the names of all sub files, if directly put all the names following cat, it might cause argument list too long error
    tmp_lst = prepare_tmp_file(suffix=".txt").name
    with open(tmp_lst, "w") as tf:
        total_lines = []
        for bed_file in [os.path.join(output_folder, subdir, subdir + ".bed") for subdir in first_level_dirs]:
            with open(bed_file, "r") as bf:
                total_lines = total_lines + bf.readlines()
        tf.write("\n".join([l.strip("\n") for l in total_lines if len(l) > 1]))

    executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total}".format(all = tmp_lst,
                                                                                                        total = os.path.join(output_folder, "all_PC_related_homo_regions.bed")))

    executeCmd("cat {all} | bedtools sort -i stdin | bedtools merge -i stdin > {total} && ls -lht {total} && rm {all}".format(all = tmp_lst,
                                                                                                        total = os.path.join(output_folder, "all_PC_regions.bed")))


    # Fifth, extract all intrinsic vcfs and use bcftools to concat them together
    intrinsic_vcfs = []
    for root, dirs, files in os.walk(output_folder):
        for file in files:
            if file.endswith(f'{int(avg_frag_size)}.vcf.gz'):
                intrinsic_vcfs.append(os.path.join(root, file))

    vcf_list_file = os.path.join(output_folder, "intrinsic_vcf.lst")
    final_intrinsic_vcf = os.path.join(output_folder, "all_pc_region_intrinsic_variants.vcf.gz")
    tmp_intrinsic_vcf = os.path.join(output_folder, "all_pc_region_intrinsic_variants.tmp.vcf.gz")
    with open(vcf_list_file, "w") as f:
        for v in intrinsic_vcfs: f.write(v + "\n")

    executeCmd(f"bcftools concat -o {tmp_intrinsic_vcf} -a --no-version -Oz -f {vcf_list_file} && bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp -o {final_intrinsic_vcf} -Oz {tmp_intrinsic_vcf} && ls -lh {final_intrinsic_vcf} && rm {tmp_intrinsic_vcf}")
    return



def graph_to_sorted_bedtool(G):
    # Extract nodes from the graph and convert them into bed format
    bed_format = ["\t".join(map(str, node)) for node in G.nodes]

    # Create a BedTool object
    bedtool = BedTool("\n".join(bed_format), from_string=True).sort()

    return bedtool



def is_enclosed(intervals, new_interval):
    from bisect import bisect_right
    # Sort intervals by the first element
    intervals.sort()

    # Find insertion point for new interval
    index = bisect_right(intervals, new_interval)

    # Check the interval before the insertion point
    if index > 0 and intervals[index - 1][0] <= new_interval[0] and intervals[index - 1][1] >= new_interval[1]:
        return True

    # Check the interval after the insertion point
    if index < len(intervals) and intervals[index][0] <= new_interval[0] and intervals[index][1] >= new_interval[1]:
        return True

    return False


def calculate_interval_overlaps(interval1, interval2, fraction_select = None):
    overlap_span = min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])
    if overlap_span < 0:
        return 0
    else:
        interval1_size = interval1[1] - interval1[0]
        interval2_size = interval2[1] - interval2[0]
        if fraction_select:
            return [overlap_span/interval1_size, overlap_span/interval2_size][fraction_select]
        else:
            return max(overlap_span/interval1_size, overlap_span/interval2_size)


def calculate_internal_overlaps(intervals):
    raw_bed = BedTool("\n".join(["\t".join([str(x) for x in interval]) for interval in intervals]))
    raw_size = raw_bed.to_dataframe(disable_auto_names=True, names=["chr", "start", "end"]).apply(lambda x: (float(x["end"]) - float(x["start"])), axis=1).sum()
    merged_bed = raw_bed.sort().merge()
    merged_size = merged_bed.to_dataframe(disable_auto_names=True, names=["chr", "start", "end"]).apply(lambda x: (float(x["end"]) - float(x["start"])), axis=1).sum()
    internal_overlap = raw_size - merged_size
    return internal_overlap, merged_size


def calculate_overlap_to_central_interval(central_interval, peripheral_intervals):
    central_bed = BedTool("\t".join([str(x) for x in central_interval]), from_string=True).sort().merge()
    peripheral_bed = BedTool("\n".join(["\t".join([str(x) for x in interval]) for interval in peripheral_intervals]), from_string=True).sort().merge()
    intersect_bed = central_bed.intersect(peripheral_bed, wo=True)
    total_overlap_size = intersect_bed.to_dataframe(disable_auto_names=True, names=["chrom", "start", "end",
                                                                                    "pchr", "pstart", "pend",
                                                                                    "overlap_len"]).loc[:, "overlap_len"].astype(float).sum()
    return total_overlap_size


def solve_lp_for_coverage(central_interval, peripheral_intervals):
    # First try to find any peripheral intervals that nearly completely enclose the central interval
    peripheral_interval_dict = {"peripheral_interval": peripheral_intervals,
                                "central_interval_overlap": [calculate_overlap_to_central_interval(central_interval, [peripheral_interval]) for peripheral_interval in peripheral_intervals],
                                "peripheral_interval_size": [peripheral_interval[2] - peripheral_interval[1] for peripheral_interval in peripheral_intervals]}
    peri_df = pd.DataFrame(peripheral_interval_dict).sort_values(by=["central_interval_overlap",
                                                                     "peripheral_interval_size"], ascending=[False, True]).reset_index(drop=True)
    if peri_df.loc[:, "central_interval_overlap"].tolist()[0] >= 0.95:
        return [peri_df.iat[0, 0]]

    # Now start to find the best combination of peripheral intervals
    total_combs = []
    for i in range(1, len(peripheral_intervals) + 1):
        list_of_comb = itertools.combinations(peripheral_intervals, i)
        total_combs = total_combs + list_of_comb

    tups = [calculate_internal_overlaps(comb) for comb in total_combs]
    comb_dict = {"peripheral_interval": list(range(0, len(total_combs))),
                 "central_interval_overlap": [calculate_overlap_to_central_interval(central_interval, comb) for comb in total_combs],
                 "peripheral_interval_size": [t[1] for t in tups],
                 "peripheral_internal_overlap": [t[0] for t in tups]}

    comb_df = pd.DataFrame(comb_dict).sort_values(by=[ "central_interval_overlap",
                                                       "peripheral_interval_size",
                                                       "peripheral_internal_overlap" ], ascending=[False, True, True]).reset_index(drop=True)

    comb_df["final_sort_index"] = comb_df["central_interval_overlap"].astype(float) * 100 - comb_df["peripheral_interval_size"].astype(float) * 10 - comb_df["peripheral_internal_overlap"].astype(float) * 1
    comb_df = comb_df.sort_values(by=["final_sort_index"], ascending=False).reset_index(drop=True)
    select_comb_ind = comb_df.iat[0, 0]
    return total_combs[select_comb_ind]



def imap_pickup_regions_per_exon(tup_args):
    return pickup_regions_per_exon(*tup_args)


def pickup_regions_per_exon(groupdf, group_interval_chr = "chr_query",
                                     group_interval_start = "start_query",
                                     group_interval_end = "end_query",
                                     avg_frag_size = 560,
                                     std_frag_size = 140):
    '''
    This is a function that:
    1. select intervals (recorded in the first three columns, recording the intervals of SDs) that covers the most region of the fixed interval
    2. The later three columns are intervals of target region (by d efault exons and splicing sites), they are also recording the fixed interval
    3. keep the selected intervals have the minimum combined span
    4. This function assumes the input groupdf has all recorded intervals in the same chromosome
    '''
    from operator import itemgetter
    # Fixed interval in this context are usually exons
    fixed_chr, fixed_start, fixed_end = groupdf.loc[:, group_interval_chr].tolist()[0], int(groupdf.loc[:, group_interval_start].to_list()[0]), int(groupdf.loc[:, group_interval_end].to_list()[0])  # Ensure these are integers

    # Create a list of tuples, where each tuple represents a interval (start, end)
    intervals = sorted([(chrom, int(start), int(end), strand) for chrom, start, end, strand in groupdf.iloc[:, :4].itertuples(index=False)], key=itemgetter(2))

    intervals_coverage = [(it, calculate_interval_overlaps((int(fixed_start), int(fixed_end)), it[1:3], 0)) for it in intervals]
    total_covered_intervals = [it for it, cov in intervals_coverage if cov > 0.95]

    if len(total_covered_intervals) > 0:
        selected_intervals = [sorted(total_covered_intervals, key=lambda t:abs(t[2] - t[1] - (avg_frag_size + 1.96 * std_frag_size)))[0]]
    else:
        logger.warning(f"Cant find any overlapping intervals that largely covers the fixed interval {fixed_chr}:{fixed_start}-{fixed_end}, so we will try to find the best combination of intervals")
        selected_intervals = solve_lp_for_coverage((fixed_chr, fixed_start, fixed_end), intervals)

    # Now we slice the groupdf to only keep the selected intervals
    selected_rows = groupdf.loc[[interval in selected_intervals for interval in groupdf.iloc[:, :4].itertuples(index=False)], :]
    return selected_rows

def string_to_tuple(string):
    try:
        return ast.literal_eval(string)
    except ValueError:
        raise ValueError(f"Failed to convert string to tuple: {string}")


def read_graphml(graph_path):
    literal_graph = nx.read_graphml(graph_path)
    # We need to convert the string back to tuple
    graph = literal_graph.__class__()
    for node, data in literal_graph.nodes(data=True):
        graph.add_node(string_to_tuple(node), **data)

    for node1, node2, data in literal_graph.edges(data=True):
        graph.add_edge(string_to_tuple(node1), string_to_tuple(node2), **data)

    return graph


def read_pickle_graph(graph_path):
    with open(graph_path, "rb") as f:
        graph = pickle.load(f)
    return graph


def write_pickle_graph(graph_path, graph):
    with open(graph_path, "wb") as f:
        pickle.dump(graph, f)

def extract_FC_NFC_pairs_from_graph(query_nodes, directed_graph, graph_path = "", avg_frag_size = 500, std_frag_size = 150, threads=12):
    # First we need to annotate the graph about which nodes are query_nodes
    # query nodes are dataframe with chr, start, end ,three columns
    query_nodes = list(zip(query_nodes["chr"], query_nodes["start"], query_nodes["end"], query_nodes["strand"]))
    logger.info(f"How many query nodes we have ? {len(query_nodes)}, stored as a list of {type(query_nodes[0])}")
    logger.info(f"We got a graph with {len(directed_graph.nodes())} nodes and {len(directed_graph.edges())} edges")

    for node in directed_graph.nodes():
        if node in query_nodes:
            directed_graph.nodes[node]["query_node"] = True
        else:
            directed_graph.nodes[node]["query_node"] = False
    logger.info("Now we have marked the query nodes in the graph")

    unfiltered_graph = directed_graph.copy()
    # Then we need to filter out nodes that are too small
    cutoff = max(140, avg_frag_size - .675 * std_frag_size)
    tobe_removed_nodes = []
    for node in directed_graph.nodes(data=True):
        size = float(node[0][2]) - float(node[0][1])
        directed_graph.nodes[node[0]]["size"] = int(size)
        if size <= cutoff:
            # Dont modify collections during the iteration
            tobe_removed_nodes.append(node[0])
    tobe_removed_nodes = list(dict.fromkeys(tobe_removed_nodes))
    directed_graph.remove_nodes_from(tobe_removed_nodes)
    logger.info("Now we have removed nodes that are too small, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))

    # Then we filter out the edges that the overlap size is too small
    tobe_removed_edges = []
    for edge in directed_graph.edges(data=True):
        overlap = edge[2].get("overlap", False)
        if overlap:
            smaller_node = edge[1]
            smaller_node_size = smaller_node[2] - smaller_node[1]
            weight = float(edge[2]["weight"])
            overlap_size = smaller_node_size * weight
            if overlap_size < cutoff and weight < 0.6:
                tobe_removed_edges.append(edge[:2])
        if edge[0] == edge[1]:
            # Remove self loop edges
            tobe_removed_edges.append(edge[:2])
    tobe_removed_edges = list(dict.fromkeys(tobe_removed_edges))
    directed_graph.remove_edges_from(tobe_removed_edges)  # This directed_graph also needs to be returned for further saving as a file
    logger.info("Now we have removed edges that implies an overlap with a size too small, the graph now has {} nodes and {} edges".format(len(directed_graph.nodes()), len(directed_graph.edges())))

    """
    Then we can handle the BFS search to another function and perform this on each query node in parallel
    """
    assert len(query_nodes) > 0, "No query nodes are found in the graph"
    undirected_graph = directed_graph.to_undirected()
    with Pool(threads) as pool:
        hom_results = pool.imap_unordered(imap_traverse,
                                          zip(query_nodes,
                                              repeat(undirected_graph),
                                              repeat(avg_frag_size),
                                              repeat(std_frag_size)))
        # Now we need to identify the FC-NFC pairs and annotate them on the directed graph
        i = 0
        filtered_results = set([])
        for tup in hom_results:
            if len(tup) == 0:
                continue
            i += 1
            filtered_results.add(tup)
            qnode, counterparts_nodes = tup
            # Here the tuple should contain (tuple_node, [list of nodes in self-defined class HOMOSEQ_REG])
            assert type(qnode) == tuple
            try:
                assert isinstance(counterparts_nodes[0], HOMOSEQ_REGION)
            except IndexError:
                # There is one interval (chrX, 70902050, 71018311) in the WGAC database, that it only corresponds with itself, so it will be filtered out
                logger.warning(f"No counterparts nodes are found for query node {qnode}, so we will skip this query node")
                continue

            # First we need to annotate the query node
            unfiltered_graph.nodes[qnode]["FRA_dest"] = str(i)
            # Then we need to annotate the counterparts nodes
            for cnode in counterparts_nodes:
                unfiltered_graph.nodes[cnode]["FRA_source"] = (unfiltered_graph.nodes[cnode].get("FRA_source", "") + "," + str(i)).lstrip(",")

    # After the annotation with the FC-NFC pairs in the graph, we need to save it to a graphml file
    if graph_path:
        nx.write_graphml(unfiltered_graph, graph_path)

    return filtered_results


class HOMOSEQ_REGION:
    def __init__(self, chrom, start, end, strand, ref_fasta=""):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.ref_fasta = ref_fasta
        self.rela_start = 0
        self.rela_end = end - start
        self.data = tuple([self.chrom, self.start, self.end, self.strand])

    def __hash__(self):
        return hash(self.data)

    def __repr__(self):
        return ", ".join(str(x) for x in [self.chrom, self.start, self.end, self.strand])

    def __eq__(self, other):
        return self.chrom == other[0] and self.start == other[1] and self.end == other[2] and self.strand == other[3]

    def __iter__(self):
        for item in [self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand]:
            yield item

    def __getitem__(self, index):
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self.data))
            return [self.data[i] for i in range(start, stop, step)]
        else:
            return self.data[index]

    def __len__(self):
        return len(self.data)

    def fix_coord(self):
        return tuple([self.chrom, self.start + self.rela_start, self.start + self.rela_end, self.strand])


def evaluate_neighbor_nodes_per_node(qnode, graph, ori_qnode,
                                     visited_nodes = [],
                                     extend_type = None,
                                     avg_frag_size = 500,
                                     std_frag_size = 150):
    # Here the input qnode should be as a type of self-defined class obj instead of a pure tuple
    # Here the graph is already an undirected graph so bfs_edges can neglect the directions of edges.
    # visited_nodes = set(visited_nodes)
    extend_nodes = []
    counter_nodes = []

    neighbor_edges = nx.bfs_edges(graph, source = qnode, depth_limit=1)
    for edge in neighbor_edges:
        # Here u,v are tuples
        u, v = edge
        cnode = v if u == qnode else u
        cnode_reg = HOMOSEQ_REGION(*cnode)
        if cnode_reg in visited_nodes:
            continue

        edge_attr = graph.get_edge_data(u, v, None)
        if not edge_attr:
            edge_attr = graph.get_edge_data(v, u, None)
        if not edge_attr:
            continue

        qnode_attr = graph.nodes[qnode]
        cnode_attr = graph.nodes[cnode]
        logger.info("The attributes of query node and counterpart node are: {} and {}".format(qnode_attr, cnode_attr))

        if edge_attr.get("type", None) == "segmental_duplication":
            if isinstance(qnode, HOMOSEQ_REGION):
                cnode_reg.rela_start = qnode.rela_start
                cnode_reg.rela_end = qnode.rela_end
            counter_nodes.append(cnode_reg)
            extend_nodes.append((cnode_reg, "SD"))
        elif edge_attr.get("overlap", None):
            # The cnode is overlapping with the qnode.
            if cnode_attr.get("size") <= qnode_attr.get("size"):
                # Qnode is larger
                if edge_attr.get("weight") * cnode_attr.get("size") >= avg_frag_size - 0.675 * std_frag_size:
                    if cnode_reg.start >= qnode[1] and cnode_reg.end <= qnode[2]:
                        # Qnode completely enwraps the cnode
                        cnode_rela_start = 0
                        cnode_rela_end = cnode_reg.end - cnode_reg.start
                    elif cnode_reg.start >= qnode[1]:
                        # Qnode is larger but cnode end is larger than the qnode end
                        cnode_rela_start = 0
                        cnode_rela_end = qnode[2] - cnode_reg.start
                    else:
                        # Qnode is larger but the cnode start is smaller than the qnode start
                        cnode_rela_start = qnode[1] - cnode_reg.start
                        cnode_rela_end = cnode_reg.end - cnode_reg.start
                    cnode_reg.rela_start = cnode_rela_start
                    cnode_reg.rela_end = cnode_rela_end
                    if extend_type == "overlap":
                        visited_nodes.append(cnode_reg)
                    else:
                        extend_nodes.append((cnode_reg, "overlap"))
                else:
                    visited_nodes.append(cnode_reg)
            else:
                # Counter node is larger and qnode and cnode are physically overlapping with each other
                # It does not mean that cnode compeletely enclose the qnode
                # First decide the relative position of cnode
                if edge_attr.get("weight") * qnode_attr.get("size") >= avg_frag_size - 0.675 * std_frag_size:
                    if cnode_reg.start <= qnode[1] and cnode_reg.end >= qnode[2]:
                        cnode_rela_start = qnode[1] - cnode_reg.start
                        cnode_rela_end = cnode_rela_start + qnode_attr.get("size")
                    elif cnode_reg.start <= qnode[1]:
                        cnode_rela_start = qnode[1] - cnode_reg.start
                        cnode_rela_end = cnode_attr.get("size")
                    else:
                        cnode_rela_start = 0
                        cnode_rela_end = qnode[2] - cnode_reg.start
                    cnode_reg.rela_start = cnode_rela_start
                    cnode_reg.rela_end = cnode_rela_end
                    if extend_type == "overlap":
                        visited_nodes.append(cnode_reg)
                    else:
                        extend_nodes.append((cnode_reg, "overlap"))
                else:
                    visited_nodes.append(cnode_reg)

    return counter_nodes, extend_nodes, visited_nodes



def traverse_network_to_get_homology_counterparts(qnode, directed_graph, avg_frag_size = 500, std_frag_size = 150):
    """
    The directed graph input here should be handled in advance to remove some nodes that are too small
    Or remove some edges that their weights are too low

    Then we need to extract the FC-NFC pairs by performing BFS for each query node
    For each pair of FC-NFC pairs, FC is single genomic interval and NFC can be a list of genomic intervals
    1. We need to perform BFS layer by layer, each BFS only performs at depth equals to 1
    2. We need to make sure the following BFS search does not traverse back to the nodes that have been visited
    3. We need to filter out nodes based on the connecting edge type and attributes value

    Note that the qnode is from the graph bed that has not been filtered by the node size
    Meanwhile, the directedg graph input here has already been filtered by the node size.
    """
    if type(directed_graph) == str:
        if re.search(r"\.graphml$", directed_graph):
            directed_graph = read_graphml(directed_graph)
        else:
            directed_graph = read_pickle_graph(directed_graph)

    qnode_size = BedTool("\t".join([str(x) for x in qnode]), from_string=True).sort().total_coverage()
    if qnode not in directed_graph.nodes():
        logger.warning("The query node {} is not in the directed graph".format(qnode))
        return tuple([])

    qnode_reg = HOMOSEQ_REGION(*qnode)
    qnode_size = qnode[2] - qnode[1]
    visited_nodes = [qnode_reg]
    logger.info("Start to search the homology counterparts for {} (first bfs with visited nodes set as qnode)".format(qnode))
    counterparts_nodes, extend_nodes, visited_nodes = evaluate_neighbor_nodes_per_node(qnode_reg,
                                                                                       directed_graph,
                                                                                       qnode,
                                                                                       visited_nodes = visited_nodes,
                                                                                       avg_frag_size = avg_frag_size,
                                                                                       std_frag_size = std_frag_size)

    visited_nodes = list(dict.fromkeys(visited_nodes + counterparts_nodes + extend_nodes))
    extend_nodes = {1: extend_nodes}
    new_counter = len(counterparts_nodes)
    n = 1
    while len([v for k,vs in extend_nodes.items() for v in vs]) > 0:
        '''
        The only conditions to break the while loop is:
        1. No new counterparts nodes are found
        &
        2. No new extend nodes are to be checked
        '''
        n+=1
        counter_len = len(counterparts_nodes)
        extend_nodes = { k:sorted(vs, key = lambda t: abs(t[0].end - t[0].start - qnode_size)) for k, vs in extend_nodes.items() }
        closest_layer = min([k for k in extend_nodes.keys() if len(extend_nodes[k]) > 0])
        closest_layer_enodes = extend_nodes[closest_layer]
        enode, enode_type = closest_layer_enodes.pop(0)
        iter_counter_nodes, iter_extend_nodes, visited_nodes = evaluate_neighbor_nodes_per_node(enode,
                                                                                                directed_graph,
                                                                                                qnode,
                                                                                                visited_nodes = visited_nodes,
                                                                                                extend_type = enode_type,
                                                                                                avg_frag_size = avg_frag_size,
                                                                                                std_frag_size = std_frag_size)
        counterparts_nodes = counterparts_nodes + iter_counter_nodes
        extend_nodes[closest_layer + 1] = extend_nodes.get(closest_layer + 1, []) + iter_extend_nodes
        visited_nodes = list(dict.fromkeys(visited_nodes + iter_counter_nodes + iter_extend_nodes))
        new_counter = len(counterparts_nodes) - counter_len
        assert all([type(cnode) == HOMOSEQ_REGION for cnode in counterparts_nodes]), "The counterparts nodes are not of type HOMOSEQ_REGION"
        counter_bed_size = BedTool("\n".join(["\t".join([str(x) for x in cnode]) for cnode in counterparts_nodes]), from_string=True).sort().merge().total_coverage()
        logger.info(f"Doing {n}th bfs for node {enode}, getting {new_counter} new counterparts nodes: {iter_counter_nodes} and number of extend nodes is {len([v for k,vs in extend_nodes.items() for v in vs])} by the end of current bfs")
        if len(counterparts_nodes) > 9 or counter_bed_size/qnode_size > 9 or closest_layer > 3:
            logger.info(f"Stop the bfs search for node {qnode} because the number of counterparts nodes has reached 10. We might already collect enough reads from other regions to the current target region")
            break

    # Now we have qnode and its counterparts nodes, we can return them
    if len(counterparts_nodes) == 0:
        logger.warning(f"This query node {qnode} cannot identify any region in the graph that shares homologous sequences\n\n")
    else:
        assert all([type(cnode) == HOMOSEQ_REGION for cnode in counterparts_nodes]), "The counterparts nodes are not of type HOMOSEQ_REGION"
        counter_bed_str = "\n".join(["\t".join([str(x) for x in cnode]) for cnode in counterparts_nodes])
        logger.info(f"This query node {qnode} has {len(counterparts_nodes)} counterparts nodes that share homologous sequences and they are:\n{counter_bed_str}\n\n")

    return (qnode, tuple(counterparts_nodes))

def deploy_PCs_for_SDrecall_main(ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                                 base_dir = "/paedyl01/disk1/yangyxt/indexed_genome/HG002_SD_priority_component_pairs",
                                 target_bed = "",
                                 input_bam = "",
                                 fraction_cutoff = 0.7,
                                 conf_level = 0.05,
                                 threads = 12,
                                 aggregation_resolution = .5,
                                 mq_cutoff = 20,
                                 reference_sd_map = "/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/hg19/test_BISER/WGAC.hg19.cigar.trimmed.homo.expanded.bed",
                                 target_tag = None,
                                 parameter_tuning = False,
                                 extract_read_SD = False):

    try:
        os.mkdir(base_dir)
    except FileExistsError:
        pass


    # Step 0: Calculate the avg frag size and std frag size
    avg_frag_size, std_frag_size = get_bam_frag_size(input_bam)

    logger.info(f"BAM {input_bam} has an average fragment size of {avg_frag_size}")

    # Step 0: Define names of intermediate files
    if not target_tag:
        target_tag = str(os.path.basename(target_bed)).split(".")[0]

    qname_lst = input_bam.replace(".bam", f".{target_tag}.qname.lst")
    bam_json = input_bam.replace(".bam", f".{target_tag}.json")
    tmp_bed = prepare_tmp_file(suffix=".bed").name
    multi_align_bed = input_bam.replace(".bam", f".{target_tag}.multialign.bed")

    # Step 1 : Pick the multi_aligned regions within the target regions and save as a bed file
    if os.path.exists(multi_align_bed) and is_file_up_to_date(multi_align_bed, [input_bam, target_bed, os.path.abspath(__file__)]):
        multi_align_bed_obj = BedTool(multi_align_bed)
    else:
        # Original code to generate multi_align_bed
        multi_align_bed = main_func_pick_region(input_bam, multi_align_bed, MQ_threshold=40, depth_threshold=10, minimum_depth=3, target_region=target_bed, target_tag="FCRs")
        multi_align_bed_obj = BedTool(multi_align_bed)

    multi_align_bed_obj.saveas(tmp_bed)
    logger.info("The targeted bed to slice from {} is {} and it covers {} bp.".format(input_bam,
                                                                                      tmp_bed,
                                                                                      calculate_bed_size(tmp_bed)))

    # Step 2: Slice the BAM by the bed file generated in the step above.
    # Then extract the XA read_pairs with low MAPQ and output the results as bam json file (exclusively from sambamba)
    cmd = f"sambamba slice -q -L {tmp_bed} {input_bam} | \
            sambamba view -q -F \"[SA] == null and mapping_quality <= {mq_cutoff + 10}\" -f sam -t {threads} /dev/stdin | \
            cut -f 1 | sort - | uniq - > {qname_lst} && \
            samtools view -N {qname_lst} -@ {threads} -u {input_bam} | \
            sambamba view -q -f json -o {bam_json} -t {threads} /dev/stdin && \
            ls -lh {bam_json}"

    executeCmd(cmd)

    # Step 3: Extract the SD pair table str from the json file
    sd_map_str = main_parse_json_and_process(bam_json,
                                             avg_frag_size=avg_frag_size,
                                             std_frag_size=std_frag_size,
                                             conf_level=conf_level,
                                             threads=threads,
                                             read_SD_extract=extract_read_SD)

    # Now we need to filter out the geonmic intervals that not in the SD pair list using bedtools intersect.
    regex_pattern = r'(chr)*(X|Y|MT*|[0-9][0-9]*)$'
    sd_map_df = pd.read_csv(StringIO(sd_map_str), low_memory=False, header=None, names=["chr", "start", "end", "chr_counterparts", "start_counterparts", "end_counterparts"]).drop_duplicates()
    all_involved_chrs = set(sd_map_df.loc[:, "chr"].drop_duplicates())
    all_involved_chrs.update(set(sd_map_df.loc[:, "chr_counterparts"].drop_duplicates()))
    logger.info(f"{input_bam} extracted SD maps from intersected regions between multi-aligned regions and target regions specified by {target_bed} contain {sd_map_df.shape[0]} binary SD maps. All the included chromosomes are {all_involved_chrs}")
    sd_map_df = sd_map_df.loc[sd_map_df.loc[:, "chr"].str.match(regex_pattern) & sd_map_df.loc[:, "chr_counterparts"].str.match(regex_pattern), :]
    all_involved_chrs = set(sd_map_df.loc[:, "chr"].drop_duplicates())
    all_involved_chrs.update(set(sd_map_df.loc[:, "chr_counterparts"].drop_duplicates()))
    logger.info(f"After removing the SD relationships involving alternative chromosomes, {input_bam} extracted SD maps from intersected regions between multi-aligned regions and target regions specified by {target_bed} contain {sd_map_df.shape[0]} binary SD maps.  All the included chromosomes are {all_involved_chrs}")
    assert len(all_involved_chrs) <= 25

    # Since the first 3-column set in the sd_map_df is storing the primary alignment region of the reads. And we only care the region of target and poorly aligned region.
    # So we can filter on the first 3 columns region by doing bedtools intersect to pickout
    sd_map_bed = BedTool.from_dataframe(sd_map_df).intersect(multi_align_bed_obj,
                                                             wa=True)
    sd_map_df = sd_map_bed.to_dataframe(disable_auto_names=True, names = ["chr", "start", "end",
                                                                          "chr_counterparts", "start_counterparts", "end_counterparts"]).drop_duplicates()

    # Now we need to expand the SD map df to include the counterparts regions
    reverse_sd_df = sd_map_df.loc[:, ["chr_counterparts", "start_counterparts", "end_counterparts", "chr", "start", "end"]].rename(columns={"chr_counterparts":"chr", "start_counterparts":"start", "end_counterparts":"end", "chr":"chr_counterparts", "start":"start_counterparts", "end":"end_counterparts"})
    merged_sd_df = pd.concat([sd_map_df, reverse_sd_df], axis=0, ignore_index=True).drop_duplicates()
    sd_map_bed = BedTool.from_dataframe(merged_sd_df)
    logger.info("The total targeted SD region size is {}bp according to the BAM file {}. The BAM extracted SD map table looks like:\n{}\n".format(sd_map_bed.sort().merge().total_coverage(),
                                                                                                                                                  input_bam,
                                                                                                                                                  sd_map_df[:10].to_string(index=False)))

    # Step 4: Now we have a bed file recording all the region involved in multialignment by the mapper.
    # We need to intersect this bed file with the reference SD map to get the SDs that really causing mapping ambiguity for this sample
    ref_bed_obj = BedTool(reference_sd_map)
    logger.info("The reference SD map total coverage size is {}bp".format(ref_bed_obj.sort().merge().total_coverage()))
    total_bin_sd_bed = ref_bed_obj.intersect(sd_map_bed,
                                             wo=True)
    total_bin_sd_df = total_bin_sd_bed.to_dataframe(disable_auto_names=True,
                                                    names = ["chr_1", "start_1", "end_1",
                                                             "chr_2", "start_2", "end_2",
                                                             "strand1", "strand2",
                                                             "cigar", "mismatch_rate",
                                                             "chr_bam1", "start_bam1", "end_bam1",
                                                             "chr_bam2", "start_bam2", "end_bam2",
                                                             "overlap_len"]).loc[:,["chr_1", "start_1", "end_1", "strand1",
                                                                                    "chr_2", "start_2", "end_2", "strand2",
                                                                                    "chr_bam1", "start_bam1", "end_bam1"]].drop_duplicates().dropna()
    logger.info("The total binary SD dataframe has shape of {} and it looks like: \n{}\n".format(total_bin_sd_df.shape, total_bin_sd_df[:5].to_string(index=False)))
    both_neg_strand = (total_bin_sd_df.loc[:, "strand1"] == "-") & (total_bin_sd_df.loc[:, "strand2"] == "-")
    reverse_strand_dict = {"-": "+", "+": "-"}
    total_bin_sd_df.loc[both_neg_strand, "strand1"] = total_bin_sd_df.loc[both_neg_strand, "strand1"].map(reverse_strand_dict)
    total_bin_sd_df.loc[both_neg_strand, "strand2"] = total_bin_sd_df.loc[both_neg_strand, "strand2"].map(reverse_strand_dict)

    # Now we need to filter out the intervals on the alternative contigs
    main_contigs = total_bin_sd_df.loc[:, "chr_1"].str.match(regex_pattern) & \
                   total_bin_sd_df.loc[:, "chr_2"].str.match(regex_pattern) & \
                   total_bin_sd_df.loc[:, "chr_bam1"].str.match(regex_pattern)

    total_bin_sd_df = total_bin_sd_df.loc[main_contigs, :]
    total_bin_sd_df.to_csv(os.path.join(base_dir, "raw_SD_binary_map.tsv"), sep="\t", index=False)

    # Then we need to filter out some SDs overlapped with the same XA region
    # Given an XA region, multiple SDs overlapped with them. So just keep the minimal SD interval with sufficient overlap
    by_bam_region = total_bin_sd_df.groupby(["chr_bam1", "start_bam1", "end_bam1"], as_index=False)
    groups = [group for _, group in by_bam_region]
    with Pool(threads) as pool:
        results = pool.imap_unordered(select_minimal_sd_span, groups)
        total_bin_sd_df = pd.concat(results, axis=0, ignore_index=True).loc[:, ["chr_1", "start_1", "end_1", "strand1",
                                                                                "chr_2", "start_2", "end_2", "strand2"]].drop_duplicates().dropna()

    logger.info("Now the total binary SD table looks like this:\n{}\n\n".format(total_bin_sd_df[:5].to_string(index=False)))

    # Remove the duplicated rows (Here duplicates means the same SDs but with different order of intervals)
    total_bin_sd_df.loc[:, "frozenset_indx"] = total_bin_sd_df.apply(lambda row: frozenset({ row["chr_1"]+":"+str(row["start_1"])+"-"+str(row["end_1"])+":"+row["strand1"], row["chr_2"]+":"+str(row["start_2"])+"-"+str(row["end_2"])+":"+row["strand2"]}), axis=1)
    total_bin_sd_df = total_bin_sd_df.drop_duplicates(subset="frozenset_indx").drop(columns=["frozenset_indx"])
    total_bin_sd_df.to_csv(os.path.join(base_dir, "filtered_SD_binary_map.tsv"), sep="\t", index=False)
    logger.info("After intersecting with known SDs, the total targeted SD region size for sample {} is {}bp, The total binary SD dataframe has shape of {} and it looks like: \n{}\n".format(input_bam,
                                                                                                                                                                                             total_bin_sd_bed.sort().merge().total_coverage(),
                                                                                                                                                                                             total_bin_sd_df.shape,
                                                                                                                                                                                             total_bin_sd_df[:5].to_string(index=False)))

    # Step 5: Create a graph on the total_bin_sd_df
    graph_path = os.path.join(base_dir, "multiplexed_homologous_sequences.graphml")
    graph = create_multiplex_graph( total_bin_sd_df, graph_path,
                                        overlap_frac_min = fraction_cutoff,
                                        avg_frag_size = avg_frag_size,
                                        std_frag_size = std_frag_size,
                                        resolution = aggregation_resolution,
                                        threads = threads,
                                        target_bed = target_bed,
                                        base_dir = base_dir )

    # Step 6: Pick out the nodes that overlaps with target region, these intervals come from filtered reference SD map. Intervals might tangle(overlap) with each other
    graph_bed = graph_to_sorted_bedtool(graph)
    # Filter out the graph bed where the interval length should be larger than 150 bp (common read length)
    graph_bed = graph_bed.filter(lambda l: len(l) >= 140)

    # This is to identify which graph nodes overlap with the target regions to pool the reads to
    intersect_df = graph_bed.intersect(multi_align_bed_obj.sort().merge(), wo=True).to_dataframe(disable_auto_names=True, names=["chr", "start", "end", "strand",
                                                                                                                                 "chr_query", "start_query", "end_query", "overlap_len"])
    intersect_df = intersect_df.loc[(intersect_df["end"] > intersect_df["start"]) & (intersect_df["end_query"] > intersect_df["start_query"]), :]

    # Here we can filter out the FC-overlapping intervals that are completely enwrapped by other FC-overlapping intervals
    fc_intervals = intersect_df.loc[:, ["chr", "start", "end", "strand"]].drop_duplicates()
    fc_interval_bed = BedTool.from_dataframe(fc_intervals)
    fc_internal_overlap = fc_interval_bed.intersect(fc_interval_bed, wao=True).to_dataframe(disable_auto_names=True, names=[ "chr1", "start1", "end1", "strand1",
                                                                                                                       "chr2", "start2", "end2", "strand2",
                                                                                                                       "overlap_len" ])
    enwrapped_fcs = set([])
    for idx, row in fc_internal_overlap.iterrows():
        interval_1_size = int(row["end1"]) - int(row["start1"])
        interval_2_size = int(row["end2"]) - int(row["start2"])
        interval_1 = (row["chr1"], int(row["start1"]), int(row["end1"]))
        interval_2 = (row["chr2"], int(row["start2"]), int(row["end2"]))
        if interval_1 == interval_2:
            continue
        if row["overlap_len"] == interval_1_size:
            enwrapped_fcs.add(interval_1)
        elif row["overlap_len"] == interval_2_size:
            enwrapped_fcs.add(interval_2)
    logger.info("There are {} FC-overlapping intervals that are enwrapped by other FC-overlapping intervals".format(len(enwrapped_fcs)))

    by_query_interval = intersect_df.groupby(["chr_query", "start_query", "end_query"], as_index=False)
    query_groups = [group for _, group in by_query_interval]
    with Pool(threads) as pool:
        results = pool.imap_unordered(imap_pickup_regions_per_exon, zip(query_groups,
                                                                        repeat("chr_query"),
                                                                        repeat("start_query"),
                                                                        repeat("end_query"),
                                                                        repeat(avg_frag_size),
                                                                        repeat(std_frag_size)))
        intersect_df = pd.concat(results, axis=0, ignore_index=True).drop_duplicates()
    logger.info("These intervals overlaps with the filtered SD intervals: \n{}\n".format(intersect_df[:20].to_string(index=False)))

    query_nodes = intersect_df.loc[:, ["chr", "start", "end", "strand"]].drop_duplicates()
    final_query_bed = BedTool.from_dataframe(query_nodes)
    final_query_bed.intersect(BedTool(tmp_bed), wo=True).saveas(os.path.join(base_dir, "coding_query_nodes.bed"))
    logger.info("Save the coding overlap query nodes to this file: {}\n".format(os.path.join(base_dir, "coding_query_nodes.bed")))

    logger.info("Finally after picking up the merged the SD intervals, there are {} merged SDs from graph overlapped with exon, and it covers a region of {}bp.".format(query_nodes.shape[0],
                                                                                                                                                                        final_query_bed.sort().merge().total_coverage()))

    # Step 7. Now use a function to generate FC and NFC pairs for all the query nodes
    final_graph_path = graph_path.replace(".graphml", ".trim.annoPC.graphml")
    fc_nfc_pairs = extract_FC_NFC_pairs_from_graph(query_nodes, graph, graph_path = final_graph_path, avg_frag_size = avg_frag_size, std_frag_size = std_frag_size, threads=threads)
    logger.info("After extracting the FC-NFC pairs, we got {} FC-NFC pairs".format(len(fc_nfc_pairs)))

    # Now we need to test which FC-NFC pairs can be merged together (meaning FCs cannot appear in NFCs covered region)
    # Step 8: Pass the query nodes to function and get the cluster of SD intervals (The query nodes here is a pandas dataframe object)
    '''
    2023-12-15:
    Lets try not to merge them. This way wo can have better control of the remapping of the reads.
    This way we can also have better estimation of the possible haplotypes.

    2023-12-18:
    Failed to use separate FC-NFC pairs due to the number of them.
    '''
    # fc_nfc_dict = {n:ns for n, ns in fc_nfc_pairs}
    # all_results = [ {"PCs": [n], "SD_counterparts": ns} for n, ns in fc_nfc_pairs ]
    # disconnected_result = all_results[0]
    # connected_results = all_results[1:]
    pair_graph_path = graph_path.replace(".graphml", ".pair.graphml")
    disconnected_result, connected_results, fc_nfc_dict = query_connected_nodes(fc_nfc_pairs,
                                                                                multiplex_graph = graph,
                                                                                threads = threads,
                                                                                avg_frag_size = avg_frag_size,
                                                                                std_frag_size = std_frag_size,
                                                                                graph_path = pair_graph_path)

    # Step 8.5: We need to tag the nodes with FC-NFC pair ID
    all_results = [disconnected_result] + connected_results
    final_graph = graph.copy()
    for i, result in enumerate(all_results):
        tag = "PC" + str(i)
        for node in result["PCs"]:
            # Add an attribute "FC" to the node
            final_graph.nodes[node]["FC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("FC", "").split(",") + [tag])) if len(e) > 0])
        for node in result["SD_counterparts"]:
            # Add an attribute "NFC" to the node
            assert type(node) == HOMOSEQ_REGION, f"The counterpart node {node} type is not HOMOSEQ_REGION, the FC nodes are {result['PCs']}"
            final_graph.nodes[node]["NFC"] = ",".join([e for e in list(dict.fromkeys(final_graph.nodes[node].get("NFC", "").split(",") + [tag])) if len(e) > 0])
    nx.write_graphml(final_graph, final_graph_path)


    # Step 9: Create beds and masked genomes
    convert_nodes_into_hierachical_beds(disconnected_result = disconnected_result,
                                        connected_results = connected_results,
                                        fc_nfc_dict = fc_nfc_dict,
                                        output_folder = base_dir,
                                        ref_genome = ref_genome,
                                        nthreads = threads,
                                        avg_frag_size = avg_frag_size,
                                        std_frag_size = std_frag_size)


