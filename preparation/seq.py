from tempfile import NamedTemporaryFile

from src.utils import executeCmd

def getRawseq(bedf_path: str, fastq_path: str, ref_genome: str, padding = 0):
    if padding > 0:
        with NamedTemporaryFile(dir = "/tmp") as tmp_bedf:
            tmp_bedf_path = tmp_bedf.name
            cmd = f"bedtools slop -i {bedf_path} -g {ref_genome}.fai -b {padding} > {tmp_bedf_path}"
            executeCmd(cmd)
            cmd = f"seqtk subseq -s {ref_genome} {tmp_bedf_path} | \
                    seqtk seq -F 'F' - | \
                    mawk 'FNR %4 == 2{{ print toupper($0); }} FNR % 4 !=2 {{ print; }}' > {fastq_path}"
            executeCmd(cmd)
    else:
        cmd = f"seqtk subseq -s {ref_genome} {bedf_path} | \
                seqtk seq -F 'F' - | \
                mawk 'FNR %4 == 2{{ print toupper($0); }} FNR % 4 !=2 {{ print; }}' > {fastq_path}"     
        executeCmd(cmd)
    return fastq_path

def get_bam_frag_size(input_bam):
    cmd = f"samtools stats {input_bam}"
    output = executeCmd(cmd)
    for line in output.split('\n'):
        if line.startswith('SN'):
            fields = line.split('\t')
            if fields[1].startswith('insert size average:'):
                avg_insert_size = float(fields[2])
            elif fields[1].startswith('insert size standard deviation:'):
                std_insert_size = float(fields[2])
    return avg_insert_size, std_insert_size

