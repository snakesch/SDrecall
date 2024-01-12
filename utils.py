def get_bam_frag_size(input_bam):
    
    import subprocess
    
    p = subprocess.run(f"samtools stats {input_bam} | grep ^SN", shell=True, capture_output=True, text=True)
    for line in p.stdout.split('\n'):
        if line.startswith('SN'):
            fields = line.split('\t')
            if fields[1].startswith('insert size average:'):
                avg_insert_size = float(fields[2])
            elif fields[1].startswith('insert size standard deviation:'):
                std_insert_size = float(fields[2])
    return avg_insert_size, std_insert_size