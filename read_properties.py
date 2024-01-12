def get_fragment_dist(bam_file):
    
    import subprocess
    
    proc = subprocess.run(["samtools", "stats", bam_file], text=True, capture_output=True)

    for line in proc.stdout.split('\n'):
        if line.startswith('SN'):
            fields = line.split('\t')
            if fields[1].startswith('insert size average:'):
                avg_insert_size = float(fields[2])
            elif fields[1].startswith('insert size standard deviation:'):
                std_insert_size = float(fields[2])
        elif line.startswith("FF"):
            break
    return avg_insert_size, std_insert_size
