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

def calculate_tlen_per_pair(read, mate): 
    '''
    This function expects two tuples in the format ('chrX', '-1776046', '148M', '3'), ('chrX', '-1776046', '148M', '3') and computes the TLEN and the genomic span.
    '''
    import re
    
    read_chr, mate_chr = read[0], mate[0]
    if read_chr != mate_chr:
        return 0, ""
    
    read_pos, mate_pos = int(read[1]), int(mate[1])
    if read_pos * mate_pos >= 0:
        return 0, "" # Both reads mapped to the same strand?
    
    genomic_span = read_chr + ":"
    
    ## read_pos and mate_pos correspond to the leftmost nt position
    if read_pos < mate_pos:
        positive_strand_leftmost = mate_pos
        negative_strand_leftmost = - read_pos
        if negative_strand_leftmost >= positive_strand_leftmost:
            tlen = negative_strand_leftmost - positive_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', read[2]))
            genomic_span += str(positive_strand_leftmost) + "-" + str(positive_strand_leftmost + tlen)
        else: # negative_strand_leftmost < positive_strand_leftmost
            tlen = positive_strand_leftmost - negative_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', mate[2]))
            genomic_span += str(negative_strand_leftmost) + "-" + str(negative_strand_leftmost + tlen)
    else: # read_pos > mate_pos
        positive_strand_leftmost = read_pos
        negative_strand_leftmost = - mate_pos
        if negative_strand_leftmost >= positive_strand_leftmost:
            tlen = negative_strand_leftmost - positive_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', mate[2]))
            genomic_span += str(positive_strand_leftmost) + "-" + str(positive_strand_leftmost + tlen)
        else: # negative_strand_leftmost < positive_strand_leftmost
            tlen = positive_strand_leftmost - negative_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', read[2]))
            genomic_span += str(negative_strand_leftmost) + "-" + str(negative_strand_leftmost + tlen)
    
    return tlen, genomic_span
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
