def getIntersect(all_blocks: list):
    '''
    This function takes in a list of lists and groups sub-lists if their intersects is not NULL.
    '''
    # Base case
    if len(all_blocks) == 0:
        return None
    
    intersect = []
    tmp = []
    base = set()
    for cur in range(len(all_blocks)):
        # Does it intersect?
        if base.intersection(all_blocks[cur]) == set() and cur != 0:
            intersect.append(tuple(tmp))
            tmp.clear()
            tmp.append(all_blocks[cur])
            base = set(all_blocks[cur])
        else:
            tmp.append(all_blocks[cur])
            # Find new members
            for _idx in all_blocks[cur]:
                if _idx not in base:
                    base.add(_idx)
    intersect.append(tuple(tmp))
    
    return intersect

def read_vcf_gz(path=""):
    '''
    Returning a tuple
    Return vcf head as a list of lines(each line is a str object) as the first element of the tuple
    Return vcf body as pandas DataFrame Object as the second element of the tuple
    '''
    with bgzf.open(path, "r") as bgzf_handle:
        line = ""
        header_lines = []
        while True:
            reader_pos = bgzf_handle.tell()
            line = bgzf_handle.readline()
            if re.search(r'#CHROM', line):
                break
            header_lines.append(line)
            
        headers = line.strip("\n").split("\t")
        logging.info("This vcf ({}) body content start at line {}".format(path, reader_pos))
        logging.info("This vcf ({}) body header are {}".format(path, headers))

        total_records = []
        while True:
            reader_pos = bgzf_handle.tell()
            logging.debug("Current reader pos is {}".format(reader_pos))
            line = bgzf_handle.readline()
            total_records.append(tuple(line.strip('\n').split("\t"))) 
            if reader_pos == bgzf_handle.tell():
                del total_records[-1]
                break
        bgzf_handle.close()
        return (header_lines, pd.DataFrame(total_records, columns=headers))


def read_vcf_gz_chunkwise(path="", chunksize=10000):
    '''
    Return an iterator of tuples
    Return vcf head as a list of lines(each line is a str object) as the first element of the tuples (the same content across all the tuples)
    Return vcf body as pandas DataFrame Object as the second element of the tuple (chunked dfs)
    '''
    with bgzf.open(path, "r") as bgzf_handle:
        line = ""
        header_lines = []
        while True:
            reader_pos = bgzf_handle.tell()
            line = bgzf_handle.readline()
            header_lines.append(line)
            if re.search(r'#CHROM', line):
                break
        headers = line.strip("\n").split("\t")
        logging.info("This vcf ({}) body content start at pos {}".format(path, reader_pos))
        logging.info("This vcf ({}) body header are {}".format(path, headers))
        
        while True:
            chunk_records = []
            for i in range(0, chunksize):
                reader_pos = bgzf_handle.tell()
                chunk_records.append(tuple(bgzf_handle.readline().strip("\n").split("\t")))
                if reader_pos == bgzf_handle.tell():
                    break
            if reader_pos == bgzf_handle.tell():
                del chunk_records[-1]
            if len(chunk_records) > 0:
                yield (header_lines, pd.DataFrame(chunk_records, columns=headers))
            else:
                break
