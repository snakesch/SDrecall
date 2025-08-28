import pysam
from collections import defaultdict

from src.utils import prepare_tmp_file, executeCmd
from src.log import logger

def get_supporting_tags(bam_file, chrom, pos, ref, alts, tag, min_mapq=10, min_bq=15):
    supporting_tags = defaultdict(set)
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for pileup_column in samfile.pileup(chrom, pos-1, pos, truncate=True, min_base_quality=min_bq, ignore_orphans=False):
            for pileup_read in pileup_column.pileups:
                if pileup_read.alignment.mapping_quality < min_mapq:
                    logger.debug(f"Skipping read {pileup_read.alignment.query_name}:{pileup_read.alignment.flag} with mapping quality {pileup_read.alignment.mapping_quality}")
                    continue

                if pileup_read.is_del:
                    logger.debug(f"Found deletion at {chrom}:{pos} on read {pileup_read.alignment.query_name}:{pileup_read.alignment.flag}")
                    # Deletion
                    if len(ref) > 1 and ref[1:] in alts:
                        alt = ref[1:]
                        if pileup_read.alignment.has_tag(tag):
                            tag_value = pileup_read.alignment.get_tag(tag)
                            supporting_tags[alt].add(tag_value)
                elif pileup_read.indel > 0:
                    logger.debug(f"Found insertion at {chrom}:{pos} on read {pileup_read.alignment.query_name}:{pileup_read.alignment.flag}")
                    # Insertion
                    ins_seq = pileup_read.alignment.query_sequence[pileup_read.query_position:pileup_read.query_position + pileup_read.indel]
                    for alt in alts:
                        if alt.startswith(ref[0]) and ins_seq in alt:
                            if pileup_read.alignment.has_tag(tag):
                                tag_value = pileup_read.alignment.get_tag(tag)
                                supporting_tags[alt].add(tag_value)
                            break
                elif not pileup_read.is_refskip:
                    # SNV
                    read_base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    logger.debug(f"Found SNV at {chrom}:{pos} with ALT base {read_base} on read {pileup_read.alignment.query_name}:{pileup_read.alignment.flag}")
                    if read_base in alts:
                        if pileup_read.alignment.has_tag(tag):
                            tag_value = pileup_read.alignment.get_tag(tag)
                            supporting_tags[read_base].add(tag_value)
                else:
                    logger.debug(f"Skipping read {pileup_read.alignment.query_name}:{pileup_read.alignment.flag} at {chrom}:{pos} due to is_refskip")

    return supporting_tags



def annotate_vcf(input_vcf, output_vcf, bam_file, tag, min_mapq=10, min_bq=15, logger=logger):
    assert output_vcf.endswith(".vcf.gz"), "Output VCF must be a .vcf.gz file"
    tmp_output = prepare_tmp_file(suffix = ".vcf").name
    with pysam.VariantFile(input_vcf) as vcf_in:
        # Add new FORMAT field for the tag
        vcf_in.header.formats.add(f'{tag}SUP', 'A', 'String', f'Supporting {tag} values for each variant allele')
        
        with pysam.VariantFile(tmp_output, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                supporting_tags = get_supporting_tags(bam_file, record.chrom, record.pos, record.ref, record.alts, tag, min_mapq=min_mapq, min_bq=min_bq)
                
                # Format the supporting tags for VCF output
                tag_strings = []
                for alt in record.alts:
                    tag_string = ';'.join(map(str, supporting_tags[alt])) if alt in supporting_tags else '.'
                    tag_strings.append(tag_string)
                
                # Add the new field to the record
                record.samples[0][f'{tag}SUP'] = tag_strings
                
                vcf_out.write(record)

    cmd = f"bcftools sort -Oz -o {output_vcf} {tmp_output} && \
            tabix -f -p vcf {output_vcf} && \
            rm {tmp_output}"
    executeCmd(cmd, logger=logger)
    logger.info(f"Annotation {tag} complete. Output written to {output_vcf}")
    return output_vcf


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate VCF with BAM tag information")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")
    parser.add_argument("-t", "--tag", required=True, help="BAM tag to extract (e.g., 'HP')")
    parser.add_argument("-m", "--min_mapq", type=int, default=10, help="Minimum mapping quality for reads")
    parser.add_argument("-b", "--min_bq", type=int, default=15, help="Minimum base quality for reads")
    args = parser.parse_args()

    annotate_vcf(input_vcf = args.vcf, 
                 output_vcf = args.output, 
                 bam_file = args.bam, 
                 tag = args.tag, 
                 min_mapq = args.min_mapq, 
                 min_bq = args.min_bq)
