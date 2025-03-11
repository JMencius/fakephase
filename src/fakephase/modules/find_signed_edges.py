import logging
import pysam


def is_snv_in_read(read, target_pos, target_snv) -> int:
    pass


def is_indel_in_read(read, target_pos, target_indel) -> int:
    pass



def process_read(read) -> dict:
    pass


def find_signed_edges(bamfile: str, maxlen: int, min_mapq: int, variants: dict, ref_len: dict, working_chr: str) -> dict:
    edges_dict = dict()
    
    # region in pysam is [start, end)
    regions = [(0, maxlen + 1), (ref_len[working_chr] - maxlen - 1, ref_len[working_chr])])]
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for start, end in regions:
            for read in bam.fetch(working_chr, start, end):
                # filter unmapped read
                if not read.is_unmapped:
                    if read.mapping_quality > min_mapq:
                        

    return edges_dict
