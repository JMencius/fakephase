import logging
import pysam
import sys


def is_snv_in_read(read, target_pos: int, variant) -> int:
    """
    return 0 means not in read
    return -1 means left match
    return 1 means right match
    """
    variant_pool = list(variant.ref) + variant.alt
    read_base = read.query_sequence[target_pos]
    if read_base == variant_pool[variant.left]:
        return -1
    elif read_base == variant_pool[variant.right]:
        return 1
    else:
        return 0



def get_cigar_type(ref_pos: int, target_pos: int, cigartuples: list) -> str:
    # index for inside read
    query_pos = 0

    for cigartype, length in cigartuples:
        # Matches (M)
        if cigartype == 0:
            if ref_pos <= target_pos < ref_pos + length:
                return 'M'
            ref_pos += length
            query_pos += length

        # Insertion (I)
        elif cigartype == 1:
            if ref_pos == target_pos:
                return 'I'
            query_pos += length
        
        # Deletion (D)
        elif cigartype == 2:
            if ref_pos <= target_pos < ref_pos + length:
                return 'D'
            ref_pos += length
    
        # skip for soft/hard clip
        elif cigartype in (4, 5):
            pass
    
    # Not found
    return None



def is_indel_in_read(read, target_pos: int, variant) -> int:
    """
    return 0 means not in read
    return -1 means left match
    return 1 means right match
    """
    category = get_cigar_type(read.reference_start, target_pos, read.cigartuples)
    if not category:
        logging.error("Target position {target_pos} not in {read.query_name}")
        sys.exit(1)
    
    variant_pool = list(variant.ref) + variant.alt   
    if category == 'M':
        if variant.left == 0:
            return -1
        if variant.right == 0:
            return 1
    
    else:
        left_seq = variant_pool[variant.left]
        insert_seq = read.query_sequence[target_pos : target_pos + len(left_seq)]
        if left_seq == insert_seq:
            return -1        

        right_seq = variant_pool[variant.right]
        insert_seq = read.query_sequence[target_pos : target_pos + len(right_seq)]
        if right_seq == insert_seq:
            return 1
    


def process_pairs(edges_list: list) -> list:
    processed_pairs = list()
    for i in range(len(edges_list) - 1):
        for j in range(i + 1, len(edges_list)):
            a = edges_list[i][0]
            b = edges_list[j][0]
            sign = edges_list[i][1] * edges_list[j][1]
            processed_pairs.append((min(a, b), max(a, b), sign))
        
    return processed_pairs



def process_read(read, variants: dict) -> list:
    edges_list = list()
    for site in range(read.reference_start, read.reference_end + 1):
        if site in variants:
            subject_variant = variants[site]
            if subject_variant.category == "SNV":
                edges_list.append((site, is_snv_in_read(read, site, variants)))
            if subject_variant.category == "INDEL":
                edges_list.append((site, is_indel_in_read(read, site, variants)))               
    
    signed_edges_in_read = process_paris(edges_list)
    
    return signed_edges_in_read



def find_signed_edges(bamfile: str, maxlen: int, min_mapq: int, variants: dict, ref_len: dict, working_chr: str) -> dict:
    # region in pysam is [start, end)
    if (maxlen + 1) < (ref_len[working_chr] - maxlen - 1):
        regions = [(0, maxlen + 1), (ref_len[working_chr] - maxlen - 1, ref_len[working_chr])]
    else:
        regions = [(0, ref_len[working_chr])]    
    

    edges_dict = dict()
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for start, end in regions:
            for read in bam.fetch(working_chr, start, end):
                # filter unmapped read
                if not read.is_unmapped:
                    # filter based on mapping quality threshold
                    if read.mapping_quality > min_mapq:
                        edges_in_read = process_read(read, variants)
                        for edg in edges_in_read
                            if (edg[0], edg[1]) not in edges_dict:
                                edges_dict[(edg[0], edg[1])] = [0, 0]
                            if edg[2] == 1:
                                edges_dict[(edg[0], edg[1])][0] += 1
                            else:
                                edges_dict[(edg[0], edg[1])][1] += 1

    return edges_dict


