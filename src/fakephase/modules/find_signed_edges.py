import logging
import pysam
import sys


def reverse_complement(seq: str) -> str:
    if seq == '':
        return ''

    complement = str.maketrans(
        "ACGTacgtRYMKBDHVSWNrymkbdhvswn",  # original
        "TGCAtgcaYRKMVHDBSWNyrkmvhdbswn"   # complement
    )
    result_seq = (seq[::-1]).translate(complement)
    
    return result_seq


def is_snv_in_read(read, pos, current_variant, ref) -> int:
    """
    return 0 means not in read
    return -1 means left match
    return 1 means right match
    """
    if not read.is_reverse:
        query_sequence = read.query_sequence
    else:
        query_sequence = reverse_complement(read.query_sequence)
    

    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos == pos - 1:  # pysam 是 0-based，pos 是 1-based
            query_base = query_sequence[read_pos]

            variant_pool = list(current_variant.ref) + current_variant.alt
            if variant_pool[current_variant.left] == query_base:
                return -1
            elif variant_pool[current_variant.right] == query_base:
                return 1           

    return None  # SNV 不在 read 中



def is_variant_in_read(read, target_pos: int, variant, ref, cigar_dict) -> int:
    """
    return 0 means not in read
    return -1 means left match
    return 1 means right match
    """
    category = cigar_dict[target_pos]
    if not category:
        logging.error(f"Target position {target_pos} not in {read.query_name}, read starts at {read.reference_start}, span with Reference Length: {read.reference_length}")
        sys.exit(1)
    
    if not read.is_reverse:
        query_sequence = read.query_sequence
    else:
        query_sequence = reverse_complement(read.query_sequence)

    variant_pool = list(variant.ref) + variant.alt   
    if category == 'M':
        if variant.left == 0:
            return -1
        if variant.right == 0:
            return 1
    elif category == 'X':
        if variant.left != 0:
            return -1
        if variant.right != 0:
            return 1    

    else:
        left_seq = variant_pool[variant.left]
        insert_seq = query_sequence[target_pos : target_pos + len(left_seq)]
        if left_seq == insert_seq:
            return -1        

        right_seq = variant_pool[variant.right]
        insert_seq = query_sequence[target_pos : target_pos + len(right_seq)]
        if right_seq == insert_seq:
            return 1
    


def process_pairs(edges_list: list) -> list:
    processed_pairs = list()
    if len(edges_list) <= 1:
        return list()
    for i in range(len(edges_list) - 1):
        for j in range(i + 1, len(edges_list)):
            a = edges_list[i][0]
            b = edges_list[j][0]
            sign = edges_list[i][1] * edges_list[j][1]
            processed_pairs.append((min(a, b), max(a, b), sign))
        
    return processed_pairs



def process_read(read, variants: dict, ref) -> list:
    edges_list = list()

    count = 0
    for site in range(read.reference_start, read.reference_end):
        if site in variants:
            subject_variant = variants[site]
            if subject_variant.category == "SNV":
                ##print(read.reference_start, site, read.reference_end)
                pos = is_snv_in_read(read, site, variants[site], ref)
                if pos:
                    edges_list.append((site, pos))
            if subject_variant.category == "INDEL":
                pos = is_snv_in_read(read, site, variants[site], ref)
                if pos:
                    edges_list.append((site, pos))               
            

    signed_edges_in_read = process_pairs(edges_list)
    
    return signed_edges_in_read



def find_signed_edges(bamfile: str, maxlen: int, min_mapq: int, variants: dict, ref_len: dict, working_chr: str, ref: str, start_end: tuple) -> dict:
    # region in pysam is [start, end)
    s, e = start_end
    regions = [(s, s + maxlen + 1), (e - maxlen - 1, e)]

    edges_dict = dict()
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for start, end in regions:
            for read in bam.fetch(working_chr, start, end):
                # filter unmapped read
                if (not read.is_unmapped) and (not read.is_secondary):
                    # filter based on mapping quality threshold
                    if read.mapping_quality > min_mapq:
                        edges_in_read = process_read(read, variants, ref)
                        for edg in edges_in_read:
                            if (edg[0], edg[1]) not in edges_dict:
                                edges_dict[(edg[0], edg[1])] = [0, 0]
                            if edg[2] == 1:
                                edges_dict[(edg[0], edg[1])][0] += 1
                            else:
                                edges_dict[(edg[0], edg[1])][1] += 1

    return edges_dict


