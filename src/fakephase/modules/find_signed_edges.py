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


def is_snv_in_read(read, pos, current_variant, aligned_pairs) -> int:
    """
    return 0 means not in read
    return -1 means left match
    return 1 means right match
    """
    if not read.is_reverse:
        query_sequence = read.query_sequence
    else:
        query_sequence = reverse_complement(read.query_sequence)
    

    for read_pos, ref_pos in aligned_pairs:
        if ref_pos == pos - 1:  # pysam 是 0-based，pos 是 1-based
            query_base = query_sequence[read_pos]

            variant_pool = list(current_variant.ref) + current_variant.alt
            if variant_pool[current_variant.left] == query_base:
                return -1
            elif variant_pool[current_variant.right] == query_base:
                return 1
            else:
                return None      

    return None  # SNV 不在 read 中




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
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    for site in range(read.reference_start, read.reference_end):
        if site in variants:
            subject_variant = variants[site]
            if subject_variant.category == "SNV":
                pos = is_snv_in_read(read, site, variants[site], aligned_pairs)
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

    logging.info(f"{working_chr} signed edges extract finished")    

    return edges_dict


