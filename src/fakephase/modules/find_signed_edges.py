import logging
import pysam
import sys
from multiprocessing import Pool


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
        if ref_pos == pos - 1:  # pysam is 0-based，pos is 1-based
            query_base = query_sequence[read_pos]

            variant_pool = list(current_variant.ref) + current_variant.alt
            if variant_pool[current_variant.left] == query_base:
                return -1
            elif variant_pool[current_variant.right] == query_base:
                return 1
            else:
                return None      

    return None 




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



def find_signed_edges(bamfile: str, maxratio: float, min_mapq: int, variants: dict, ref_len: dict, chrom: list, ref: str, start_end: dict, chunk_size: int, threads: int) -> dict:
    # only phase to a proportion of the long arm
    # The phase length is defined by `maxratio`
    to_process_regions = list()
    count = 0
    for c in chrom:
        s, e = start_end[count]
        phase_len = ref_len[c] * maxratio

        start_site = e - phase_len - 1
        while start_site + chunk_size < (e):
            to_process_regions.append((c, start_site, start_site + chunk_size))
            start_site += chunk_size
        if start_site != e:
            to_process_regions.append((c, start_site, e))

        count += 1

    logging.info(f"In total {len(to_process_regions)} chunks to process")

    with Pool(threads) as p:
        temp_results = p.starmap(worker, [(bamfile, maxratio, variants[chrom.index(i[0])], i[0], ref, (i[1], i[2])) for i in to_process_regions])

    assembled_results = [dict() for i in range(len(chrom))]
    for i in range(len(to_process_regions)):
        idx = chrom.index(to_process_regions[i][0])
        for e, f in temp_results[i].items():
            if e not in assembled_results[idx]:
                assembled_results[idx][e] = f
            else:
                assembled_results[idx][e] += f            

    return assembled_results



def worker(bamfile: str, min_mapq: int, variants: dict, working_chr: str, ref: str, start_end: tuple) -> dict:
    edges_dict = dict()
    max_bytes = int(1.8 * (2**30))

    start, end = start_end

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        if True:
            for read in bam.fetch(working_chr, start, end):
                # filter unmapped, secondary, supplementary read
                if (not read.is_unmapped) and (not read.is_secondary) and (not read.is_supplementary):
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
                
                if sys.getsizeof(edges_dict) > max_bytes:
                    logging.warning(f"{working_chr} {start}-{end} signed edges dict exceeded {max_bytes} bytes, stopping early")
                    return edges_dict
    

    return edges_dict





