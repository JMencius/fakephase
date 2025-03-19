import pyfastx
import logging

def find_first_last_non_N(sequence: str) -> tuple:
    first_non_N = None
    last_non_N = None
    
    for i, base in enumerate(sequence):
        if base != 'N':
            first_non_N = i
            break
    
    for i in range(len(sequence)-1, -1, -1):
        if sequence[i] != 'N':
            last_non_N = i
            break
    
    return (first_non_N, last_non_N)



def find_start_pos(reference_file: str, chrom: list) -> list:
    result_dict = dict()
    for name, seq in pyfastx.Fasta(reference_file, build_index = False):
        if name in chrom:
            result_dict[name] = find_first_last_non_N(seq)
    
    if len(result_dict.keys()) != len(chrom):
        for i in chrom:
            if i not in result_dict.keys():
                logging.error(f"{i} not in reference file")
        sys.exit(1)
    else:
        start_end_list = list()
        for i in chrom:
            start_end_list.append(result_dict[i])

    return start_end_list
