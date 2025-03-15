import pyfastx
import logging
import sys


# This function operate in single thread
def get_ref_len(reference_file: str, chrom: list) -> list:
    chr_len = dict()
    exist = set()
    for name, seq in pyfastx.Fasta(reference_file, build_index = False):
        if name in chrom:
            chr_len[name] = len(seq)
            exist.add(name)

    if len(chr_len.keys()) != len(chrom):
        for i in chrom:
            if i not in exist:
                logging.error(f"{i} not in reference file")
        sys.exit(1)
    else:
        return chr_len




