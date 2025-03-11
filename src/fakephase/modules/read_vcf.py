import logging
from cyvcf import VCF
import sys
from fakephase.classes.variant import myvariant


def get_category(ref: str, alt: list) -> str:
    max_len = max(len(ref), max([len(i) for i in alt]))
    if max_len == 1:
        ctype = "SNV"
    elif 1 < max_len <= 30:
        ctype = "INDEL"
    else:
        ctype = "SV"

    return ctype


def read_vcf(filename: str, working_chr: str, no_indel: bool) -> dict:
    chr_variants = dict()
    for variant in VCF(filename, threads = 1):
        if variant.CHROM == working_chr:
            left, right, _ = variant.genotypes[0]
            # only read heterozyous site
            if (left != right) and ('.' not in variant.gt_bases[0]):
                category = get_category(variant.REF, variant.ALT)
                # set filter
                flag = 0
                if no_indel:
                    if category == "SNV":
                        flag = 1
                else:
                    if (category == "SNV") or (category == "INDEL"):
                        flag = 1

                if flag:
                    current_variant = myvariant(variant.CHROM, variant.POS, variant.REF, variant.ALT, left, right, category)
                    if variant.POS not in chr_variants:
                        chr_variant[variant.POS] = current_variant
                    else:
                        logging.WARNING(f"Site overlap in {variant.CHROM} at {variant.POS}, discard the second variant at this site")
    
    return chr_variants
