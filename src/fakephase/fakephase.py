import click
import sys
import os
import logging
from multiprocessing import Pool
from fakephase.modules.read_vcf import read_vcf
from fakephase.modules.get_ref_len import get_ref_len


@click.command()
@click.option("-i", "--invcf", required = True, help = "Input vcf/vcf.gz file for haplotype phasing")
@click.option("-b", "--bam", required = True, help = "Input sorted bam file for haplotype phasing")
@click.option("-r", "--ref", required = True, help = "Reference file fasta")
@click.option("-o", "--output", required = True, help = "Output vcf file")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads [default: 24]")
@click.option("--maxlen", default = 10**6, help = "Maximum length to process [default: 10**6]")
@click.option("--min-mapq", default = 0, help = "Minimum mapping quality (MAPQ)")
@click.option("--chrom", default = ','.join(["chr" + str(i) for i in range(1, 23)]), help = "Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default: chr1-22]")
@click.option("--no-indel", is_flag = True, help = "Do not phase insertion and deletion (INDEL)")
@click.option("--verbose", is_flag = True, help = "Verbose mode print progess to standard output")
@click.version_option(version = "0.1.0", prog_name = r"Fake haplotype phasing algorithm deceive most evalation metrics, based on Python 3.10+")
def main(invcf, bam, ref, output, threads, maxlen, min_mapq, chrom, no_indel, verbose):
    start_time = time.time()
    
    # clean parameters
    ## Need to add check parameters function ##
    invcf = os.path.abspath(invcf)
    output = os.path.abspath(output)

    # set logging level
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)
    
    # get reference length
    logging.info("Reading reference file")
    ref_len = get_ref_len(ref, chrom)

    # read vcf file
    logging.info("Reading input VCF file")
    with Pool(threads) as p:
        subject_variants = p.starmap(read_vcf, [(invcf, c) for c in chrom])
    
    # find signed edges in bam
    logging.info("Finding connections in BAM file")
    with Pool(threads) as p:
        signed_edges = p.starmap(find_signed_edges, [(, c) for c in chrom])

    # triangle consistent verify


    # output to vcf


    # all done
    logging.info("ALL DONE")
    end_time = time.time()
    logging.info(f"Total processing time is {end_time - start_time} seconds")


if __name__ == "__main__":
    main()


