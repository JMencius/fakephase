import click
import sys
import os
import logging
from multiprocessing import Pool


@click.command()
@click.option("-i", "--input", required = True, help = "Input vcf/vcf.gz file for haplotype phasing")
@click.option("-b", "--bam", required = True, help = "Input sorted bam file for haplotype phasing")
@click.option("-r", "--ref", required = True, help = "Reference file fasta")
@click.option("-o", "--output", required = True, help = "Output vcf file")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads")
@click.option("--verbose", is_flag = True, help = "Verbose mode print progess to standard output")
@click.version_option(version = "0.1.0", prog_name = r"Fake haplotype phasing algorithm deceive most evalation metrics, based on Python 3.10+")
def main(input, bam, ref, output, threads, verbose):
    start_time = time.time()
    
    # set logging
    logging.basicConfig(level = logging.DEBUG, format = "%(asctime)s - %(levelname)s - %(message)s")
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARNING)

    



if __name__ == "__main__":
    main()
