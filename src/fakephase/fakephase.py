import time
import click
import sys
import os
import logging
from multiprocessing import Pool
from fakephase.modules.read_vcf import read_vcf
from fakephase.modules.get_ref_len import get_ref_len
from fakephase.modules.find_signed_edges import find_signed_edges
from fakephase.modules.build_phase_blocks import build_phase_blocks
from fakephase.modules.fake_blocks import fake_blocks
from fakephase.modules.output_vcf import output_vcf



@click.command()
@click.option("-i", "--invcf", required = True, help = "Input vcf/vcf.gz file for haplotype phasing")
@click.option("-b", "--bam", required = True, help = "Input sorted bam file for haplotype phasing")
@click.option("-r", "--ref", required = True, help = "Reference file fasta")
@click.option("-o", "--output", required = True, help = "Output vcf file")
@click.option("-t", "--threads", default = 24, help = "Maximum numbers of parallel threads [default: 24]")
@click.option("--mincoverage", default = 10, help = "Minimal signed edges to phase two variants [default: 10]")
@click.option("--conf", default = 0.8, help = "Minimal confidence to phase two variants [default: 0.9]")
@click.option("--maxlen", default = 10**6, help = "Maximal length to process [default: 10**6]")
@click.option("--min-mapq", default = 0, help = "Minimal mapping quality (MAPQ)")
@click.option("--chrom", default = ','.join([f"chr{i}" for i in range(1, 23)]), help = "Chromosome to evaluate,use comma to join chromosome name e.g. --chrom chr1,chr2,chr3 [default: chr1-22]")
@click.option("--no-indel", is_flag = True, help = "Do not phase insertion and deletion (INDEL)")
@click.option("--no-low-qual", is_flag = True, help = "Do not phase low quality variant (LowQual in filter column)")
@click.option("--verbose", is_flag = True, help = "Verbose mode print progess to standard output")
@click.version_option(version = "0.1.0", prog_name = r"Fake haplotype phasing algorithm deceive most evalation metrics, based on Python 3.10+")
def main(invcf, bam, ref, output, threads, mincoverage, conf, maxlen, min_mapq, chrom, no_indel, no_low_qual, verbose):
    start_time = time.time()
    
    # clean parameters
    ## Need to add check parameters function ##
    invcf = os.path.abspath(invcf)
    output = os.path.abspath(output)    
    chrom = chrom.split(',')
    
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
        subject_variants = p.starmap(read_vcf, [(invcf, c, no_indel, no_low_qual) for c in chrom])
    
    # find signed edges in bam
    logging.info("Finding signed edges from reads")
    with Pool(threads) as p:
        signed_edges = p.starmap(find_signed_edges, [(bam, maxlen, min_mapq, subject_variants[c], ref_len, chrom[c], ref) for c in range(len(chrom))])

    # triangle consistent verify
    logging.info("Verfying signed edges")
    with Pool(threads) as p:
        phased_blocks = p.starmap(build_phase_blocks,[(signed_edges[c], mincoverage, conf) for c in range(len(chrom))])

    # fake phased blocks 
    logging.info("Frabricating phased blocks")
    with Pool(threads) as p:
        fake_blocks = p.starmap(fake_blocks, [(ref_len, chrom[c], subject_variants[c], phased_blocks[c], signed_edges[c], max_len) for c in range(len(chrom))]

    # output to vcf
    logging.info("Writing output VCF")
    write_vcf(invcf, output, fake_blocks, chrom)    

    # all done
    logging.info("ALL DONE")
    end_time = time.time()
    logging.info(f"Total processing time is {end_time - start_time} seconds")



if __name__ == "__main__":
    main()


