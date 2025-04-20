# fakephase
## Critical warning
:warning: <mark>**Do not use this software for haplotype phasing in any general scientific research, clinical applications, industrial production, or other medical, commercial, and political applications.**</mark>

## Introduction
Fake telomere-to-telomere haloptype phasing algorithm. 

## Installation
### Option1. local installation
1. Create new virtual environment
```
conda create -n fakephase python=3.7;
conda activate fakephase;
```

2. Download fakephase and install with `pip`
```
git clone https://github.com/JMencius/fakephase;
cd fakephase;
pip install .;
```


## Usage
```
Usage: fakephase [OPTIONS]

Options:
  -i, --invcf TEXT       Input vcf/vcf.gz file for haplotype phasing
                         [required]
  -b, --bam TEXT         Input sorted bam file for haplotype phasing
                         [required]
  -r, --ref TEXT         Reference file fasta  [required]
  -o, --output TEXT      Output vcf file  [required]
  -t, --threads INTEGER  Maximum numbers of parallel threads [default: 24]
  --mincoverage INTEGER  Minimal signed edges to phase two variants [default:
                         30]
  --conf FLOAT           Minimal confidence to phase two variants [default: 0.9]
  --maxratio FLOAT       Maximal ratio of haplotype to process [default: 0.05]
  --min-mapq INTEGER     Minimal mapping quality (MAPQ)
  --chrom TEXT           Chromosome to evaluate,use comma to join chromosome
                         name e.g. --chrom chr1,chr2,chr3 [default: chr1-22]
  --chunk-size INTEGER   Chunk size length for parallel processing
  --no-low-qual          Do not phase low quality variant (LowQual in filter
                         column)
  --verbose              Verbose mode print progess to standard output
  --version              Show the version and exit.
  --help                 Show this message and exit.
```

## Examples
```
fakephase --verbose \
-i variants.vcf \
-b ./HG00733.nanopore.R9.hac.sorted.bam \
-r GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
-o fake_output.vcf;
```

## Resource consumption
In our tests with 176G nanopore BAM file, Fakephase completed processing within 3 hour using 48 threads with an [AMD EPYC 7K62 CPU](https://openbenchmarking.org/s/AMD+EPYC+7K62+48-Core).

The actual performance may vary depending on factors such as I/O speed, memory speed, and CPU capabilities.

