# FinaleToolkit
A package and standalone program to extract fragmentation patterns of cell-free
DNA from paired-end sequencing data. FinaleToolkit refers to FragmentatIoN
AnaLysis of cEll-free DNA Tools.

FinaleToolkit is in active development, and all API is subject to change and
should be considered unstable.

## Installation
Instructions:
- (Optional) create a conda or venv environment to use FinaleToolkit in.
- Run `pip install finaletoolkit`

To verify FinaleToolkit has been successfully installed, try
```
$ finaletoolkit -h
usage: finaletoolkit [-h]
                   {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
                   ...

Calculates fragmentation features given a CRAM/BAM/SAM file

options:
  -h, --help            show this help message and exit

subcommands:
  {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
```

## Usage
Documentation can be found at https://epifluidlab.github.io/finaletoolkit-docs/

FinaleToolkit functions generally accept reads in a few file formats:
- Binary Alignment Map (BAM) Files
- Compressed Reference-oriented Alignment Map
- FinaleDB Frag.gz Files

Frag.gz files are block-gzipped BED3+2 files with the following format:
`chrom  start  stop  mapq  strand(+/-)`

The below script can be used to convert from bam to frag.gz:
```
INPUT=input.bam
OUTPUT=output.frag.gz

samtools sort -n -o qsorted.bam -@ 16 input.bam;
samtools view -h -f 3 -F 3852 -G 48 --incl-flags 48 \
  qsorted.bam |\
  bamToBed -bedpe -mate1 -i stdin |\
  awk -F'\t' -v OFS="\t" '{if ($1!=$4) next; if ($9=="+") {s=$2;e=$6} else {s=$5;e=$3} if (e>s) print $1,s,e,$8,$9}' |\
  sort -k1,1V -k2,2n |\
  bgzip > $OUTPUT;
tabix -p bed $OUTPUT;
```

Frag.gz files can be retrieved from http://finaledb.research.cchmc.org/

Because FinaleToolkit uses pysam, BAM files should be bai-indexed and Frag.gz files should be tabix-indexed.

To view fragment length distribution
```
$ finaletoolkit frag-length-bins --contig 22 --histogram sample.bam
Fragment Lengths for 22:-
10.61%                            ▇                              mean      :169.28
09.85%                           ▆█▁                             median    :169.00
09.09%                           ███                             stdev     :25.52
08.34%                           ████                            min       :67.00
07.58%                          ▁████                            max       :289.00
06.82%                          █████▂                          
06.06%                          ██████                          
05.31%                         ▆██████▂                         
04.55%                        ▄████████▁                        
03.79%                       ▃██████████                        
03.03%                     ▂████████████▆                       
02.27%                     ██████████████▇▃                     
01.52%                   ▇█████████████████▅▂                   
00.76%     ▂▂▂▂▂▂▃▃▄▅▄████████████████████████▆▅▄▃▂▂▂▂▂▂▂▁▁▂▂▂▁▁
len (nt)067   091   115   139   163   187   211   235   259   283
```

## Testing

To run unit tests, navigate to the root directory of your local copy of this
repo and run `pytest`. You may have to download pytest first.

## FAQ
Q: When running on an ARM64 Mac, I can install FinaleToolkit without errors.
However, I get an `ImportError` when I run it.

A: Try `brew install curl`. Otherwise, email me and I will try to help you.