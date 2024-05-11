# <img alt="dna with letters FT" src="https://bananasrlowkeygood.github.io/images/finaletools_logo_rounded.png" height="60"> ‎ ‎ ‎FinaleToolkit
<summary><h3>Table of Contents</h2></summary>
<ol>
  <li><a href="#about-the-project">About The Project</a></li>
  <li><a href="#installation">Installation</a></li>
  <li>
    <a href="#usage">Usage</a>
    <ul>
      <li><a href="#documentation">Documentation</a></li>
      <li><a href="#compatible-file-formats">Compatible File Formats</a></li>
      <li><a href="#retrieving-fragment-files">Generating Fragment Data</a></li>
    </ul>
  </li>
  <li><a href="#testing">Testing</a></li>
  <li><a href="#license">License</a></li>
  <li><a href="#contact">Contact</a></li>
</ol>




## About The Project
FinaleToolkit (**F**ragmentat**I**o**N** **A**na**L**ysis of c**E**ll-free DNA **Toolkit**) is a package and standalone program to extract fragmentation features of cell-free DNA from paired-end sequencing data.

## Installation
You can install the package using `pip`.
```
$ pip install finaletoolkit
```

## Usage

### Documentation
Documentation for FinaleToolkit can be found [here](https://epifluidlab.github.io/finaletoolkit-docs/).

### Compatible File Formats


### Retrieving Fragment Files

## Testing
## Contact
- Ravi Bandaru: ravi.bandaru@northwestern.edu
- Jiayan Ma: jiayan.ma@northwestern.edu
- Yaping Liu: yaping@northwestern.edu




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



To run unit tests, navigate to the root directory of your local copy of this
repo and run `pytest`. You may have to download pytest first.

## FAQ
Q: When running on an ARM64 Mac, I can install FinaleToolkit without errors.
However, I get an `ImportError` when I run it.

A: Try `brew install curl`. Otherwise, email me and I will try to help you.
