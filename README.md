# FinaleTools

Lightweight Python library and standalone program for analysis of cell free DNA fragment data. FinaleTools refers to FragmentatIoN AnaLysis of cEll-free DNA Tools.

## Installation
This package is currently work-in-progress. Here is how to install (wip):

- Clone repo from https://github.com/epifluidlab/FinaleTools/
- Create a conda environment with all of the dependencies. If only using pip to install dependencies, create venv/conda env and go to next step.
  Pip might not install all non-python dependencies. Currently only samtools needs to be installed in addition.
    - For conda, run the following:
    ```
    conda create -n <env_name> -c bioconda -c conda-forge python=3.<version>
    <dependencies>
    - see pyproject.toml for dependencies
    ```
    - Alternatively if using Linux, you can also use the attached YAML file by navigating to the project directory and running:
      ```
      conda create -f conda_env/environment.yml
      ```
      The resulting environment should be called `samenv3`.
- cd to project directory
- Run `pip install -e .`
- Run `finaletools -h` to see if you did this right

## Usage
Documentation can be found at https://epifluidlab.github.io/finaletools-docs/

FinaleTools functions generally accept reads in two file formats:
- Binary Alignment Map (BAM) Files
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

Because FinaleTools uses pysam, BAM files should be bai-indexed and Frag.gz files should be tabix-indexed.

To verify FinaleTools has been successfully installed, try
```
$ finaletools -h
usage: finaletools [-h]
                   {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
                   ...

Calculates fragmentation features given a CRAM/BAM/SAM file

options:
  -h, --help            show this help message and exit

subcommands:
  {coverage,frag-length,frag-length-bins,frag-length-intervals,wps,delfi,filter-bam,adjust-wps,agg-wps,delfi-gc-correct,end-motifs,mds}
```

To view fragment length distribution
```
$ finaletools frag-length-bins --contig 22 --histogram sample.bam
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
