
Input Data
=========================================

**FinaleToolkit** is compatible with almost any paired-end sequence data.

SAM
^^^^^^^^^^^

A sequence alignment map (SAM) file is a human-readable file format that stores the results of sequence alignment. It contains information about the alignment of each read to the reference genome, as well as information about the read itself. 

BAM
^^^^^^^^^^^

A binary alignment file (BAM) provides the same information as a SAM file, but in a binary format. This can be useful for saving space on disk, but is not human-readable. 

**FinaleToolkit** requires that BAM files be BAI indexed. Therefore, you should have an associated ``.bam.bai`` file in the same directory of your input data.

CRAM
^^^^^^^^^^^

A compressed read alignment map file is a compressed version of a SAM file. It is a binary file that is smaller than a BAM file, but still contains all of the same information. 

**FinaleToolkit** requires that CRAM files be CRAI indexed. Therefore, you should have an associated ``.cram.crai`` file in the same directory of your input data.

Fragment File
^^^^^^^^^^^^^^^^

A fragment file (`.frag.gz`) file that is derived from information in a BAM file. A fragment file is a block-gzipped BED3+2 file (similar to a tab-separated value file) that contains the following columns (with one row entry for each fragment): ``chrom``, ``start``, ``stop``, ``mapq``, and ``strand``.

Here, ``mapq`` is the mapping quality of the fragment, and ``strand`` is the strand of the fragment. The ``strand`` column can be either ``+`` or ``-``.

**FinaleToolkit** requires that fragment files be Tabix indexed. Therefore, you should have an associated ``.frag.gz.tbi`` file in the same directory of your input data.

For your reference, here is an example fragment file::

        #chrom    start    stop    mapq    strand
        chr1    10000    10050    60    +
        chr1    10050    10100    60    -
        chr1    10100    10150    60    +
        chr1    10150    10200    60    -
        chr1    10200    10250    60    +
        chr1    10250    10300    60    -
        chr1    10300    10350    60    +
        chr1    10350    10400    60    -
        chr1    10400    10450    60    +
        chr1    10450    10500    60    -

We encourage you to use our comprehensive database, **FinaleDB**, to access relevant fragment files. Learn more about **FinaleDB** `here <http://finaledb.research.cchmc.org>`_ .