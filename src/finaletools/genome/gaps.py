from __future__ import annotations
from typing import Union
import gzip
try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import finaletools.genome as genome

HG19GAPS: Path = (files(genome) / 'data' / 'hg19.gap.txt.gz')
HG38GAPS: Path = (files(genome) / 'data' / 'hg19.gap.txt.gz')
HG38CENTROMERES: Path = (files(genome) / 'data' / 'hg38.centromeres.txt.gz')


class GenomeGaps:
    """
    Reads telomere, centromere, and short_arm intervals from a bed file
    or generates these intervals from UCSC gap and centromere tracks for
    hg19 and hg38.
    """
    def __init__(self, gaps_bed: Union[Path, str]=None):
        self.centromeres: NDArray
        self.telomeres: NDArray
        self.short_arms: NDArray
        if gaps_bed is None:
            pass
        else:
            gaps = np.genfromtxt(
                gaps_bed,
                dtype=[
                    ('contig', '<U32'),
                    ('start', '<i8'),
                    ('stop', '<i8'),
                    ('type', '<U32'),
                ]
            )
            self.centromeres = gaps[gaps['type'] == 'centromere']
            self.telomeres = gaps[gaps['type'] == 'telomere']
            self.short_arms = gaps[gaps['type'] == 'short_arm']

    @classmethod
    def ucsc_hg19(cls):
        """
        Creates a GenomeGaps for the UCSC hg19 reference genome. This
        sequences uses chromosome names that start with 'chr' and is
        based on a version of the GRCh37 reference genome.

        Returns
        -------
        gaps : GenomeGaps
            GenomeGaps for the UCSC hg19 reference genome.
        """
        genome_gaps = cls()
        gaps = np.genfromtxt(
            HG19GAPS,
            usecols=[1, 2, 3, 7],
            dtype=[
                ('contig', '<U32'),
                ('start', '<i8'),
                ('stop', '<i8'),
                ('type', '<U32'),
            ]
        )
        genome_gaps.centromeres = gaps[gaps['type'] == 'centromere']
        genome_gaps.telomeres = gaps[gaps['type'] == 'telomere']
        genome_gaps.short_arms = gaps[gaps['type'] == 'short_arm']

        return genome_gaps

    @classmethod
    def b37(cls):
        """
        Creates a GenomeGaps for the Broad Institute GRCh37 reference
        genome i.e b37. This reference genome is also based on GRCh37,
        but differs from the UCSC hg19 reference in a few ways,
        including the absence of the 'chr' prefix. We generate this
        GenomeGap using an ad hoc method where we take the UCSC hg19
        gap track and drop 'chr' from the chromosome names. Because
        there are other differences between hg19 and b37, this is not
        a perfect solution.

        Returns
        -------
        gaps : GenomeGaps
            GenomeGaps for the b37 reference genome.
        """
        genome_gaps = cls()
        gaps = np.genfromtxt(
            HG19GAPS,
            usecols=[1, 2, 3, 7],
            dtype=[
                ('contig', '<U32'),
                ('start', '<i8'),
                ('stop', '<i8'),
                ('type', '<U32'),
            ]
        )
        gaps['contig'] = np.char.replace(gaps['contig'], 'chr', '')
        genome_gaps.centromeres = gaps[gaps['type'] == 'centromere']
        genome_gaps.telomeres = gaps[gaps['type'] == 'telomere']
        genome_gaps.short_arms = gaps[gaps['type'] == 'short_arm']

        return genome_gaps

    @classmethod
    def ucsc_hg38(cls):
        return NotImplemented

    def in_tcmere(self, contig: str, start: int, stop: int) -> bool:
        """
        Checks if specified interval is in a centromere or telomere

        Parameters
        ----------
        contig : str
            Chromosome name
        start : int
            Start of interval
        stop : int
            End of interval

        Returns
        -------
        in_telomere_or_centromere : bool
            True if in a centromere or telomere
        """
        # get centromere and telomeres for contig
        centromere = self.centromeres[self.centromeres['contig'] == contig]
        telomeres = self.telomeres[self.telomeres['contig'] == contig]
        if not (centromere.shape[0] or telomeres.shape[0]):
            return None
        else:
            in_centromere = np.logical_and(
                start >= centromere['start'],
                stop <= centromere['stop'],
            )
            in_telomeres = np.sum(np.logical_and(
                start >= telomeres['start'],
                stop <= telomeres['stop'],
            )) > 0
            return in_centromere or in_telomeres

    def get_arm(self, contig: str, start: int, stop: int) -> str:
        """
        Returns the chromosome arm the interval is in. If in
        the short arm of an acrocentric chromosome or intersects a
        centromere, returns an empty string.

        contig : str
            Chromosome of interval.
        start : int
            Start of interval.
        stop : int
            End of interval.

        Returns
        -------
        arm : str
            Arm that interval is in.

        Raises
        ------
        ValueError
            Raised for invalid coordinates
        """
        if stop < start:
            raise ValueError('start must be less than stop')

        # get centromere and short_arm for contig
        centromere = self.centromeres[self.centromeres['contig'] == contig]
        short_arm = self.short_arms[self.short_arms['contig'] == contig]
        has_short_arm = short_arm.shape[0] > 0
        if stop < centromere['start'][0]:
            if not has_short_arm:
                return f"p{contig.replace('chr', '')}"
            else:
                return ''
        elif start > centromere['stop'][0]:
            return f"q{contig.replace('chr', '')}"
        else:
            return ''


    def to_bed(self, output_file: str):
        """
        Prints gap intervals in GenomeGaps to a BED4 file where the name
        is the type of gap interval.

        Parameters
        ----------
        output_file : str
            File to write to. Optionally gzipped.
        """
        gaps = np.sort(np.append(np.append(self.centromeres, self.telomeres), self.short_arms))
        if output_file.endswith('.gz'):
            with gzip.open(output_file, 'w') as output:
                for interval in gaps:
                    output.write(
                        f"{interval['contig']}\t{interval['start']}\t"
                        f"{interval['stop']}\t{interval['type']}\n"
                    )
        else:
            with open(output_file, 'w') as output:
                for interval in gaps:
                    output.write(
                        f"{interval['contig']}\t{interval['start']}\t"
                        f"{interval['stop']}\t{interval['type']}\n"
                    )


def ucsc_hg19_gap_bed(output_file: str):
    return GenomeGaps.ucsc_hg19().to_bed(output_file)

def b37_gap_bed(output_file: str):
    return GenomeGaps.b37().to_bed(output_file)

def ucsc_hg38_gap_bed(output_file: str):
    return NotImplemented