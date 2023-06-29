#!/usr/bin/env python3

from __future__ import annotations
import argparse
import importlib
import sys

from finaletools.frag.frag_length import frag_length
from finaletools.frag.coverage import frag_center_coverage
from finaletools.frag.wps import wps
from finaletools.frag.agg_wps import aggregate_wps
from finaletools.frag.delfi import delfi


# TODO: implement subcommands read from stdin
# TODO: implement pipelining
def main_cli():
    parser = argparse.ArgumentParser(
        description='Calculates fragmentation features given a CRAM/BAM/SAM '
        'file',
        epilog='')
    subparsers = parser.add_subparsers(title='subcommands',
                                       dest='subcommand')

    # Common arguments

    # Subcommand 1: frag-coverage
    parser_command1 = subparsers.add_parser('frag-center-coverage',
                                            description=(
                                                'Calculates fragmentation '
                                                'coverage over a region given '
                                                'a CRAM/BAM/SAM file')
                                            )
    parser_command1.add_argument('input-file')
    parser_command1.add_argument('contig')
    # inclusive location of region start in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command1.add_argument('--start', default=0, type=int)
    # exclusive location of region end in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command1.add_argument('--stop', type=int)

    parser_command1.add_argument('-o', '--output_file')
    # parser_command1.add_argument('--method', default="frag-center")
    parser_command1.add_argument('-q', '--quality_threshold', default=30, type=int)
    parser_command1.add_argument('-v', '--verbose', action='store_true', default=0)
    parser_command1.set_defaults(func=frag_center_coverage)

    # Subcommand 2: frag-length
    parser_command2 = subparsers.add_parser(
        'frag-length', prog='frag-length',
        description='Calculates fragment lengths given a CRAM/BAM/SAM file'
        )
    parser_command2.add_argument('input_file')
    parser_command2.add_argument('-c', '--contig')
    parser_command2.add_argument('-o', '--output_file')
    parser_command2.add_argument('-w', '--workers', default=1, type=int)
    parser_command2.add_argument('-q', '--quality_threshold', default=30, type=int)
    parser_command2.add_argument('-v', '--verbose', action='count', default=0)
    parser_command2.set_defaults(func=frag_length)

    # Subcommand 3: wps()
    parser_command3 = subparsers.add_parser(
        'wps', prog='wps',
        description='Calculates Windowed Protection Score over a region given '
        'a CRAM/BAM/SAM file'
        )
    parser_command3.add_argument('input-file')
    parser_command3.add_argument('contig')
    # inclusive location of region start in 0-based coordinate system.
    # If not included, will start at the beginning of the chromosome
    parser_command3.add_argument('start', type=int)
    # exclusive location of region end in 0-based coordinate system.
    # If not included, will end at the end of the chromosome
    parser_command3.add_argument('stop', type=int)
    parser_command3.add_argument('-o', '--output-file')
    parser_command3.add_argument('--window-size', default=120, type=int)
    parser_command3.add_argument('-lo', '--fraction-low', default=120,
                                 type=int)
    parser_command3.add_argument('-hi', '--fraction-high', default=180,
                                 type=int)
    parser_command3.add_argument('--quality-threshold', default=30, type=int)
    parser_command3.add_argument('-v', '--verbose', action='count', default=0)
    parser_command3.set_defaults(func=wps)

    # Subcommand 4: aggregate-wps
    parser_command4 = subparsers.add_parser(
        'aggregate-wps',
        prog='aggregate-wps',
        description='Calculates Windowed Protection Score over a region '
        'around sites specified in a BED file from alignments in a '
        'CRAM/BAM/SAM file'
        )
    parser_command4.add_argument('input-file')
    parser_command4.add_argument('site-bed')
    parser_command4.add_argument('-o', '--output-file')
    parser_command4.add_argument('--size-around-sites', default=5000, type=int)
    parser_command4.add_argument('--window_size', default=120, type=int)
    parser_command4.add_argument('-lo', '--fraction_low', default=120,
                                 type=int)
    parser_command4.add_argument('-hi', '--fraction_high', default=180,
                                 type=int)
    parser_command4.add_argument('--quality_threshold', default=30, type=int)
    parser_command4.add_argument('--workers', default=1, type=int)
    parser_command4.add_argument('-v', '--verbose', action='count', default=0)
    parser_command4.set_defaults(func=aggregate_wps)

    # Subcommand 5: delfi
    parser_command5 = subparsers.add_parser(
        'delfi',
        prog='delfi',
        description='Calculates DELFI score over genome'
        )
    parser_command5.add_argument('input_file')
    parser_command5.add_argument('autosomes')
    parser_command5.add_argument('reference-file')
    parser_command5.add_argument('-b', '--blacklist-file')
    parser_command5.add_argument('-o', '--output-file')
    parser_command5.add_argument('--window-size', default=100000, type=int)
    parser_command5.add_argument('--quality_threshold', default=30, type=int)
    parser_command5.add_argument('-w', '--workers', default=1, type=int)
    parser_command5.add_argument('-v', '--verbose', action='count', default=0)
    parser_command5.set_defaults(func=delfi)


    args = parser.parse_args()
    function = args.func
    funcargs = vars(args)
    funcargs.pop('func')
    funcargs.pop('subcommand')
    # print(funcargs)
    function(**funcargs)

if __name__ == '__main__':
    main_cli()