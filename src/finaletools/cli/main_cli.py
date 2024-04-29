import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], max_content_width=200)

@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)

def main_cli():
    pass

@main_cli.command(name='coverage-intervals', short_help="Calculates fragmentation coverage over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.option('-o', '--output_file', default='-', type=str, help='Output BED file path.')
@click.option('-s', '--scale_factor', default=1, show_default=True, type=float, help='Scale factor for coverage values.')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in coverage.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in coverage.')
@click.option('-i', '--intersect_policy', default='midpoint', type=str, show_default=True, help='Specifies what policy is used to include fragments in the given interval.')
@click.option('-q', '--quality', default=30, show_default=True, help='Minimum mapping quality threshold.')
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-n', '--normalize', is_flag=True, help='Normalize (divide) coverage values by the total number of fragments.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')

def coverage_cmd(input_file, interval_file, output_file, scale_factor, min_length, max_length, intersect_policy, quality, workers, verbose, normalize):
    """
    Calculates fragmentation coverage over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.

        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to calculate coverage over.

    \033[91mOUTPUT:\033[0m A BED file containing coverage values over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.coverage import coverage_intervals
    coverage_intervals(input_file, interval_file, output_file=output_file, scale_factor=scale_factor, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, quality_threshold=quality, workers=workers, normalize=normalize, verbose=verbose)

@main_cli.command(name='frag-length-bins', short_help="Retrieves fragment lengths given a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-c', '--contig', type=str, help='Specify the contig or chromosome to select fragments from. (Required if using --start or --stop.)')
@click.option('-s', '--start', type=int, help='Specify the 0-based left-most coordinate of the interval to select fragments from. (Must also specify --contig.)')
@click.option('-e', '--stop', type=int, help='Specify the 1-based right-most coordinate of the interval to select fragments from. (Must also specify --contig.)')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-i', '--intersect_policy', default='midpoint', type=str, show_default=True, help='Specifies what policy is used to include fragments in the given interval.')
@click.option('-b', '--bin_size', default=1, show_default=True, type=int, help='Specify the size of the bins to group fragment lengths into.')
@click.option('-o', '--output_file', default='-', type=str, help='Output .tsv file path.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-H', '--hist_path', default=None, type=str, help='Specifies the directory to which to save the .png histogram plot.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def frag_length_bin(input_file, contig, start, stop, min_length, max_length, intersect_policy, bin_size, output_file, quality, hist_path, verbose):
    """
    Retrieves fragment lengths given a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.

    \033[91mOUTPUT:\033[0m A .TSV file containing containing fragment lengths binned according to the specified bin size.
    """
    from finaletools.frag.frag_length import frag_length_bins
    frag_length_bins(input_file, contig=contig, start=start, stop=stop, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, bin_size=bin_size, output_file=output_file, quality_threshold=quality, histogram_path=hist_path, verbose=verbose)

@main_cli.command(name='frag-length-intervals', short_help="Retrieves fragment length summary statistics over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.option('-o', '--output_file', default='-', type=str, help='Output BED file path.')
@click.option('-i', '--intersect_policy', default='midpoint', type=str, show_default=True, help='Specifies what policy is used to include fragments in the given interval.')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def frag_length_intervals(input_file, interval_file, output_file, min_length, max_length, quality, intersect_policy, workers, verbose):
    """
    Retrieves fragment length summary statistics over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to retrieve fragment length summary statistics over.

    \033[91mOUTPUT:\033[0m A BED file containing fragment length summary statistics (mean, median, st. dev, min, max) over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.frag_length import frag_length_intervals
    frag_length_intervals(input_file, interval_file, output_file=output_file, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, quality_threshold=quality, workers=workers, verbose=verbose)


@main_cli.command(name='cleavage-profile', short_help="Calculates cleavage proportion over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.argument('output_file', type=str)
@click.option('--five/--three', default=True, help="Calculate cleavage proportion based on 5' or 3' end of fragments.")
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-c', '--chrom_sizes', default=None, type=click.Path(exists=True), help='Path to a chrom.sizes file. This is required and only used if a fragment file (.bed.frag.gz) is provided as input.')
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def cleavage_profile(input_file, interval_file, output_file, min_length, max_length, quality, chrom_sizes, five, workers, verbose):
    """
    Retrieves fragment length summary statistics over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to calculates cleavage proportion over.

    \033[91mOUTPUT:\033[0m A bigWig file containing the cleavage proportion results over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.cleavage_profile import cleavage_profile_intervals
    cleavage_profile_intervals(input_file, interval_file, output_file=output_file, min_length=min_length, max_length=max_length, quality_threshold=quality, workers=workers, verbose=verbose, chrom_sizes=chrom_sizes, five=five)

@main_cli.command(name='to-frag', short_help="Converts a BAM file to a Fragment file.", no_args_is_help=True)
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_file', type=str)
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
def to_frag(input_bam, output_file, workers):
    """
    Converts a BAM file to a Fragment file.

    Arguments:

        \033[92mINPUT_BAM:\033[0m Path to a BAM file containing fragment data.


        \033[92mOUTPUT_FILE:\033[0m Path to the output .bed.frag.gz file.

    \033[91mOUTPUT:\033[0m A Fragment file (.bed.frag.gz) containing fragment data.
    """
    from finaletools.utils.to_frag import toFrag
    toFrag(input_bam, output_file, workers=workers)

@main_cli.command(name='wps', short_help="Calculates Windowed Protection Score over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.argument('output_file', type=str)
@click.option('-W', '--window_size', default=120, show_default=True, type=int, help='Size of the sliding window used to calculate WPS scores.')
@click.option('-i', '--intersect_policy', default='midpoint', type=str, show_default=True, help='Specifies what policy is used to include fragments in the given interval.')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-c', '--chrom_sizes', default=None, type=click.Path(exists=True), help='Path to a chrom.sizes file. This is required and only used if a fragment file (.bed.frag.gz) is provided as input.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def wps(input_file, interval_file, output_file, window_size, min_length, max_length, quality, intersect_policy, chrom_sizes, workers, verbose):
    """
    Calculates Windowed Protection Score (WPS) over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to calculate WPS over.

    \033[91mOUTPUT:\033[0m A bigWig file containing the WPS results over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.multi_wps import multi_wps
    multi_wps(input_file, interval_file, output_file, window_size=window_size, chrom_sizes=chrom_sizes, min_length=min_length, max_length=max_length, intersect_policy=intersect_policy, quality_threshold=quality, workers=workers, verbose=verbose)

@main_cli.command(name='wps', short_help="Calculates Windowed Protection Score (WPS) over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.argument('interval_size', type=int)
@click.argument('output_file', type=str)
@click.option('-W', '--window_size', default=120, show_default=True, type=int, help='Size of the sliding window used to calculate WPS scores.')
@click.option('-i', '--intersect_policy', default='midpoint', type=str, show_default=True, help='Specifies what policy is used to include fragments in the given interval.')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-c', '--chrom_sizes', default=None, type=click.Path(exists=True), help='Path to a chrom.sizes file. This is required and only used if a fragment file (.bed.frag.gz) is provided as input.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def wps(input_file, interval_file, output_file, window_size, interval_size, min_length, max_length, quality, intersect_policy, chrom_sizes, workers, verbose):
    """
    Calculates Windowed Protection Score (WPS) over intervals defined in a BED file based on alignment data from a BAM/SAM/CRAM/Fragment file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to a BAM/SAM/CRAM/Fragment file containing fragment data.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to calculate WPS over.


        \033[92mINTERVAL_SIZE:\033[0m Constant size of each interval specified in INTERVAL_FILE. 

    \033[91mOUTPUT:\033[0m A bigWig file containing the WPS results over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.multi_wps import multi_wps
    multi_wps(input_file, interval_file, output_file, window_size=window_size, chrom_sizes=chrom_sizes, min_length=min_length, max_length=max_length, interval_size=interval_size, intersect_policy=intersect_policy, quality_threshold=quality, workers=workers, verbose=verbose)

@main_cli.command(name='adjust-wps', short_help="Adjusts raw Windowed Protection Score (WPS) by applying a median filter and Savitsky-Golay filter.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.argument('output_file', type=str)
@click.argument('chrom_sizes', type=click.Path(exists=True))
@click.option('-m', '--median_window_size', default=1000, show_default=True, type=int, help='Size of the median filter window used to adjust WPS scores.')
@click.option('-s', '--savgol_window_size', default=21, type=int, show_default=True, help='Size of the Savitsky-Golay filter window used to adjust WPS scores.')
@click.option('-p', '--savgol_poly_deg', default=2, type=int, show_default=True, help='Degree polynomial for Savitsky-Golay filter.')
@click.option('-M', '--mean', is_flag=True, help='A mean filter is used instead of median.')
@click.option('-e', '--subtract_edges', is_flag=True, help='Take the median of the first and last --edge_size bases in a window and subtract from the whole interval.')
@click.option('-E', '--edge_size', default=500, type=int, show_default=False, help='Size of the edge to subtract from the whole interval.')
@click.option('-w', '--workers', default=1, type=int, show_default=True, help='Number of worker processes.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def adjust_wps(input_file, interval_file, output_file, chrom_sizes, median_window_size, savgol_window_size, savgol_poly_deg, mean, subtract_edges, edge_size, workers, verbose):
    """
    Adjusts raw Windowed Protection Score (WPS) by applying a median filter and Savitsky-Golay filter.

    Arguments:

        \033[92mINPUT_FILE:\033[0m A bigWig file containing the WPS results over the intervals specified in INTERVAL_FILE.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to WPS was calculated over.


        \033[92mCHROM_SIZES:\033[0m A .chrom.sizes file containing chromosome sizes.

    \033[91mOUTPUT:\033[0m A bigWig file containing the adjusted WPS results over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.adjust_wps import adjust_wps
    adjust_wps(input_file=input_file, interval_file=interval_file, output_file=output_file, chrom_sizes=chrom_sizes, median_window_size=median_window_size, savgol_window_size=savgol_window_size, savgol_poly_deg=savgol_poly_deg, mean=mean, subtract_edges=subtract_edges, edge_size=edge_size, workers=workers, verbose=verbose)


@main_cli.command(name='agg-signal', short_help="Aggregates a bigWig signal over constant-length intervals defined in a BED file.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('interval_file', type=click.Path(exists=True))
@click.argument('output_file', type=str)
@click.option('-m', '--median_window_size', default=0, show_default=True, type=int, help='Size of the median filter window used to adjust WPS scores. Only modify if aggregating WPS signals.')
@click.option('-M', '--mean', is_flag=True, help='A mean aggregation is returned rather than sum.')
@click.option('-s', '--strand_location', default=5, show_default=True, help='Index of column in --interval_file that contains the strand (+/-) information.')
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def agg_signal(input_file, interval_file, output_file, median_window_size, mean, strand_location, verbose):
    """
    Aggregates a bigWig signal over constant-length intervals defined in a BED file.

    Arguments:

        \033[92mINPUT_FILE:\033[0m A bigWig file containing signals over the intervals specified in INTERVAL_FILE.


        \033[92mINTERVAL_FILE:\033[0m Path to a BED file containing intervals to which signals were calculated over.

    \033[91mOUTPUT:\033[0m A wiggle file containing the aggregate signal over the intervals specified in INTERVAL_FILE.
    """
    from finaletools.frag.agg_bw_signal import agg_bw_signal
    agg_bw_signal(input_file=input_file, interval_file=interval_file, output_file=output_file, median_window_size=median_window_size, mean=mean, strand_location=strand_location, verbose=verbose)

@main_cli.command(name='filter-bam', short_help="Filters a BAM file for mapped read pairs meeting specified criteria.", no_args_is_help=True)
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-r', '--region_file', default=None, type=click.Path(exists=True), help='Only output alignments overlapping the intervals in this BED file will be included.')
@click.option('-o', '--output_file', default="-", type=str, help='Output BAM file path.')
@click.option('-min', '--min_length', default=None, type=int, help='Minimum length for a fragment to be included in fragment length.')
@click.option('-max', '--max_length', default=None, type=int, help='Maximum length for a fragment to be included in fragment length.')
@click.option('-q', '--quality', show_default=True, default=30, type=int, help="Minimum mapping quality threshold.")
@click.option('-v', '--verbose', is_flag=True, help='Enable verbose mode to display detailed processing information.')
def filter_bam(input_file, interval_file, output_file, median_window_size, mean, strand_location, verbose):
    """
    Filters a BAM file for mapped read pairs meeting specified criteria.

    Arguments:

        \033[92mINPUT_FILE:\033[0m Path to BAM file.


    \033[91mOUTPUT:\033[0m A BAM file containing filtered read pairs.
    """
    from finaletools.utils.filter_bam import filter_bam
    filter_bam(input_file=input_file, region_file=region_file, output_file=output_file, min_length=min_length, max_length=max_length, quality=quality, verbose=verbose)

if __name__ == '__main__':
    main_cli()
