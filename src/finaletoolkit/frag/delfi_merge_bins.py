try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files
from pathlib import Path

import numpy as np
import pandas as pd

import finaletoolkit.frag as pkg_data

# path to csv containing 5mb bins from delfi_scripts
BINS_PATH: Path = (
    files(pkg_data) /'data' / 'Jan28.hg19.5mb_features_formatted.csv'
)

def delfi_merge_bins(
        hundred_kb_bins: pd.DataFrame,
        gc_corrected: bool=True,
        add_chr: bool=False,
        verbose: bool=False
):
    #FIXME: find a way to calculate this that is not dependent on b37 format
    delfi_scripts_5mb_bins = pd.read_csv(BINS_PATH)

    if add_chr:
        chr_bins = hundred_kb_bins.copy()
        chr_bins['contig'] = hundred_kb_bins['contig'].apply(
            lambda x: f"chr{x}")
    else:
        chr_bins = hundred_kb_bins.copy()

    
    five_mb_bins = delfi_scripts_5mb_bins.copy()

    #finding overlaps
    contigs_1 = chr_bins['contig'].to_numpy()[:, np.newaxis]
    starts_1 = chr_bins['start'].to_numpy()[:, np.newaxis]
    stops_1 = chr_bins['stop'].to_numpy()[:, np.newaxis]

    contigs_2 = delfi_scripts_5mb_bins['seqnames'].to_numpy()[np.newaxis]
    starts_2 = delfi_scripts_5mb_bins['start'].to_numpy()[np.newaxis]
    stops_2 = delfi_scripts_5mb_bins['end'].to_numpy()[np.newaxis]

    contig_blind_overlaps = np.logical_and(
        (starts_1 < stops_2),
        (stops_1 > starts_2)
    )
    in_same_contig = contigs_1 == contigs_2
    overlaps = np.logical_and(contig_blind_overlaps, in_same_contig)

    # merging bins
    ft_5mb_indices = np.arange(five_mb_bins.shape[0])
    if gc_corrected:
        for i in ft_5mb_indices:
            subset = (chr_bins.loc[overlaps[:,i],:])
            assert subset.shape[0] == 50, "5mb bin does not have 50 100bins."
            five_mb_bins.iloc[[i],[4]] = subset['short_corrected'].sum()
            five_mb_bins.iloc[[i],[5]] = subset['long_corrected'].sum()
            five_mb_bins.iloc[[i],[6]] = subset['ratio_corrected'].mean()
            five_mb_bins.iloc[[i],[7]] = subset['num_frags_corrected'].sum()
    else:
        for i in ft_5mb_indices:
            subset = (chr_bins.loc[overlaps[:,i],:])
            assert subset.shape[0] == 50, "5mb bin does not have 50 100bins."
            five_mb_bins.iloc[[i],[4]] = subset['short'].sum()
            five_mb_bins.iloc[[i],[5]] = subset['long'].sum()
            five_mb_bins.iloc[[i],[6]] = subset['ratio'].mean()
            five_mb_bins.iloc[[i],[7]] = subset['num_frags'].sum()
    
    return five_mb_bins

# TODO: make binning customizable