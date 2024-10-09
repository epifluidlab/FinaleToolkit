from importlib.resources import files
from pathlib import Path

import pandas as pd

import finaletoolkit.frag as pkg_data

# path to csv containing 5mb bins from delfi_scripts
BINS_PATH: Path = (
    files(pkg_data) /'data' / 'Jan28.hg19.5mb_features_formatted.csv'
)

def delfi_merge_bins(
        hundred_kb_bins: pd.DataFrame,
        gc_corrected: bool=True,
        verbose: bool=False
) -> pd.DataFrame:

    # TODO: allow for no gc-correction
    five_mb_bins = []

    for arm in hundred_kb_bins["arm"].unique():
        arm_bins = hundred_kb_bins[hundred_kb_bins["arm"] == arm].reset_index()
        if "p" in arm:
            for i in range(0,arm_bins.shape[0],50):
                chunk = arm_bins.loc[i:i+49,:]
                if chunk.shape[0] < 50:
                    continue
                contig = arm[:-1]
                start = chunk['start'].min()
                stop = chunk['stop'].max()
                short = chunk['short'].sum()
                long = chunk['long'].sum()
                gc = chunk['gc'].mean()
                num_frags = chunk['num_frags'].sum()
                ratio = chunk['ratio'].mean()
                short_corrected = chunk['short_corrected'].sum()
                long_corrected = chunk['long_corrected'].sum()
                num_frags_corrected = chunk['num_frags_corrected'].sum()
                ratio_corrected = chunk['ratio_corrected'].mean()
                five_mb_bins.append((
                    contig,
                    start,
                    stop,
                    arm,
                    short,
                    long,
                    gc,
                    num_frags,
                    ratio,
                    short_corrected,
                    long_corrected,
                    num_frags_corrected,
                    ratio_corrected
                ))
        elif "q" in arm:
            five_mb_bins_reversed = []
            for i in range(arm_bins.shape[0]-1,0,-50):
                chunk = arm_bins.loc[i-49:i,:]
                if chunk.shape[0] < 50:
                    continue
                contig = arm[:-1]
                start = chunk['start'].min()
                stop = chunk['stop'].max()
                short = chunk['short'].sum()
                long = chunk['long'].sum()
                gc = chunk['gc'].mean()
                num_frags = chunk['num_frags'].sum()
                ratio = chunk['ratio'].mean()
                short_corrected = chunk['short_corrected'].sum()
                long_corrected = chunk['long_corrected'].sum()
                num_frags_corrected = chunk['num_frags_corrected'].sum()
                ratio_corrected = chunk['ratio_corrected'].mean()
                five_mb_bins_reversed.append((
                    contig,
                    start,
                    stop,
                    arm,
                    short,
                    long,
                    gc,
                    num_frags,
                    ratio,
                    short_corrected,
                    long_corrected,
                    num_frags_corrected,
                    ratio_corrected
                ))
            for row in five_mb_bins_reversed[-1::-1]:
                five_mb_bins.append(row)
        del arm_bins
    
    five_mb_bins_df = pd.DataFrame(
        five_mb_bins,
        columns = hundred_kb_bins.columns[hundred_kb_bins.columns!='index'])

    return five_mb_bins_df