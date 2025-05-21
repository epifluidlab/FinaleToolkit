"""
Contains the `sort_bed_file` function, which reads a BED format file and
yields entries sorted by genomic coordinates.
"""
from __future__ import annotations
from typing import Tuple, Generator
try:
    from typing_extensions import Unpack
except ImportError:
    from typing import Unpack
from .typing import Intervals, ChromSizes


def sort_bed_file(
        bed_file_path: Intervals,
        chrom_sizes_path: ChromSizes | None = None
        ) -> Generator[
            Tuple[str, int, int, Unpack[Tuple[str, ...]]], None, None]:
    """
    Reads a BED format file and yields entries sorted by chromosome and start
    coordinates.

    Args:
        bed_file_path (str or pathlike): Path to the BED file
        chrom_sizes_path (str or pathlike, optional): Path to a
        chrom.sizes file that defines chromosome ordering. If None,
        uses default chromosome sorting logic.

    Yields:
        tuple: Parsed BED entry containing all original fields
    """
    # First, establish chromosome ordering if a chrom.sizes file is provided
    chrom_order = {}
    if chrom_sizes_path:
        with open(chrom_sizes_path, 'r') as f:
            # Create a dictionary with chromosome name as key and its order as
            # value
            for i, line in enumerate(f):
                fields = line.strip().split()
                if fields:
                    chrom_name = fields[0]
                    chrom_order[chrom_name] = i
    
    # Read all entries from the BED file
    entries = []
    with open(bed_file_path, 'r') as f:
        for line in f:
            # Skip comment or header lines
            if (line.startswith('#') or line.startswith('track')
                    or line.startswith('browser')):
                continue
                
            # Parse the line
            fields = line.strip().split('\t')
            if len(fields) >= 3:  # BED format requires at least 3 fields
                chrom = fields[0]
                # Convert start and end to integers for proper sorting
                try:
                    start = int(fields[1])
                    end = int(fields[2])

                    # Store all fields as a tuple to preserve the entire entry
                    entries.append((chrom, start, end, tuple(fields)))
                except ValueError:
                    # Skip lines with non-integer start/end coordinates
                    continue
            else:
                raise ValueError(
                    "Invalid file format for interval file."
                )

    # Custom sort key function that uses the chrom.sizes order if available
    def sort_key(entry):
        chrom, start, end, _ = entry

        # If we have a chrom_order dictionary, use it for chromosome ordering
        if chrom_order:
            # If the chromosome is in our order dict, use its index
            # If not, assign a very high value so it comes after all known
            # chromosomes
            chrom_val = chrom_order.get(chrom, float('inf'))
        else:
            # Default sorting logic for chromosomes when no chrom.sizes file
            # is provided
            if chrom.startswith('chr'):
                # Handle chrX, chrY, chrM specially
                if chrom[3:] == 'X':
                    chrom_val = 23
                elif chrom[3:] == 'Y':
                    chrom_val = 24
                elif chrom[3:] in ('M', 'MT'):
                    chrom_val = 25
                else:
                    # Try to convert to integer if possible
                    try:
                        chrom_val = int(chrom[3:])
                    except ValueError:
                        # For non-standard chromosomes, use the original name
                        chrom_val = chrom
            else:
                # Handle non-prefixed chromosome names
                if chrom == 'X':
                    chrom_val = 23
                elif chrom == 'Y':
                    chrom_val = 24
                elif chrom in ('MT', 'M'):
                    chrom_val = 25
                else:
                    try:
                        chrom_val = int(chrom)
                    except ValueError:
                        chrom_val = chrom

        # Return a tuple of (chrom_val, start, end) for sorting
        # This will use end as a tie-breaker when chrom and start are identical
        return (chrom_val, start, end)

    # Sort the entries using our custom sort key
    entries.sort(key=sort_key)

    # Yield the sorted entries, keeping all original fields
    for _, _, _, fields in entries:
        # cast output start and end to ints
        fields[1] = int(fields[1])
        fields[2] = int(fields[2])
        yield fields
