
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import finaletools
import numpy as np

lengths = np.array(finaletools.frag_length("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", quality_threshold=))
with open(f'{os.getcwd()}/tests/test_frag_lengths.bin', 'wb') as f_out:
    lengths.tofile(f_out)

