
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import src.finaletools.finaletools as finaletools
import numpy as np
import pysam

# scores = np.array([[1, 2], [3, 4]])
# scores = finaletools.wps("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1", 12500, 17500, verbose=True)
# print(scores)

#lengths = finaletools.frag_length("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", contig="2" ,verbose=True)
#with open(f'{os.getcwd()}/tests/test_frag_wps.npy', 'wb') as f_out:
#    np.save(f_out, scores)
#print(lengths[0:50])

# coverage = finaletools.frag_center_coverage("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1", start=10000, stop=30000, verbose=True)
#print(coverage)
"""
for i in range(10):
    foo = finaletools.frag_center_coverage("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1", start=10000, stop=20000, verbose=True)
"""

with pysam.AlignmentFile("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", 'rb') as file:
    # coverage = finaletools.frag_center_coverage(file, "1", start=10000, stop=30000, verbose=True)
    # print(coverage)
    # lengths = finaletools.frag_length(file, contig="2", verbose=True)
    # print(lengths[0:50])
    for region_size in np.arange(1000, 10000, 1000):
        finaletools.wps(file, "1", 8000, 8000+region_size, verbose=True)

