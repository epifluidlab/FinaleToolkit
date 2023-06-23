
import sys
import os

import numpy as np
import pysam

import finaletools

# contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT', 'GL000207.1', 'GL000226.1', 'GL000229.1', 'GL000231.1', 'GL000210.1', 'GL000239.1', 'GL000235.1', 'GL000201.1', 'GL000247.1', 'GL000245.1', 'GL000197.1', 'GL000203.1', 'GL000246.1', 'GL000249.1', 'GL000196.1', 'GL000248.1', 'GL000244.1', 'GL000238.1', 'GL000202.1', 'GL000234.1', 'GL000232.1', 'GL000206.1', 'GL000240.1', 'GL000236.1', 'GL000241.1', 'GL000243.1', 'GL000242.1', 'GL000230.1', 'GL000237.1', 'GL000233.1', 'GL000204.1', 'GL000198.1', 'GL000208.1', 'GL000191.1', 'GL000227.1', 'GL000228.1', 'GL000214.1', 'GL000221.1', 'GL000209.1', 'GL000218.1', 'GL000220.1', 'GL000213.1', 'GL000211.1', 'GL000199.1', 'GL000217.1', 'GL000216.1', 'GL000215.1', 'GL000205.1', 'GL000219.1', 'GL000224.1', 'GL000223.1', 'GL000195.1', 'GL000212.1', 'GL000222.1', 'GL000200.1', 'GL000193.1', 'GL000194.1', 'GL000225.1', 'GL000192.1']
# scores = np.array([[1, 2], [3, 4]])

# print(scores)

#lengths = finaletools.frag_length("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", verbose=True)

#print(lengths[0:50])

# coverage = finaletools.frag_center_coverage("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1", start=10000, stop=30000, verbose=True)
#print(coverage)
"""
for i in range(10):
    foo = finaletools.frag_center_coverage("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1", start=10000, stop=20000, verbose=True)
"""

if __name__ == '__main__':

    #scores = finaletools.wps("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", "1",14000, 16000, output_file="/Users/lij3mv/Desktop/Jan28.hg19.1.14000-16000.wps.wig", workers=2, verbose=True)
    #with open(f'{os.getcwd()}/tests/test_frag_wps.npy', 'wb') as f_out:
    #    np.save(f_out, scores)

    contig_dict = finaletools.contig_site_bams('/users/lij3mv/Desktop/Homo_sapiens.b37.75.TSS.commonChr.TJ_GG_Irizarr.sort.bed', '/users/lij3mv/Desktop/human.hg19.nc.genome')
    for line in contig_dict['21'].readlines()[:10]:
        print(line)

    #with pysam.AlignmentFile("/users/lij3mv/Desktop/Jan28.hg19.mdups.bam", 'rb') as file:
        # coverage = finaletools.frag_center_coverage(file, "1", start=10000, stop=30000, verbose=True)
        # print(coverage)
        # lengths = finaletools.frag_length(file, contig="2", verbose=True)
        # print(lengths[0:50])


        #for region_size in range(1000, 6000, 1000):
        #    finaletools.wps(file, "2", 8000, 8000+region_size, fraction_low=120, fraction_high=180,  workers=2, verbose=2)
            # finaletools.frag_center_coverage(file, "1", 10000, 10000+region_size, verbose=True)
        #for pool_size in [1, 2, 4, 8]:
        #    finaletools.wps(file, "1", 14000, 16000, workers=pool_size, verbose=True)