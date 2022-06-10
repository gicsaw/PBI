#!/usr/bin/env python
import sys
import numpy as np
from pbi.nwa import NWalign, Blosum62, read_fa
sys.setrecursionlimit(10000)

USAGE = """
NWalign.py list_file
"""

def main():

    if len(sys.argv) < 3:
        print('nwalign.py query.fasta template.fasta')
        sys.exit()

    file_q = sys.argv[1]
    file_d = sys.argv[2]

    q, q_header = read_fa(file_q)
    d, d_header = read_fa(file_d)

    nwalign = NWalign(scoring_matrix=Blosum62, g_extend=1.0, g_open2=11.0)
    align, score = nwalign.NW(q, d)
#    print(len(align))
#    for q3,d3,q4,d4,count in align:
#        print(count/len(q),count/len(d),count,len(q),len(d))
#        print(q4)
#        print(d4)

    q3, d3, q4, d4, count = align[0]
    print(q4)
    print(d4)

#    query_aligned_sequence, ini, fin = find_query_aligned_sequence(q4, d4)
#    line_out = '>query:%s|template:%s|%d-%d' % (file_q, file_d, ini, fin)
#    print(line_out)
#    print(query_aligned_sequence)


if __name__ == '__main__':
    main()
