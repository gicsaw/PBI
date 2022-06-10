#!/usr/bin/env python
import sys
import os
from urllib import request

usage = '''
dw_pdb.py list_file pdb_dir
'''


def main():
    if len(sys.argv) < 2:
        print(usage)
        sys.exit()

    list_file = sys.argv[1]
    pdb_dir = 'pdb'
    if len(sys.argv) >= 3:
        pdb_dir = sys.argv[2]

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)

    for line in lines:
        lis = line.strip().split(';')
        pdb_id = lis[0]
#        method = lis[1]
#        resolution = lis[2]
#        chain_positions = lis[3]

        print(pdb_id)
        line_pdb = 'https://files.rcsb.org/download/%s.pdb' % pdb_id
        pdb_file = '%s/%s.pdb' % (pdb_dir, pdb_id)
        request.urlretrieve(line_pdb, pdb_file)

        line_fasta = 'https://www.rcsb.org/fasta/entry/%s' % pdb_id
        fasta_file = '%s/%s.fasta' % (pdb_dir, pdb_id)
        request.urlretrieve(line_fasta, fasta_file)

#        with request.urlopen(url) as r:
#            lines = r.read().decode('utf-8')
#            with open(fasta_file, 'w') as fp:
#                fp.write(lines)


if __name__ == "__main__":
    main()
