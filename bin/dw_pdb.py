#!/usr/bin/env python
import sys
import os
from urllib import request
import gzip
import shutil

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
        line_pdb = 'https://files.rcsb.org/download/%s.pdb.gz' % pdb_id
        pdb_gz_file = '%s/%s.pdb.gz' % (pdb_dir, pdb_id)
        pdb_file = '%s/%s.pdb' % (pdb_dir, pdb_id)

        request.urlretrieve(line_pdb, pdb_gz_file)

        with gzip.open(pdb_gz_file, 'rb') as fp_gz:
            with open(pdb_file, 'wb') as fp_pdb:
                shutil.copyfileobj(fp_gz, fp_pdb)

        line_fasta = 'https://www.rcsb.org/fasta/entry/%s' % pdb_id
        fasta_file = '%s/%s.fasta' % (pdb_dir, pdb_id)
        request.urlretrieve(line_fasta, fasta_file)

#        with request.urlopen(url) as r:
#            lines = r.read().decode('utf-8')
#            with open(fasta_file, 'w') as fp:
#                fp.write(lines)


if __name__ == "__main__":
    main()
