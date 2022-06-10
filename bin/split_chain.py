#!/usr/bin/env python
from pbi.pdbtools import PDBtools


def main():

    import argparse
    title_line = 'split chain from pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-d', '--out_dir', type=str, required=False,
                        default='.', help='output directory')
    parser.add_argument('-e', '--exclude_water', action='store_true',
                        required=False,
                        default=False, help='if exist, exclude water atoms')

    args = parser.parse_args()
    input_file = args.input_file
    out_dir = args.out_dir
    exclude_water = args.exclude_water
    pdb_code = input_file.split('/')[-1].split('.')[0]
    PDBtools.split_chain(input_file, out_dir, pdb_code,
                         exclude_water=exclude_water)


if __name__ == "__main__":
    main()
