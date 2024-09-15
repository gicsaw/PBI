#!/usr/bin/env python
from pbi.ligandtools import pdb_to_mol2



def main():

    import argparse
    title_line = 'pdb2mol2'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input ligand pdb file')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.mol2',  help='output ligand mol2 file')
    parser.add_argument('-k', '--kekule', action='store_true',
                        required=False, help='write molecule to kekule form')

    args = parser.parse_args()
    pdb_file = args.input_file
    mol2_file = args.output_file
    is_kekule = args.kekule


    pdb_to_mol2(pdb_file, mol2_file, is_kekule)


if __name__ == "__main__":
    main()
