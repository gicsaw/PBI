#!/usr/bin/env python
from pbi.pdbtools import PDBtools


def main():

    import argparse
    title_line = 'Fixer for protein pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdb', help='output protein pdb file')

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    out_format = output_file.split('.')[-1]
    if out_format == 'pdbqt':
        tmp_file = '.'.join(output_file.split('.')[:-1]) + '.pdb'
    else:
        tmp_file = output_file

    PDBtools.fix_protein(input_file, tmp_file,
                         add_missing_residue=True, pH=7.4)
    if out_format == 'pdbqt':
        e = PDBtools.protein_to_pdbqt(tmp_file, output_file)
        if e is not None:
            print(e)


if __name__ == "__main__":
    main()
