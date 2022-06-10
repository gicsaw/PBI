#!/usr/bin/env python
import sys
from pbi.pdbtools import PDBtools


def main():

    import argparse
    title_line = 'convert pdb to pdbqt'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdbqt', help='output protein pdb file')
    parser.add_argument('-r', '--receptor', action='store_true',
                        default=False, help='is receptor?')
    parser.add_argument('-l', '--ligand', action='store_true',
                        default=False, help='is ligand?')

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    receptor = args.receptor
    ligand = args.ligand
    if (receptor and ligand) or (not receptor and not ligand):
        print('receptor or ligand?')
        sys.exit()

    if receptor:
        e = PDBtools.protein_to_pdbqt(input_file, output_file)
        if e is not None:
            print(e)
    if ligand:
        e = PDBtools.ligand_to_pdbqt(input_file, output_file)
        if e is not None:
            print(e)


if __name__ == "__main__":
    main()
