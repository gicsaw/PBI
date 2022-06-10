#!/usr/bin/env python
from pbi import ligandtools


def main():

    import argparse
    title_line = 'convert pdbqt to pdb using reference pdb file'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', required=True,
                        help='input ligand pdbqt file')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output ligand pdb file')
    parser.add_argument('-r', '--ref_file', required=True,
                        help='reference ligand pdb file')

    args = parser.parse_args()
    ligand_input_file = args.input_file
    ligand_output_file = args.output_file
    ref_file = args.ref_file

    e = ligandtools.pdbqt_to_pdb_ref(ligand_input_file, ligand_output_file,
                                      ref_file)
    if e is not None:
        print(e)


if __name__ == "__main__":
    main()
