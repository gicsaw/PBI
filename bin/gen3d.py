#!/usr/bin/env python
from pbi import ligandtools


def main():

    import argparse
    title_line = 'generate 3d conformation from SMILES'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_smi', required=True,
                        help='input SMILES')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output ligand pdb file')
    parser.add_argument('-m', '--mol_id', required=False, default=None,
                        help='write molecule id to file.')

    args = parser.parse_args()
    smi = args.input_smi
    output_file = args.output_file
    mol_id = args.mol_id

    e = ligandtools.gen_3d(smi, output_file, mol_id)
    if e is not None:
        print(e)


if __name__ == "__main__":
    main()
