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
    parser.add_argument('-r', '--random_coor', required=False, default=False,
                        action='store_true',
                        help='use random coordinates, defalut=False.')

    parser.add_argument('-s', '--use_simple', required=False, default=False,
                        action='store_true',
                        help='use short option, defalut=False.')

    parser.add_argument('-u', '--use_uff', required=False, default=False,
                        action='store_true',
                        help='use UFFOptimization, defalut=False.')


    args = parser.parse_args()
    smi = args.input_smi
    output_file = args.output_file
    mol_id = args.mol_id
    random_coor = args.random_coor
    use_simple = args.use_simple
    use_uff = args.use_uff


    e = ligandtools.gen_3d(smi, output_file, mol_id, random_coor=random_coor,
                           use_simple=use_simple, use_uff=use_uff)
    if e is not None:
        print(e)


if __name__ == "__main__":
    main()
