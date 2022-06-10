#!/usr/bin/env python
import sys
from pbi.pdbtools import PDBtools


def main():

    import argparse
    title_line = 'Fixer for ligand pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input ligand pdb file')
    parser.add_argument('-r', '--ref_file', type=str, required=True,
                        help='ref cif file from https://files.rcsb.org/ligands')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdb',  help='output ligand pdb file')
    parser.add_argument('-n', '--neutralize', action='store_true',
                        required=False, help='neutralize mol ')
    parser.add_argument('-p', '--pH', type=float, required=False, default=None,
                        help='protonation state')
    parser.add_argument('--is_peptide', action='store_true',
                        required=False, help='is peptide?')
    parser.add_argument('--fix_atom_idx', action='store_true',
                        required=False, help='fix_atom_idx')

    args = parser.parse_args()
    input_file = args.input_file
    ref_file = args.ref_file
    output_file = args.output_file
    neutralize = args.neutralize
    pH = args.pH
    is_peptide = args.is_peptide

    is_fix_atom_idx = args.fix_atom_idx
    out_format = output_file.split('.')[-1]
    ref_format = ref_file.split('.')[-1]
    if ref_format != 'cif':
        print('reference file format is not cif')
        sys.exit()
    if out_format == 'pdbqt':
        tmp_file = '.'.join(output_file.split('.')[:-1]) + '.pdb'
    else:
        tmp_file = output_file
    PDBtools.fix_ligand_ref(input_file, ref_file, output_file,
                            neutralize=neutralize,
                            pH=pH, is_peptide=is_peptide,
                            is_fix_atom_idx=is_fix_atom_idx)
    if out_format == 'pdbqt':
        PDBtools.ligand_to_pdbqt(tmp_file, output_file)


if __name__ == "__main__":
    main()
