#!/usr/bin/env python
from pbi.pdbtools import PDBtools


def main():

    import argparse
    title_line = 'split ligand in pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-r', '--ref_file', type=str, required=False,
                        default=None, help='reference protein pdb file')
    parser.add_argument('-t', '--tmp_dir', type=str, required=False,
                        default='.', help='temp directory, defalut=.')
    parser.add_argument('-d', '--out_dir', type=str, required=False,
                        default='.', help='output directory')
    parser.add_argument('-e', '--exclude', type=str, nargs='+',
                        required=False, default=list(),
                        help='exclude ligand name in PDB, ex: -e HOH AAA' +
                        ', default: no exclude ligand ')

    args = parser.parse_args()
    input_file = args.input_file
    ref_file = args.ref_file
    is_align = False
    if ref_file is not None:
        is_align = True
    tmp_dir = args.tmp_dir

    out_dir = args.out_dir
    exclude_list = args.exclude

    pdb_code = input_file.split('/')[-1].split('.')[0].split('_')[0]

    result = PDBtools.read_pdb_protein(input_file, remain_remark=True)
    (protein_chain_residues, pdb_info_lines, protein_dict,
     ligand_dict, conect_dict) = result

    if is_align:
        result = PDBtools.align_3d(input_file, protein_dict,
                                   ligand_dict, ref_file, tmp_dir)
        (protein_dict, ligand_dict) = result

    model_dict = dict()
    model_dict[0] = (pdb_info_lines, protein_dict, conect_dict)
    receptor_file = '%s/%s_receptor.pdb' % (out_dir, pdb_code)
    PDBtools.write_model_pdb(model_dict, receptor_file)

    ligand_mol_dict, conect_mol_dict = PDBtools.split_ligand(
        ligand_dict, conect_dict)
    mol_code_list = ligand_mol_dict.keys()

    pdb_info_lines0 = list()
    for mol_code in mol_code_list:
        if mol_code in exclude_list:
            continue
        ligand_dict0 = ligand_mol_dict[mol_code]
        conect_dict0 = conect_mol_dict[mol_code]

        model_dict = dict()
        model_dict[0] = (pdb_info_lines0, ligand_dict0, conect_dict0)
        ligand_file = '%s/%s_%s.pdb' % (out_dir, pdb_code, mol_code)
        PDBtools.write_model_pdb(model_dict, ligand_file)


if __name__ == "__main__":
    main()
