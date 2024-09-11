#!/usr/bin/env python
from openbabel import pybel
from pbi.pdbtools import PDBtools


def gen_bond_info(conect):
    bond_info = dict()
    for i in conect.keys():
        bond_info[i] = dict()
        b = conect[i]
        for j in b:
            if j not in bond_info[i]:
                bond_info[i][j] = 1
            else:
                bond_info[i][j] += 1
    return bond_info


def mol2_to_kekule(mol2_line_list, bond_info):
    mol2_line_list_kekule = list()
    bond_check = False
    for line in mol2_line_list:

        if line == '':
            mol2_line_list_kekule.append(line)
            continue
        if line.startswith('@'):
            if line == '@<TRIPOS>BOND':
                bond_check = 1
            else:
                bond_check = 0
            mol2_line_list_kekule.append(line)
            continue

        if not bond_check:
            mol2_line_list_kekule.append(line)
            continue
#        bond_idx = int(line[0:6])
        atom1_idx = int(line[6:12])
        atom2_idx = int(line[12:18])
        bond_name = line[18:23]

        if bond_name == '   ar':
            bond_new = bond_info[atom1_idx][atom2_idx]
            line_new = '%s%5d' % (line[0:18], bond_new)
            mol2_line_list_kekule.append(line_new)
        else:
            mol2_line_list_kekule.append(line)
    return mol2_line_list_kekule


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

    m_pdb = PDBtools.read_pdb_ligand(pdb_file)
    p_keys = sorted(m_pdb.keys())

    ms = pybel.readfile('pdb', pdb_file)
    fp_out = open(mol2_file, 'w')

    for model_idx, m in enumerate(ms):
        mol2_line = m.write('mol2')
        mol2_line_list = mol2_line.split('\n')

        if is_kekule:
            conect = m_pdb[p_keys[model_idx]][2]
            bond_info = gen_bond_info(conect)
            mol2_line_list_kekule = mol2_to_kekule(mol2_line_list, bond_info)
            mol2_line_list = mol2_line_list_kekule
        for line in mol2_line_list:
            line += '\n'
            fp_out.write(line)

    fp_out.close()


if __name__ == "__main__":
    main()
