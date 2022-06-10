#!/usr/bin/env python
import sys
import numpy as np
import argparse

usage = '''
select_water.py receptor_file lig_file wat_file
'''


def read_pdb(pdb_file):
    fp = open(pdb_file)
    lines = fp.readlines()
    fp.close()

    atom_dict = dict()
    for line in lines:
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            atom_name = line[12:16]
#            residue_name = line[17:20].strip()
            atom_symbol = atom_name[0:2].strip()
            coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            coor = np.array(coor)
            atom_number = int(line[6:11])
            atom_dict[atom_number] = (atom_symbol, coor, line)

    return atom_dict


def main():
    if len(sys.argv) < 3:
        print(usage)
        sys.exit()
    receptor_file = sys.argv[1]
    ligand_file = sys.argv[2]
    water_file = sys.argv[3]

    parser = argparse.ArgumentParser(description='filter_pdb_list')
    parser.add_argument('-r', '--receptor', type=str, required=True,
                        help=' --input receptor pdb file')
    parser.add_argument('-l', '--ligand', type=str, required=True,
                        help=' --input ligand pdb file')
    parser.add_argument('-w', '--water', type=str, required=True,
                        help=' --input water pdb file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help=' --output water pdb file')
    parser.add_argument('-c', '--dist_cutoff', type=float, required=False,
                        default=3.4, help=' --input_water pdb file')

    args = parser.parse_args()
    receptor_file = args.receptor
    ligand_file = args.ligand
    water_file = args.water
    output_file = args.out
    dist_cutoff = args.dist_cutoff

    receptor_dict = read_pdb(receptor_file)
    ligand_dict = read_pdb(ligand_file)
    water_dict = read_pdb(water_file)

    ligand_keys = ligand_dict.keys()

    receptor_keys = receptor_dict.keys()

    water_keys = water_dict.keys()

    fp_out = open(output_file, 'w')
    for water_num in water_keys:
        water = water_dict[water_num]
        water_coor = water[1]
        check_bind_ligand = False
        check_bind_receptor = False

        for ligand_num in ligand_keys:
            ligand = ligand_dict[ligand_num]
#            ligand_atom_symbol = ligand[0]
            ligand_coor = ligand[1]
            dist = np.linalg.norm(water_coor-ligand_coor)
            if dist <= dist_cutoff:
                check_bind_ligand = True
                break

        if not check_bind_ligand:
            continue

        for receptor_num in receptor_keys:
            receptor = receptor_dict[receptor_num]
#            receptor_atom_symbol = receptor[0]
            receptor_coor = receptor[1]
            dist = np.linalg.norm(water_coor-receptor_coor)
            if dist <= dist_cutoff:
                check_bind_receptor = True
                break

        if not check_bind_receptor:
            continue
        line_out = water[2]
        fp_out.write(line_out)
    fp_out.close()


if __name__ == "__main__":
    main()
