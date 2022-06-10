#!/usr/bin/env python
import sys
import numpy as np

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
    if len(sys.argv) < 2:
        print(usage)
        sys.exit()
    water_file = sys.argv[1]
    consensus_file = sys.argv[2]

    water_dict = read_pdb(water_file)
    fp = open(consensus_file)
    lines = fp.readlines()
    fp.close()
    consensus_coor_list = list()
    for line in lines:
        lis = line.strip().split()
        coor = np.array((float(lis[2]), float(lis[3]), float(lis[4])))
        consensus_coor_list += [coor]

    dist_cutoff = 0.3
    water_keys = water_dict.keys()
    for water_num in water_keys:
        water = water_dict[water_num]
        water_coor = water[1]
        for consensus_coor in consensus_coor_list:
            dist = np.linalg.norm(water_coor - consensus_coor)
            if dist <= dist_cutoff:
                line_out = water[2][:-1]
                print(line_out)
                break


if __name__ == "__main__":
    main()
