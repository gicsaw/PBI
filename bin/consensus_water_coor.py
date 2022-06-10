#!/usr/bin/env python
import sys
import os
import numpy as np

usage = '''
water_consensus.py water_list
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
    list_file = sys.argv[1]

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    water_con_avg = list()
    water_con_list = list()

    pdb_dir = 'select'
    water_cutoff = 0.3
    for line in lines:
        code = line[6:11]
        water_file = '%s/%s_HOH.pdb' % (pdb_dir, code)
        if not os.path.exists(water_file):
            continue

        water_dict = read_pdb(water_file)
        water_keys = water_dict.keys()
        for key in water_keys:
            water = water_dict[key]
            water_coor = water[1]
            dist_list = list()

            if len(water_con_avg) == 0:
                water_con_avg += [water_coor]
                water_con_list += [[water_coor]]
                continue
            for water_con_coor in water_con_avg:
                dist = np.linalg.norm(water_con_coor-water_coor)
                dist_list += [dist]
            dist_list = np.array(dist_list)
            j = dist_list.argmin()
            dist_min = dist_list[j]
            if dist_min > water_cutoff:
                water_con_avg += [water_coor]
                water_con_list += [[water_coor]]
            else:
                water_con_list[j] += [water_coor]
                water_con_avg[j] = np.array(water_con_list[j]).mean(axis=0)

    num_con_water = len(water_con_avg)
    for i in range(num_con_water):
        avg_i = water_con_avg[i]
        water_list_i = water_con_list[i]
        num_w = len(water_list_i)
        line_out = '%d %d %.3f %.3f %.3f' % (i, num_w,
                                             avg_i[0], avg_i[1], avg_i[2])
        print(line_out)


if __name__ == "__main__":
    main()
