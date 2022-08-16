#!/usr/bin/env python
import argparse


def main():

    parser = argparse.ArgumentParser(description='filter_pdb_list')
    parser.add_argument('list_file', type=str,
                        help=' --input_pdb_list list_file')
    parser.add_argument('--include_nmr', action='store_true',
                        required=False, default=False,
                        help='if --include_nmr, use nmr, default=False ')
    parser.add_argument('--resolution', type=float,
                        required=False, default=None,
                        help='--resolution ')
    parser.add_argument('-r', '--residue_range', type=str, required=False,
                        default='-',
                        help='--residue_range 4-100')

    args = parser.parse_args()
    list_file = args.list_file
    include_nmr = args.include_nmr
    resolution_cutoff = args.resolution
    ini_fin = args.residue_range

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    ini, fin = ini_fin.strip().split('-')
    for line in lines:
        lis = line.strip().split(';')
#        pdb_id = lis[0]
        method = lis[1]
        resolution = lis[2]
        chain_positions = lis[3]
        if not include_nmr and method == 'NMR':
            continue
        if resolution_cutoff is not None and resolution.strip() != '-':
            resol = float(resolution.strip().split(' ')[0])
            if resol > resolution_cutoff:
                continue

        chain_position_list = chain_positions.strip().split(',')
        ini_p = None
        fin_p = None
        for chain_position in chain_position_list:
            chain, positions = chain_position.strip().split('=')
            ini_p0, fin_p0 = positions.strip().split('-')
            if ini_p is None:
                ini_p = int(ini_p0)
            elif ini_p > int(ini_p0):
                ini_p = int(ini_p0)
            if fin_p is None:
                fin_p = int(fin_p0)
            elif fin_p < int(fin_p0):
                fin_p = int(fin_p0)

        if ini != '':
            if int(ini) > fin_p:
                continue
        if fin != '':
            if int(fin) < ini_p:
                continue
        print(line[:-1])


if __name__ == "__main__":
    main()
