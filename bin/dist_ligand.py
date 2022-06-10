#!/usr/bin/env python
import numpy as np
from pbi import ligandtools

def main():

    import argparse
    title_line = 'dist_ligand'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-r', '--ref_file', type=str, required=True,
                        help='reference protein pdb file')
    parser.add_argument('-c', '--radius_cut', type=float,
                        required=False, default=None,
                        help='exclude ligand which radius is smaller than cutoff.' +
                        ' default is None')
    parser.add_argument('--include_Hs', action='store_false',
                        default=True, help='if --include_Hs: use H atoms')

    args = parser.parse_args()
    input_file = args.input_file
    ref_file = args.ref_file
    exclude_Hs = args.include_Hs
    radius_cut = args.radius_cut

    ref_model_dict = ligandtools.read_coor_pdb(
        ref_file, exclude_Hs=exclude_Hs)
    ref_dict = ref_model_dict[1]
    keys = sorted(ref_dict.keys())
    ref_cmin, ref_cmax = ligandtools.cal_ligand_size(ref_dict[keys[0]])
    ref_center = (ref_cmin + ref_cmax)/2.0
    ref_size = ref_cmax - ref_cmin
    ref_radius = np.linalg.norm(ref_size/2)

    ligand_model_dict = ligandtools.read_coor_pdb(
        input_file, exclude_Hs=exclude_Hs)
    ligand_dict = ligand_model_dict[1]
    keys = sorted(ligand_dict.keys())
    for key in keys:
        cmin, cmax = ligandtools.cal_ligand_size(ligand_dict[key])
        center = (cmin + cmax)/2.0
        size = cmax - cmin
        radius = np.linalg.norm(size/2)
        if radius_cut is not None:
            if radius < radius_cut:
                continue
        dist = np.linalg.norm(center - ref_center)
        rr = radius + ref_radius

        if dist < rr:
            line_out = '%s %d %.3f %.3f' % (input_file, key, dist, rr)
            print(line_out)


if __name__ == "__main__":
    main()
