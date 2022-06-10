#!/usr/bin/env python

def main():

    import argparse
    title_line = 'prepare_gpf.py'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-r', '--receptor_file', type=str, required=True,
                        help='receptor pdbqt file')
    parser.add_argument('-o', '--gpf_file', type=str, required=False,
                        default=None, help='output gpf file')
    parser.add_argument('-n', '--num_grid', type=int, required=False,
                        default=40, help='number of grid')
    parser.add_argument('-s', '--spacing', type=float, required=False,
                        default=0.375, help='spacing (A)')
    parser.add_argument('--smooth', type=float, required=False,
                        default=0.5, help='store minimum energy w/in rad(A)')
    parser.add_argument('-y', '--grid_center', required=False,
                        default=None, help='grid center')
    parser.add_argument('--dielectric', type=float, required=False,
                        default=-0.1465,
                        help='<0, AD4 distance-dep.diel;>0, constant')

    args = parser.parse_args()
    receptor_file = args.receptor_file
    gpf_file = args.gpf_file
    num_grid = args.num_grid
    grid_center = args.grid_center
    smooth = args.smooth
    spacing = args.spacing
    dielectric = args.dielectric

    receptor = '.'.join(receptor_file.strip().split('/')[0].split('.')[:-1])
#    directory = '/'.join(receptor_file.strip().split('/')[0:-1]

    atom_type_list = ['A', 'C', 'NA', 'OA', 'N', 'S', 'SA', 'HD']
    if gpf_file is None:
        gpf_file = receptor + '.gpf'

    fp = open(gpf_file, 'w')

    line = 'npts %d %d %d' % (num_grid, num_grid, num_grid)
    line_out = '%-40s # num.grid points in xyz\n' % line
    fp.write(line_out)

    line = 'gridfld %s.maps.fld' % receptor
    line_out = '%-40s # grid_data_file\n' % line
    fp.write(line_out)

    line = 'spacing %.3f' % spacing
    line_out = '%-40s # spacing(A)\n' % line
    fp.write(line_out)

    line = 'receptor_types A C NA OA N SA HD'
    line_out = '%-40s # receptor atom types\n' % line
    fp.write(line_out)

    line = 'ligand_types %s' % (' '.join(atom_type_list))
    line_out = '%-40s # ligand atom types\n' % line
    fp.write(line_out)

    line = 'receptor %s' % (receptor_file)
    line_out = '%-40s # receptor_file\n' % line
    fp.write(line_out)

    if grid_center is None:
        line = 'gridcenter auto'
    else:
        lis = grid_center.strip().split(',')
        line = 'gridcenter %s %s %s' % (lis[0], lis[1], lis[2])
    line_out = '%-40s # xyz-coordinates or auto\n' % line
    fp.write(line_out)

    line = 'smooth %.3f' % (smooth)
    line_out = '%-40s # store minimum energy w/in rad(A)\n' % line
    fp.write(line_out)

    for atom_type in atom_type_list:
        line = 'map %s.%s.map ' % (receptor, atom_type)
        line_out = '%-40s # atom-specific affinity map\n' % line
        fp.write(line_out)

    line = 'elecmap %s.e.map ' % (receptor)
    line_out = '%-40s # electrostatic potential map\n' % line
    fp.write(line_out)

    line = 'dsolvmap %s.d.map ' % (receptor)
    line_out = '%-40s # desolvation potential map\n' % line
    fp.write(line_out)

    line = 'dielectric %.4f ' % (dielectric)
    line_out = '%-40s # <0, AD4 distance-dep.diel;>0, constant\n' % line
    fp.write(line_out)


if __name__ == "__main__":
    main()
