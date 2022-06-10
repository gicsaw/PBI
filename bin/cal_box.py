#!/usr/bin/env python
import sys
import argparse
from pbi import ligandtools


class LoadFromConfig(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            parser.parse_args(f.read().split(), namespace)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


def main():

    title_line = 'calculate docking box from ligand file'
    parser = argparse.ArgumentParser(description=title_line)
    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('-i', '--autobox_file_list', type=str, action='extend',
                        nargs='+', required=False, default=[],
                        help='--autobox_file_list a.pdb b.pdb\n' +
                             '--autobox_file_list c.pdb')
    parser.add_argument('-m', '--autobox_margin', type=float, required=False,
                        default=4.0, help='margin distance, default: 4.0')
    parser.add_argument('--include_Hs', action='store_false',
                        default=True, help='if --include_Hs: use H atoms')
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('-c', '--cubic_box', action='store_true',
                        default=False,
                        help='if -c: box size is (2*margin,2*margin,2*margin)')

    args = parser.parse_args()
    autobox_file_list = args.autobox_file_list
    margin = args.autobox_margin
    exclude_Hs = args.include_Hs
    cubic_box = args.cubic_box
    if len(autobox_file_list) < 1:
        print('the following arguments are required: --autobox_file_list')
        parser.print_usage()
        sys.exit()

    (cmin, cmax) = ligandtools.cal_box(
        autobox_file_list, exclude_Hs=exclude_Hs)

    center = (cmin + cmax)/2.0
    if cubic_box:
        size = (margin*2, margin*2, margin*2)
    else:
        size = cmax - cmin + margin*2
#    print(cmin, cmax)
#    print(center, size)
    box_center = tuple(center)
    box_size = tuple(size)
    line_out = 'center_x=%.3f\ncenter_y=%.3f\ncenter_z=%.3f\n' % (box_center)
    line_out += 'size_x=%.3f\nsize_y=%.3f\nsize_z=%.3f\n' % (box_size)
    print(line_out.strip())


if __name__ == "__main__":
    main()
