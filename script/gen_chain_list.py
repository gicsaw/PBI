#!/usr/bin/env python
import sys
import os
import argparse


def main():

    parser = argparse.ArgumentParser(description='gen_chain_list')
    parser.add_argument('-i', '--input_list', type=str, required=False,
                        default='list_mutation.txt',
                        help=' --input_list list_mutation.txt')
    parser.add_argument('-o', '--output_list', type=str, required=False,
                        default='list_chain.txt',
                        help=' --output_list list_chain.txt')

    parser.add_argument('--use_all_chain', action='store_true',
                        required=False, default=False,
                        help='if --use_all_chain, use all chain,' +
                        'else use first chain. default=False ')

    args = parser.parse_args()
    list_file = args.input_list
    output_file = args.output_list

    use_all_chain = args.use_all_chain

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    fp_out = open(output_file,'w')
    for line in lines:
        lis = line.strip().split(';')
        pdb_code = lis[0]
        chain_list = lis[1].strip().split('/')
        if use_all_chain is True:
            for chain_id in chain_list:
                line_out = '%s %s\n' %(pdb_code, chain_id)
                fp_out.write(line_out)
        else:
            chain_id = chain_list[0]
            line_out = '%s %s\n' %(pdb_code, chain_id)
            fp_out.write(line_out)

if __name__ == "__main__":
    main()
