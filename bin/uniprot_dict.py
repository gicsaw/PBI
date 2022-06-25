#!/usr/bin/env python
import sys
import os
import gzip
import pickle
from urllib import request

usage = '''
uniprot_dict.py make_pickle
uniprot_dict.py search uniprot_ID
example: uniprot_dict.py search EGFR_HUMAN
uniprot_dict.py pdb uniprot_ID ini-fin
example:  uniprot_dict.py pdb EGFR_HUMAN 650-

'''


def read_uniprot(uniprot_file):
    fp = gzip.open(uniprot_file, 'rb')

    uniprot_dict = dict()
#    check_signal = False
    check_chain = False
    check_domain = False
    check_region = False
    check_seq = False
    for i, line_b in enumerate(fp):
        line = line_b.decode('ascii')
        #        if i>100000:
        #            break
        if line[0:2] == "ID":
            uniprot_ID = line[5:29].strip()
            name, OS0 = uniprot_ID.split('_')
#            if OS0 !='HUMAN' :
#                continue
            uniprot_dict[uniprot_ID] = dict()
            uniprot_dict[uniprot_ID]['name'] = name
            uniprot_dict[uniprot_ID]['PDB'] = []
            uniprot_dict[uniprot_ID]['BindingDB'] = []
            uniprot_dict[uniprot_ID]['ChEMBL'] = []
            uniprot_dict[uniprot_ID]['DrugBank'] = []
            uniprot_dict[uniprot_ID]['signal'] = []
            uniprot_dict[uniprot_ID]['chain'] = []
            uniprot_dict[uniprot_ID]['domain'] = []
            uniprot_dict[uniprot_ID]['region'] = []
            uniprot_dict[uniprot_ID]['gene_id'] = []

            sequence = ''
            keyword = ''
            ll = line[30:].strip().rstrip('.').split('; ')
#            is_reviewed = ll[0]
            seq_len = ll[1].strip()
            uniprot_dict[uniprot_ID]['sequence_length'] = seq_len
            continue

#        if OS0 !='HUMAN' :
#            continue
        if line[0:2] == 'AC':
            ll = line[5:].strip().split('; ')
            uniprot_code = ll[0].strip().strip(';')
            if 'uniprot_code' not in uniprot_dict[uniprot_ID]:
                uniprot_dict[uniprot_ID]['uniprot_code'] = uniprot_code
            continue

        if line[0:2] == 'GN':
            j = line.find('Name=')
            if j != -1:
                gene_name = line[j+5:].strip().split(';')[0].split(' ')[0]
                uniprot_dict[uniprot_ID]['gene_name'] = gene_name
            j = line.find('Synonyms=')
            if j != -1:
                synonyms = line[j+9:].strip().split(';')[0].split(' ')[0]
                uniprot_dict[uniprot_ID]['synonyms'] = synonyms
            j = line.find('ORFNames=')
            if j != -1:
                orfnames = line[j+9:].strip().split(';')[0].split(' ')[0]
                uniprot_dict[uniprot_ID]['ORFNames'] = orfnames
            continue

        if line[0:2] == 'OS':
            OS = line[5:].strip().rstrip('.')
            uniprot_dict[uniprot_ID]['OS'] = OS
            continue

        if line[0:2] == 'DR':
            ll = line[5:].strip().rstrip('.').split('; ')
            if ll[0] == 'PDB':
                uniprot_dict[uniprot_ID]['PDB'] += [';'.join(ll[1:])]
            elif ll[0] == 'BindingDB':
                uniprot_dict[uniprot_ID]['BindingDB'] = ll[1]
            elif ll[0] == 'ChEMBL':
                uniprot_dict[uniprot_ID]['ChEMBL'] = ll[1]
            elif ll[0] == 'DrugBank':
                uniprot_dict[uniprot_ID]['DrugBank'] += [ll[1:]]
            elif ll[0] == 'GeneID':
                uniprot_dict[uniprot_ID]['gene_id'] += ll[1:]
            continue

        if line[0:2] == 'KW':
            keyword += ' ' + line[5:].strip()
            continue

        if line[0:2] == 'FT':
            if line[5:20].strip() != '':
                check_chain = False
#                check_signal = False
                check_domain = False
                check_region = False
            elif line[5:11] == 'SIGNAL':
                inifin = line[21:].strip()
#                check_signal = True
                uniprot_dict[uniprot_ID]['signal'] += [inifin]
            elif line[5:10] == 'CHAIN':
                inifin = line[21:].strip()
                check_domain = True
            elif line[5:11] == 'DOMAIN':
                inifin = line[21:].strip()
                check_domain = True
            elif line[5:11] == 'REGION':
                inifin = line[21:].strip()
                check_region = True

            elif line[5:11].strip() == '' and check_chain:
                if line[22:26] == 'note':
                    note = line[27:].strip()
                    uniprot_dict[uniprot_ID]['chain'] += [[inifin, note]]
            elif line[5:11].strip() == '' and check_domain:
                if line[22:26] == 'note':
                    note = line[27:].strip()
                    uniprot_dict[uniprot_ID]['domain'] += [[inifin, note]]
                j = line[20:].find('Protein kinase')
                if j != -1:
                    uniprot_dict[uniprot_ID]['kinase_domain'] = inifin
            elif line[5:11].strip() == '' and check_region:
                if line[22:26] == 'note':
                    note = line[27:].strip()
                    uniprot_dict[uniprot_ID]['region'] += [[inifin, note]]
            continue
        if line[0:2] == 'SQ':
            check_seq = True
            continue
        if line[0:2] == '  ' and check_seq:
            sequence += ''.join(line[5:].strip().split())
            continue
        if line[0:2] == '//':
            check_seq = False
            uniprot_dict[uniprot_ID]['sequence'] = sequence
            uniprot_dict[uniprot_ID]['keyword'] = keyword.lstrip().rstrip(
                '.').split('; ')
            if 'uniprot_code' not in uniprot_dict[uniprot_ID]:
                uniprot_dict[uniprot_ID]['uniprot_code'] = ''
            continue

    fp.close()
    return uniprot_dict


def read_human(uniprot_dict):

    keys = sorted(uniprot_dict.keys())
    uniprot_human = dict()
    for uniprot_ID in keys:
        name, OS0 = uniprot_ID.split('_')
        if OS0 != 'HUMAN':
            continue
        uniprot_human[uniprot_ID] = uniprot_dict[uniprot_ID]
    return uniprot_human


def main():
    if len(sys.argv) < 2:
        print(usage)
        sys.exit()

    db_dir = '.'
    db_con_file = './db_dir.txt'
    if not os.path.exists(db_con_file):
        print('does not found db_dir.txt')
        print('default is db_dir=.')
        print('make db_dir.txt')
        fp = open('db_dir.txt', 'w')
        fp.write('db_dir=.')
        fp.close()

    fp = open(db_con_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        lis = line.strip().split('=')
        if lis[0].strip() == 'db_dir':
            db_dir = lis[1].strip()

    uniprot_file = db_dir + '/uniprot_sprot.dat.gz'
    pickle_file = db_dir + '/uniprot.pkl'
    pickle_human_file = db_dir + '/uniprot_human.pkl'

    if sys.argv[1] == 'make_pickle':
        uniprot_dict = read_uniprot(uniprot_file)
        with open(pickle_file, 'wb') as f:
            pickle.dump(uniprot_dict, f)

        uniprot_human = read_human(uniprot_dict)
        with open(pickle_human_file, 'wb') as f:
            pickle.dump(uniprot_human, f)
        sys.exit()
    elif sys.argv[1] == 'search' or sys.argv[1].lower() == 'pdb':
        if len(sys.argv) < 3:
            print('uniprot_ID?')
            sys.exit()
        uniprot_ID = sys.argv[2]
        name, OS0 = uniprot_ID.strip().split('_')
        if OS0 == 'HUMAN':
            with open(pickle_human_file, 'rb') as f:
                uniprot_dict = pickle.load(f)
        else:
            with open(pickle_file, 'rb') as f:
                uniprot_dict = pickle.load(f)
    else:
        print('strange option')
        print(usage)
        sys.exit()

    uniprot = uniprot_dict[uniprot_ID]
    uniprot_code = uniprot['uniprot_code']
    gene_id = uniprot['gene_id']
    gene_id0 = '-'
    gene_name = '-'
    if 'gene_name' in uniprot:
        gene_name = uniprot['gene_name']
    if len(gene_id) > 1:
        gene_id0 = gene_id[0]
#    sequence = uniprot['sequence']

    print(uniprot_code, name, gene_name, gene_id0)

    fasta_file = '%s.fasta' % uniprot_code
    url = 'https://www.uniprot.org/uniprot/%s.fasta' % uniprot_code
    with request.urlopen(url) as r:
        lines = r.read().decode('utf-8')
        with open(fasta_file, 'w') as fp:
            fp.write(lines)

    if sys.argv[1].lower() == 'pdb':
        if len(sys.argv) == 3:
            ini_fin = '-'
        elif len(sys.argv) > 3:
            ini_fin = sys.argv[3]
        ini, fin = ini_fin.strip().split('-')
        pdb_list_file = '%s_pdb_list.txt' % uniprot_code
        fp = open(pdb_list_file, 'w')
        pdb_list = uniprot['PDB']
        for pdb_line in pdb_list:
            lis = pdb_line.strip().split(';')
#            pdb_id = lis[0]
#            method = lis[1]
#            resolution = lis[2]
            chain_positions = lis[3]
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
#            print(ini, fin, ini_p, fin_p)
            line_out = pdb_line + '\n'
            fp.write(line_out)
        fp.close()


if __name__ == "__main__":
    main()
