#!/usr/bin/env python
import sys
import os
from pbi import nwa
from pbi.nwa import Res31

from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue

usage = '''
mut_pdb.py list_file fasta_file pdb_dir
'''


Amino_type = "ARNDCQEGHILKMFPSTWYVX"
Ntype = 21
Amino_dict = dict()
for i in range(0, Ntype):
    amino_acid = Amino_type[i]
    Amino_dict[amino_acid] = i


def read_pdb(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()
    chain_dict = dict()
    uniprot_chain_dict = dict()
    ligand_chain_dict = dict()
    ligand_name_dict = dict()

    check_m = False
    for line in lines:
        if line[0:10] == 'REMARK 465':
            if line[15:27] == 'RES C SSSEQI':
                check_m = True
                continue
            if check_m:
                res_name3 = line[15:18]
                res_idx = int(line[22:26])
                res_inser = line[22:27]
                chain_id = line[19]
                res_name1 = Res31[res_name3]
                if chain_id not in chain_dict:
                    chain_dict[chain_id] = dict()
                chain_dict[chain_id][res_inser] = res_name1
        elif line[0:6] == 'DBREF ':
            chain_id = line[12]
#            uniprot_code = line[33:41].strip()
            uniprot_id = line[42:54].strip()
            if uniprot_id not in uniprot_chain_dict:
                uniprot_chain_dict[uniprot_id] = list()
            uniprot_chain_dict[uniprot_id] += [chain_id]

        elif line[0:6] == 'MODRES':
            #            res_name3_mod = line[12:15]
            chain_id = line[16]
#            res_idx = int(line[18:22])
            res_inser = line[18:23]
            res_name3 = line[24:27]
            res_name1 = Res31[res_name3]
            if chain_id not in chain_dict:
                chain_dict[chain_id] = dict()
            chain_dict[chain_id][res_inser] = res_name1

        elif line[0:6] == 'HET   ':
            ligand_code = line[7:10]
            chain_id = line[12]
            res_idx = int(line[13:17])
            num_atoms = int(line[20:25])
            if chain_id not in ligand_chain_dict:
                ligand_chain_dict[chain_id] = dict()
            if ligand_code not in ligand_chain_dict[chain_id]:
                ligand_chain_dict[chain_id][ligand_code] = list()
            ligand_chain_dict[chain_id][ligand_code] += [res_idx]

            if ligand_code not in ligand_name_dict:
                ligand_name_dict[ligand_code] = dict()
                ligand_name_dict[ligand_code]['name'] = ''
                ligand_name_dict[ligand_code]['num_atoms'] = num_atoms
        elif line[0:6] == 'HETNAM':
            ligand_code = line[11:14]
            lname = line[15:].strip()
            ligand_name_dict[ligand_code]['name'] += lname

        elif line[0:6] == 'ATOM  ':
            res_name3 = line[17:20]
#            res_idx = int(line[22:26])
            res_inser = line[22:27]

            chain_id = line[21]
            if chain_id not in chain_dict:
                chain_dict[chain_id] = dict()

            res_name1 = Res31[res_name3]
            chain_dict[chain_id][res_inser] = res_name1

    return uniprot_chain_dict, chain_dict, ligand_chain_dict, ligand_name_dict


def read_fasta_uniprot(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    fasta_dict = dict()
    for line in lines:
        if line[0] == '>':
            lis = line.strip().split('|')
            key = lis[1]
            fasta_dict[key] = ''
            uniprot_id = lis[2].strip().split()[0]
        else:
            fasta_dict[key] += line.strip()
    return fasta_dict, uniprot_id


def read_fasta_chain(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    fasta_dict = dict()
    for line in lines:
        if line[0] == '>':
            lis = line.strip().split('|')
            chains_line = lis[1]

            chains = list()
            idx = chains_line.find('Chains')

            if idx >= 0:
                chain_id = chains_line[idx+7]
                chains += chain_id
                idx = 8
                while True:
                    idx0 = chains_line[idx:].find(',')
                    if idx0 < 0:
                        break

                    chain_id = chains_line[idx+idx0+2]
                    chains += chain_id
                    idx = idx + idx0 + 3
            else:
                idx = chains_line.find('Chain ')
                chain_id = chains_line[idx+6]
                chains += chain_id

            for chain in chains:
                fasta_dict[chain] = ''
        else:
            for chain in chains:
                fasta_dict[chain] += line.strip()
    return fasta_dict


def creator(q, data, num_sub_proc):
    """
        put data to queue
        input: queue
            data = [(idx1, molid1, smi1), (idx2, molid2, smi2), ...]
            num_sub_proc (for end signal)
    """
    for d in data:
        idx = d[0]
        q.put((idx, d[1]))

    for i in range(0, num_sub_proc):
        q.put('DONE')


def find_m(d, para):

    pdb_code, method, resolution, chain_ini_fin = d
    nw, uniprot_id, seq_uniprot, pdb_dir = para

    line_out = pdb_code
    chain_ids, ini_fin = chain_ini_fin.strip().split('=')

    ini_fin2 = ini_fin.strip().split('-')
    ini = int(ini_fin2[0])
    fin = int(ini_fin2[1])

#    fasta_file = 'pdb/%s.fasta' % (pdb_code)
#    fasta_dict = read_fasta(fasta_file)
    pdb_file = '%s/%s.pdb' % (pdb_dir, pdb_code)
    if not os.path.exists(pdb_file):
        line_out = pdb_file + ' does not exist'
        return line_out
    (uniprot_chain_dict, chain_dict, ligand_chain_dict,
     ligand_name_dict) = read_pdb(pdb_file)
    if uniprot_id not in uniprot_chain_dict:
        line_out =  pdb_file + ' does not have ' + uniprot_id
        return line_out
    chain_list = uniprot_chain_dict[uniprot_id]
    chain_out = ''
    for chain_id in chain_list:
        chain_out += '/%s' % chain_id

    line_out += ';%s' % chain_out.strip('/')

    ligand_list = list()
    for chain_id in ligand_chain_dict:
        if chain_id not in chain_list:
            continue
        chain = ligand_chain_dict[chain_id]
        keys = list(chain.keys())
        ligand_list += keys
    ligand_list = sorted(set(ligand_list))

    ligand_line = ''
    for ligand_id in ligand_list:
        ligand_line += '%s,' % ligand_id
    line_out += ';%s' % ligand_line.strip(',')

    chain_id = chain_list[0]
    res_dict = chain_dict[chain_id]
    res_inser_list = res_dict.keys()
    res_idx_list = list()
    seq_pdb = str()
    inser_dict = dict()
    for resi in res_inser_list:
        res = int(resi[0:4])
        ins = resi[4]
        if ins != ' ':
            if res not in inser_dict:
                inser_dict[res] = ''
            inser_dict[res] += ins
            continue
        res_idx_list += [res]

    res_idx_list = sorted(res_idx_list)
    for idx in res_idx_list:
        resi = '%4d ' % (idx)
        seq_pdb += res_dict[resi]
        if idx in inser_dict:
            inser_list = sorted(inser_dict[idx])
            for inser in inser_list:
                resi = '%4d%s' % (idx, inser)
                seq_pdb += res_dict[resi]

#    start = min(res_idx_list)
#    end = max(res_idx_list)

    seq2 = seq_uniprot[ini-1:fin]
    gg = len(seq_pdb)-len(seq2)
    seq1 = seq_pdb[gg:]

    check = seq1 == seq2
    if check:
        line_out += ';wild_type'
        return line_out
    align, score = nw.NW(seq_pdb, seq_uniprot)

    q3, d3, q4, d4, count = align[0]
#    print(q3, d3)
    s_count = 0
    check_start = False
    check_end = False
    line_mutation = ''
    gap_count = 0
    for i in range(len(d4)):
        if q4[i] == d4[i]:
            if not check_start:
                ini = d3[i]
                check_start = True
                s_count = 0
            s_count += 1

            check_end = False
            gap_count = 0
            fin = d3[i]
        if q4[i] == '-' and d4[i] != '-':
            if not check_end:
                check_end = True
            gap_count += 1
        if gap_count >= 10 and s_count <= 5:
            check_start = False
    if not check_end:
        fin = len(seq_uniprot)
    for i in range(ini-1, fin):
        if q4[i] != d4[i]:
            if q4[i].upper() == 'X':
                continue
            ii = i
            while d3[ii] == '-':
                ii -= 1
            num = d3[ii]
            line_mutation += '%s%s%s,' % (d4[i], num, q4[i].upper())

    line_mutation = line_mutation.strip(',')
    if line_mutation == '':
        line_mutation = 'wild_type'
    line_out += ';%s' % (line_mutation)
    return line_out


def worker(q, para, return_dict):
    """
        generate subprocess for docking
        input
            q (queue)
            return_dict
    """

    nw, uniprot_id, seq_uniprot, pdb_dir = para

#    pid = os.getpid()
    while True:
        qqq = q.get()
        if qqq == 'DONE':
            # print('proc =', os.getpid())
            break

        (idx, d) = qqq
        pdb_code, method, resolution, chain_ini_fin = d

        return_dict[idx] = find_m(d, para)


def main():
    if len(sys.argv) < 3:
        print(usage)
        sys.exit()
#    uniprot_code = sys.argv[1]
#    uniprot_code = 'P00533'
#    print(uniprot_code)
    list_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_dir = 'pdb'
    if len(sys.argv) >= 4:
        pdb_dir = sys.argv[3]

    num_sub_proc = 16
    if len(sys.argv) >= 5:
        num_sub_proc = int(sys.argv[4])

#    list_file = '%s_pdb_list.txt' % uniprot_code
#    fasta_file = '%s.fasta' % uniprot_code
    pdb_list = [x.strip().split(';') for x in open(list_file)]

    fasta_dict, uniprot_id = read_fasta_uniprot(fasta_file)
    uniprot_code = list(fasta_dict.keys())[0]
    seq_uniprot = fasta_dict[uniprot_code]

    nw = nwa.NWalign(g_extend=1.0, g_open2=11.0, lower=1)
    num_data = len(pdb_list)
    num_sub_proc = min(num_sub_proc, num_data)

    data = enumerate(pdb_list)
    q = Queue()
    manager = Manager()
    return_dict = manager.dict()
    proc_master = Process(target=creator,
                          args=(q, data, num_sub_proc))
    proc_master.start()

    para = (nw, uniprot_id, seq_uniprot, pdb_dir)
    procs = []
    for sub_id in range(0, num_sub_proc):
        proc = Process(target=worker, args=(q, para, return_dict))
        procs.append(proc)
        proc.start()

    q.close()
    q.join_thread()
    proc_master.join()
    for proc in procs:
        proc.join()
    keys = sorted(return_dict.keys())

    for key in keys:
        line_out = return_dict[key]
        print(line_out)


if __name__ == "__main__":
    main()
