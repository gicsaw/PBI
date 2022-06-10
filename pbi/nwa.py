#!/usr/bin/env python
import sys
import numpy as np

sys.setrecursionlimit(10000)

USAGE = """
NWalign.py list_file
"""

Res31 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
         'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
         'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
         'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
         'ASX': 'N', 'GLX': 'Q', 'UNK': 'X', 'INI': 'K', 'AAR': 'R',
         'ACE': 'X', 'ACY': 'G', 'AEI': 'T', 'AGM': 'R', 'ASQ': 'D',
         'AYA': 'A', 'BHD': 'D', 'CAS': 'C', 'CAY': 'C', 'CEA': 'C',
         'CGU': 'E', 'CME': 'C', 'CMT': 'C', 'CSB': 'C', 'CSD': 'C',
         'CSE': 'C', 'CSO': 'C', 'CSP': 'C', 'CSS': 'C', 'CSW': 'C',
         'CSX': 'C', 'CXM': 'M', 'CYG': 'C', 'CYM': 'C', 'DOH': 'D',
         'EHP': 'F', 'FME': 'M', 'FTR': 'W', 'GL3': 'G', 'H2P': 'H',
         'HIC': 'H', 'HIP': 'H', 'HTR': 'W', 'HYP': 'P', 'KCX': 'K',
         'LLP': 'K', 'LLY': 'K', 'LYZ': 'K', 'M3L': 'K', 'MEN': 'N',
         'MGN': 'Q', 'MHO': 'M', 'MHS': 'H', 'MIS': 'S', 'MLY': 'K',
         'MLZ': 'K', 'MSE': 'M', 'NEP': 'H', 'NPH': 'C', 'OCS': 'C',
         'OCY': 'C', 'OMT': 'M', 'OPR': 'R', 'PAQ': 'Y', 'PCA': 'Q',
         'PHD': 'D', 'PRS': 'P', 'PTH': 'Y', 'PYX': 'C', 'SEP': 'S',
         'SMC': 'C', 'SME': 'M', 'SNC': 'C', 'SNN': 'D', 'SVA': 'S',
         'TPO': 'T', 'TPQ': 'Y', 'TRF': 'W', 'TRN': 'W', 'TRO': 'W',
         'TYI': 'Y', 'TYN': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y',
         'YOF': 'Y', 'FOR': 'X', '---': '-', 'PTR': 'Y', 'LCX': 'K',
         'SEC': 'D', 'MCL': 'K', 'LDH': 'K'}

Res20 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

Res13 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
         'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
         'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
         'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
         'X': 'UNK'}

Amino_type = "ARNDCQEGHILKMFPSTWYVX"
Ntype = 21
Amino_dict = dict()
for i in range(0, Ntype):
    amino_acid = Amino_type[i]
    Amino_dict[amino_acid] = i

Blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}


class NWalign(object):
    def __init__(self, scoring_matrix=Blosum62, g_extend=1.0, g_open2=11.0, lower=0):
        """
            lower : if 0, None
                    elif 1, q3 is lower where q3!=d3
                    elif 2, q3 and d3 are lower where q3!=d3
        """

        self.scoring_matrix = scoring_matrix
        self.g_extend = g_extend
        self.g_open = g_open2 - g_extend
        self.lower = lower

    def g(self, l):
        return self.g_open + (l)*self.g_extend

    def R(self, chr1, chr2):

        #    if chr1==chr2 :
        #        return 1
        #    else :
        #        return 0
        if (chr1, chr2) in self.scoring_matrix:
            key = (chr1, chr2)
            score = self.scoring_matrix[key]
        elif (chr2, chr1) in self.scoring_matrix:
            key = (chr2, chr1)
            score = self.scoring_matrix[key]
        else:
            score = -9

        return score

    def cal_score_matrix(self, q, d):
        q = ' ' + q
        d = ' ' + d
        len_q = len(q)
        len_d = len(d)

        H = np.zeros([len_q, len_d], float)
        E = np.zeros([len_q, len_d], float)
        F = np.zeros([len_q, len_d], float)
        G = np.zeros([len_q, len_d], float)
        for i in range(1, len_q):
            H[i][0] = -self.g(i)
            E[i][0] = -self.g(i)
            F[i][0] = -self.g(i)
            G[i][0] = -self.g(i)

        for j in range(1, len_d):
            H[0][j] = -self.g(j)
            E[0][j] = -self.g(j)
            F[0][j] = -self.g(j)
            G[0][j] = -self.g(j)

        for i in range(1, len_q):
            for j in range(1, len_d):
                E[i][j] = max(E[i-1][j] - self.g_extend,
                              F[i-1][j] - self.g_open - self.g_extend,
                              G[i-1][j] - self.g_open - self.g_extend)
                F[i][j] = max(E[i][j-1] - self.g_open - self.g_extend,
                              F[i][j-1] - self.g_extend,
                              G[i][j-1] - self.g_open - self.g_extend)
                G[i][j] = max(E[i-1][j-1] + self.R(q[i], d[j]),
                              F[i-1][j-1] + self.R(q[i], d[j]),
                              G[i-1][j-1] + self.R(q[i], d[j]))
                H[i][j] = max(E[i][j], F[i][j], G[i][j])

        score = H[len_q-1][len_d-1]

        return E, F, G, H, score

    def backtrack(self, i, j, k, H, E, F, G, q, d, q2, d2):
        q3 = []
        d3 = []
        if(i == 0):
            while(j > 0):
                q2[k] = '-'
                d2[k] = j
                k = k+1
                j = j-1
        elif(j == 0):
            while(i > 0):
                q2[k] = i
                d2[k] = '-'
                k = k+1
                i = i-1
        else:
            if(H[i][j] == G[i][j]):
                q2[k] = i
                d2[k] = j
                self.backtrack(i-1, j-1, k+1, H, E, F, G, q, d, q2, d2)
            if(H[i][j] == E[i][j]):
                l = 0
                while H[i][j] == (E[i-l][j]-l*self.g_extend) and (i > l):
                    q2[k+l] = i-l
                    d2[k+l] = '-'
                    l = l+1
                self.backtrack(i-l, j, k+l, H, E, F, G, q, d, q2, d2)
            if(H[i][j] == F[i][j]):
                l = 0
                while H[i][j] == (F[i][j-l]-l*self.g_extend) and (j > l):
                    q2[k+l] = '-'
                    d2[k+l] = j-l
                    l = l+1
                self.backtrack(i, j-l, k+l, H, E, F, G, q, d, q2, d2)

        if(i == 0 and j == 0):
            for i in range(1, k):
                q3 += [q2[k-i]]
                d3 += [d2[k-i]]
            self.align.append((q3, d3))
        return

    def dp(self, q, d, E, F, G, H):
        len_q = len(q)+1
        len_d = len(d)+1
        q2 = [' ']*(len_q+len_d)
        d2 = [' ']*(len_q+len_d)
        k = 1
        self.backtrack(len_q-1, len_d-1, k, H, E, F, G, q, d, q2, d2)
        return

    def num_to_seq(self, q, d, q3, d3):
        q4 = ""
        d4 = ""
        count = 0
        for i in range(len(q3)):
            qi = q3[i]
            di = d3[i]
            if qi == "-":
                qc = qi
            else:
                qc = q[qi-1]

            if di == "-":
                dc = di
            else:
                dc = d[di-1]
            if qc != dc:
                if self.lower >= 1:
                    qc = qc.lower()
                if self.lower >= 2:
                    dc = dc.lower()
            q4 += qc
            d4 += dc
            if qc == dc:
                count += 1
        return q4, d4, count

    def NW(self, q, d):

        self.align = []
        E, F, G, H, score = self.cal_score_matrix(q, d)
        self.dp(q, d, E, F, G, H)

        align_seq = []
        for q3, d3 in self.align:
            q4, d4, count = self.num_to_seq(q, d, q3, d3)
            align_seq.append((q3, d3, q4, d4, count))
        return align_seq, score


def read_fa(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    seq = ""
    check = 0
    for line in lines:
        if line.startswith(">"):
            if check == 1:
                break
            check = 0
            header = line.strip()
            continue
        seq += line.strip()
    return seq, header


def find_query_aligned_sequence(q4, d4):

    new_fasta = ''
    num_aa = len(q4)
    k = -1
    st = 0
    if d4[0] != '-' and q4[0] != '-' and d4[1:4] == '---':
        st = 1
        k += st
    start_check = False
    for i in range(st, num_aa):
        if q4[i] == '-':
            continue
        k += 1
        if d4[i] == '-':
            continue
        if not start_check:
            ini = k+1
            start_check = True

        new_fasta += q4[i]
        fin = k+1

    return new_fasta, ini, fin


def main():

    if len(sys.argv) < 3:
        print('nwalign.py query.fasta template.fasta')
        sys.exit()

    file_q = sys.argv[1]
    file_d = sys.argv[2]

    q, q_header = read_fa(file_q)
    d, d_header = read_fa(file_d)

    nwalign = NWalign(scoring_matrix=Blosum62, g_extend=1.0, g_open2=11.0)
    align, score = nwalign.NW(q, d)
#    print(len(align))
#    for q3,d3,q4,d4,count in align:
#        print(count/len(q),count/len(d),count,len(q),len(d))
#        print(q4)
#        print(d4)

    q3, d3, q4, d4, count = align[0]
    print(q4)
    print(d4)

#    query_aligned_sequence, ini, fin = find_query_aligned_sequence(q4, d4)
#    line_out = '>query:%s|template:%s|%d-%d' % (file_q, file_d, ini, fin)
#    print(line_out)
#    print(query_aligned_sequence)


if __name__ == '__main__':
    main()
