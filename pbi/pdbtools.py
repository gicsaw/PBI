#!/usr/bin/env python
import os
import copy
#from simtk.openmm import app
from openmm import app
import pdbfixer
import subprocess
import numpy as np


class PDBtools(object):
    '''
        This is pdbtools for fucking pdb file.
        rotation (to reference) with TMalign
        split chains
        split protein and ligand
        residue re-indexingfix_protein(protein_pdb_file, protein_fix_file)
        fix missing atoms of protein
        and so on...
    '''
    metal_atom_list = ['LI  ', 'MG  ', ' K  ', 'MN  ', 'CA  ', 'FE  ', 'ZN  ']
    metal_residue_list = ['LI', 'MG', 'K', 'MN', 'CA', 'FE', 'ZN']
    substitutions = {'2AS': 'ASP', '3AH': 'HIS', '5HP': 'GLU', 'ACL': 'ARG',
                     'AGM': 'ARG', 'AIB': 'ALA', 'ALM': 'ALA', 'ALO': 'THR',
                     'ALY': 'LYS', 'ARM': 'ARG', 'ASA': 'ASP', 'ASB': 'ASP',
                     'ASK': 'ASP', 'ASL': 'ASP', 'ASQ': 'ASP', 'AYA': 'ALA',
                     'BCS': 'CYS', 'BHD': 'ASP', 'BMT': 'THR', 'BNN': 'ALA',
                     'BUC': 'CYS', 'BUG': 'LEU', 'C5C': 'CYS', 'C6C': 'CYS',
                     'CAS': 'CYS', 'CCS': 'CYS', 'CEA': 'CYS', 'CGU': 'GLU',
                     'CHG': 'ALA', 'CLE': 'LEU', 'CME': 'CYS', 'CSD': 'ALA',
                     'CSO': 'CYS', 'CSP': 'CYS', 'CSS': 'CYS', 'CSW': 'CYS',
                     'CSX': 'CYS', 'CXM': 'MET', 'CY1': 'CYS', 'CY3': 'CYS',
                     'CYG': 'CYS', 'CYM': 'CYS', 'CYQ': 'CYS', 'DAH': 'PHE',
                     'DAL': 'ALA', 'DAR': 'ARG', 'DAS': 'ASP', 'DCY': 'CYS',
                     'DGL': 'GLU', 'DGN': 'GLN', 'DHA': 'ALA', 'DHI': 'HIS',
                     'DIL': 'ILE', 'DIV': 'VAL', 'DLE': 'LEU', 'DLY': 'LYS',
                     'DNP': 'ALA', 'DPN': 'PHE', 'DPR': 'PRO', 'DSN': 'SER',
                     'DSP': 'ASP', 'DTH': 'THR', 'DTR': 'TRP', 'DTY': 'TYR',
                     'DVA': 'VAL', 'EFC': 'CYS', 'FLA': 'ALA', 'FME': 'MET',
                     'GGL': 'GLU', 'GL3': 'GLY', 'GLZ': 'GLY', 'GMA': 'GLU',
                     'GSC': 'GLY', 'HAC': 'ALA', 'HAR': 'ARG', 'HIC': 'HIS',
                     'HIP': 'HIS', 'HMR': 'ARG', 'HPQ': 'PHE', 'HTR': 'TRP',
                     'HYP': 'PRO', 'IAS': 'ASP', 'IIL': 'ILE', 'IYR': 'TYR',
                     'KCX': 'LYS', 'LLP': 'LYS', 'LLY': 'LYS', 'LTR': 'TRP',
                     'LYM': 'LYS', 'LYZ': 'LYS', 'MAA': 'ALA', 'MEN': 'ASN',
                     'MHS': 'HIS', 'MIS': 'SER', 'MLE': 'LEU', 'MPQ': 'GLY',
                     'MSA': 'GLY', 'MSE': 'MET', 'MVA': 'VAL', 'NEM': 'HIS',
                     'NEP': 'HIS', 'NLE': 'LEU', 'NLN': 'LEU', 'NLP': 'LEU',
                     'NMC': 'GLY', 'OAS': 'SER', 'OCS': 'CYS', 'OMT': 'MET',
                     'PAQ': 'TYR', 'PCA': 'GLU', 'PEC': 'CYS', 'PHI': 'PHE',
                     'PHL': 'PHE', 'PR3': 'CYS', 'PRR': 'ALA', 'PTR': 'TYR',
                     'PYX': 'CYS', 'SAC': 'SER', 'SAR': 'GLY', 'SCH': 'CYS',
                     'SCS': 'CYS', 'SCY': 'CYS', 'SEL': 'SER', 'SEP': 'SER',
                     'SET': 'SER', 'SHC': 'CYS', 'SHR': 'LYS', 'SMC': 'CYS',
                     'SOC': 'CYS', 'STY': 'TYR', 'SVA': 'SER', 'TIH': 'ALA',
                     'TPL': 'TRP', 'TPO': 'THR', 'TPQ': 'ALA', 'TRG': 'LYS',
                     'TRO': 'TRP', 'TYB': 'TYR', 'TYI': 'TYR', 'TYQ': 'TYR',
                     'TYS': 'TYR', 'TYY': 'TYR'}
    substitutions_list = substitutions.keys()
    atom_conect_main_dict = {' CA ': [[' N  ', 1]],
                             ' C  ': [[' CA ', 1]],
                             ' O  ': [[' C  ', 2]],
                             ' OXT': [[' C  ', 1]],
                             }
    atom_conect_heavy_dict = {'ARG': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD ': [[' CG ', 1]],
                                      ' NE ': [[' CD ', 1]],
                                      ' CZ ': [[' NE ', 1]],
                                      ' NH1': [[' CZ ', 2]],
                                      ' NH2': [[' CZ ', 1]], },
                              'HIS': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' ND1': [[' CG ', 1]],
                                      ' CD2': [[' CG ', 2]],
                                      ' NE2': [[' CD2', 1]],
                                      ' CE1': [[' ND1', 1.5],
                                               [' NE2', 1.5]], },
                              'LYS': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD ': [[' CG ', 1]],
                                      ' CE ': [[' CD ', 1]],
                                      ' NZ ': [[' CE ', 1]], },
                              'ASP': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' OD1': [[' CG ', 2]],
                                      ' OD2': [[' CG ', 1]], },
                              'GLU': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD ': [[' CG ', 1]],
                                      ' OE1': [[' CD ', 2]],
                                      ' OE2': [[' CD ', 1]], },
                              'SER': {' CB ': [[' CA ', 1]],
                                      ' OG ': [[' CB ', 1]], },
                              'THR': {' CB ': [[' CA ', 1]],
                                      ' OG1': [[' CB ', 1]],
                                      ' CG2': [[' CB ', 1]], },
                              'ASN': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' OD1': [[' CG ', 2]],
                                      ' ND2': [[' CG ', 1]], },
                              'GLN': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD ': [[' CG ', 1]],
                                      ' OE1': [[' CD ', 2]],
                                      ' NE2': [[' CD ', 1]], },
                              'CYS': {' CB ': [[' CA ', 1]],
                                      ' SG ': [[' CB ', 1]], },
                              'GLY': {},
                              'PRO': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD ': [[' CG ', 1], [' N  ', 1]], },
                              'ALA': {' CB ': [[' CA ', 1]], },
                              'VAL': {' CB ': [[' CA ', 1]],
                                      ' CG1': [[' CB ', 1]],
                                      ' CG2': [[' CB ', 1]], },
                              'ILE': {' CB ': [[' CA ', 1]],
                                      ' CG1': [[' CB ', 1]],
                                      ' CG2': [[' CB ', 1]],
                                      ' CD1': [[' CG1', 1]], },
                              'LEU': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD1': [[' CG ', 1]],
                                      ' CD2': [[' CG ', 1]], },
                              'MET': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' SD ': [[' CG ', 1]],
                                      ' CE ': [[' SD ', 1]], },
                              'PHE': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD1': [[' CG ', 2]],
                                      ' CD2': [[' CG ', 1]],
                                      ' CE1': [[' CD1', 1]],
                                      ' CE2': [[' CD2', 2]],
                                      ' CZ ': [[' CE1', 2], [' CE2', 1]], },
                              'TYR': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD1': [[' CG ', 2]],
                                      ' CD2': [[' CG ', 1]],
                                      ' CE1': [[' CD1', 1]],
                                      ' CE2': [[' CD2', 2]],
                                      ' CZ ': [[' CE1', 2], [' CE2', 1]],
                                      ' OH ': [[' CZ ', 1]]},
                              'TRP': {' CB ': [[' CA ', 1]],
                                      ' CG ': [[' CB ', 1]],
                                      ' CD1': [[' CG ', 2]],
                                      ' CD2': [[' CG ', 1]],
                                      ' NE1': [[' CD1', 1]],
                                      ' CE2': [[' CD2', 2], [' NE1', 1]],
                                      ' CE3': [[' CD2', 1]],
                                      ' CZ2': [[' CE2', 1]],
                                      ' CZ3': [[' CE3', 2]],
                                      ' CH2': [[' CZ2', 2], [' CZ3', 1]], },
                              }

    @classmethod
    def read_pdb_protein(cls, input_file, remain_remark=True,
                         exclude_water=False, ligand_split=True):
        '''
            read pdb file which has protein and ligand atoms and metal atoms.
            residue_re-numbering for pdb fixer
            pdbfixer does not identify 4, 4a...
            if residue_index_initialize == 0
                no new numbering
            if == 1
                new numbering (begin index is original begin index)
            if == 2
                new numbering begin index is 1
                example A:4, A:4a, A:5, ...
                    -> A:1, A:2, A:3, ...

        '''

        fp = open(input_file)
        lines = fp.readlines()
        fp.close()

        protein_chain_residues = dict()

        pdb_info_lines = list()
#        ligand_info_lines = list()
        protein_dict = dict()
        ligand_dict = dict()
        conect_dict = dict()

        for line in lines:
            if line[0:6] == 'HEADER':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'TITLE':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'DBREF ':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'SEQRES':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'REMARK' and remain_remark:
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'COMPND':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'AUTHOR':
                pdb_info_lines += [line]
                continue

            if line[0:6] == 'ATOM  ':
                residue_name = line[17:20].strip()
                residue_num = int(line[22:26])
                residue_num2 = line[22:27]

                chain_id = line[21]
                altLoc = line[16]
                if altLoc != ' ' and altLoc != 'A':
                    continue
                atom_number = int(line[6:11])

                line_out = line

                protein_dict[atom_number] = line_out
                if chain_id not in protein_chain_residues:
                    protein_chain_residues[chain_id] = dict()
                if residue_num2 not in protein_chain_residues[chain_id]:
                    protein_chain_residues[chain_id][residue_num2] = residue_num2

        for line in lines:
            if line[0:6] == 'HETATM':
                residue_name = line[17:20].strip()
                if residue_name == 'HOH' and exclude_water:
                    continue
                residue_num = int(line[22:26])
                residue_num2 = line[22:27]

                chain_id = line[21]
                atom_number = int(line[6:11])
                altLoc = line[16]
                if altLoc != ' ' and altLoc != 'A':
                    continue
                if altLoc == 'A':
                    line = '%s %s' % (line[:16], line[17:])

                if chain_id not in protein_chain_residues:
                    protein_chain_residues[chain_id] = dict()

                if residue_num2 in protein_chain_residues[chain_id]:
                    protein_dict[atom_number] = line
#                    protein_dict[atom_number] = 'ATOM  ' + line[6:]

                    continue

                if residue_name in cls.substitutions_list:
                    if residue_num2 not in protein_chain_residues[chain_id]:
                        protein_chain_residues[chain_id][residue_num2] = residue_num
                    protein_dict[atom_number] = line
#                    protein_dict[atom_number] = 'ATOM  ' + line[6:]
                    continue

#                residue_name2 = residue_name
                atom_name = line[12:16]
#                if residue_name[-1].isdigit():
#                    residue_name2 = residue_name[:-1]
#                if residue_name2 in cls.metal_atom_list + ['HEM'] :
                if atom_name in cls.metal_atom_list:
                    if residue_num2 not in protein_chain_residues[chain_id]:
                        protein_chain_residues[chain_id][residue_num2] = residue_num
                    protein_dict[atom_number] = line
                elif not ligand_split:
                    if residue_num2 not in protein_chain_residues[chain_id]:
                        protein_chain_residues[chain_id][residue_num2] = residue_num
                    protein_dict[atom_number] = line
                else:
                    ligand_dict[atom_number] = line

        for line in lines:
            if line[0:6] == 'CONECT':
                conect_list = []
                for i in range(0, 8):
                    ini = i * 5 + 6
                    fin = (i + 1) * 5 + 6
                    atom_number = line[ini:fin].strip()
                    if len(atom_number) > 0:
                        conect_list += [int(atom_number)]
                conect_idx = conect_list[0]
                if conect_idx not in conect_dict:
                    conect_dict[conect_idx] = conect_list[1:]
                else:
                    conect_dict[conect_idx] = conect_dict[conect_idx] + \
                        conect_list[1:]

        return (protein_chain_residues, pdb_info_lines, protein_dict,
                ligand_dict, conect_dict)

    @classmethod
    def initialize_residue_index(self, result, initialize_option=1):

        (protein_chain_residues, pdb_info_lines, protein_dict,
            ligand_dict, conect_dict) = result
        protein_dict_new = dict()
        ligand_dict_new = dict()
        protein_chain_residues_new = dict()
        protein_atom_number_list = sorted(protein_dict.keys())

        residue_num_old = ''
        residue_num2_old = ''
        chain_id_old = ''

        for atom_number in protein_atom_number_list:
            line = protein_dict[atom_number]

#            residue_name = line[17:20].strip()
            residue_num = int(line[22:26])
            residue_num2 = line[22:27]
#            insertion = line[26]

            chain_id = line[21]

            if chain_id not in protein_chain_residues_new:
                protein_chain_residues_new[chain_id] = dict()

            if chain_id_old != chain_id:
                chain_id_old = chain_id
                if initialize_option == 1:
                    residue_num_new = residue_num
                if initialize_option == 2:
                    residue_num_new = 1

                residue_num2_old = residue_num2
                residue_num_old = residue_num

            if residue_num2_old != residue_num2:
                residue_num_new += max(residue_num - residue_num_old, 1)
                residue_num2_old = residue_num2
                residue_num_old = residue_num

            if residue_num2 not in protein_chain_residues_new[chain_id]:
                protein_chain_residues_new[chain_id][residue_num2] = residue_num2
            line_out = '%s%4d %s' % (line[:22], residue_num_new, line[27:])
            protein_dict_new[atom_number] = line_out
        ligand_atom_number_list = sorted(ligand_dict.keys())

        for atom_number in ligand_atom_number_list:
            line = ligand_dict[atom_number]

#            residue_name = line[17:20].strip()
            residue_num = int(line[22:26])
            residue_num2 = line[22:27]
#            insertion = line[26]

            chain_id = line[21]
            if chain_id not in protein_chain_residues_new:
                protein_chain_residues_new[chain_id] = dict()

            if chain_id_old != chain_id:
                chain_id_old = chain_id
                residue_num2_old_list = protein_chain_residues_new[chain_id].keys(
                )
                if len(residue_num2_old_list) > 0:
                    residue_num2_old = max(
                        protein_chain_residues_new[chain_id].keys())
                    residue_num_new = protein_chain_residues_new[chain_id][residue_num2_old]
                else:
                    residue_num_new = 0
                residue_num_old = residue_num

            if residue_num2_old != residue_num2:
                residue_num_new += max(residue_num - residue_num_old, 1)

                residue_num2_old = residue_num2
                residue_num_old = residue_num

            line_out = '%s%4d %s' % (line[:22], residue_num_new, line[27:])
            ligand_dict_new[atom_number] = line_out
            if residue_num2 not in protein_chain_residues[chain_id]:
                protein_chain_residues[chain_id][residue_num2] = residue_num_new

        result_new = (protein_chain_residues_new, pdb_info_lines,
                      protein_dict_new, ligand_dict_new, conect_dict)
        return result_new

    @classmethod
    def read_pdb_ligand(cls, input_file):

        fp = open(input_file)
        lines = fp.readlines()
        fp.close()
        model_dict = dict()
        pdb_info_lines = list()
        ligand_dict = dict()
        conect_dict = dict()
        for line in lines:
            if line[0:6] == 'MODEL ':
                model_id = int(line[6:].strip())
                pdb_info_lines = list()
                ligand_dict = dict()
                conect_dict = dict()

            if line[0:6] == 'HEADER':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'TITLE':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'REMARK':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'COMPND':
                pdb_info_lines += [line]
                continue
            if line[0:6] == 'AUTHOR':
                pdb_info_lines += [line]
                continue

            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                altLoc = line[16]
                if altLoc != ' ' and altLoc != 'A':
                    continue
                atom_number = int(line[6:11])
                line_out = line
                ligand_dict[atom_number] = line_out

            if line[0:6] == 'CONECT':
                conect_list = []
                for i in range(0, 8):
                    ini = i * 5 + 6
                    fin = (i + 1) * 5 + 6
                    atom_number = line[ini:fin].strip()
                    if len(atom_number) > 0:
                        conect_list += [int(atom_number)]
                conect_idx = conect_list[0]
                if conect_idx not in conect_dict:
                    conect_dict[conect_idx] = conect_list[1:]
                else:
                    conect_dict[conect_idx] = conect_dict[conect_idx] + \
                        conect_list[1:]
            if line[0:6] == 'ENDMDL':
                model_dict[model_id] = (
                    pdb_info_lines, ligand_dict, conect_dict)
        if len(model_dict.keys()) == 0:
            model_dict[1] = (pdb_info_lines, ligand_dict, conect_dict)

        return model_dict

    @classmethod
    def write_one_model(cls, molecule_model):
        (pdb_info_lines, molecule_dict, conect_dict) = molecule_model
        line_out_total = ''
        for line in pdb_info_lines:
            line_out = line
            line_out_total += line_out

        molecule_atom_numbers = sorted(molecule_dict.keys())

        enu = list(enumerate(molecule_atom_numbers))
        num_atom = len(enu)
        for i, atom_number in enu:
            line_out = molecule_dict[atom_number]
            line_out_total += line_out
            chain_id = molecule_dict[atom_number][21]
            if molecule_dict[atom_number][0:6] == 'ATOM  ':
                j = i + 1
                if j >= num_atom:
                    line_out = 'TER\n'
                    line_out_total += line_out
                elif enu[j][1] not in molecule_dict:
                    line_out = 'TER\n'
                    line_out_total += line_out
                elif molecule_dict[enu[j][1]][0:6] == 'HETATM':
                    line_out = 'TER\n'
                    line_out_total += line_out
                elif (molecule_dict[enu[j][1]][0:6] == 'ATOM  ') and (molecule_dict[enu[j][1]][21] != chain_id):
                    line_out = 'TER\n'
                    line_out_total += line_out

        for atom_number in molecule_atom_numbers:
            if atom_number not in conect_dict:
                continue
            ans = conect_dict[atom_number]
            ans2 = list()
            for an in ans:
                if an in molecule_atom_numbers:
                    ans2 += [an]
            num_conect = len(ans2)
            line_out = ''
            for i_con in range(num_conect):
                if i_con % 4 == 0:
                    line_out += 'CONECT%5d' % (atom_number)
                line_out += '%5d' % (ans2[i_con])
                if i_con % 4 == 3:
                    line_out += '\n'
            if len(line_out.strip()) < 1:
                continue
            if line_out[-1] != '\n':
                line_out += '\n'
            line_out_total += line_out
        return line_out_total

    @classmethod
    def write_model_pdb(cls, model_dict, output_file):
        fp_p = open(output_file, 'w')
        keys = sorted(model_dict.keys())
        num_model = len(keys)
        for model_id in keys:
            molecule_model = model_dict[model_id]
            if num_model > 1:
                line_out = 'MODEL %8d\n' % model_id
                fp_p.write(line_out)
            line_out = cls.write_one_model(molecule_model)
            fp_p.write(line_out)
            if num_model > 1:
                line_out = 'ENDMDL\n'
                fp_p.write(line_out)
        line_out = 'END\n'
        fp_p.write(line_out)
        fp_p.close()

    @classmethod
    def obabel_rewrite(cls, input_file, output_file, option=None):
        #        ms = pybel.readfile(mol_format, input_file)
        #        m = list(ms)[0]
        #        m.write('pdb', output_file, overwrite=True)
        run_line = 'obabel %s -O %s' % (input_file, output_file)
        if option is not None:
            run_line += ' %s' % (option)
        subprocess.check_output(run_line.split(),
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)

    @classmethod
    def split_chain(cls, input_file, out_dir, pdb_code, exclude_water=False):
        result = cls.read_pdb_protein(input_file, remain_remark=True,
                                      exclude_water=exclude_water)
        (protein_chain_residues, pdb_info_lines, protein_dict,
         ligand_dict, conect_dict) = result
        chain_list = protein_chain_residues.keys()
        chain_dict = dict()
        for chain_id in chain_list:
            chain_dict[chain_id] = dict()
            chain_dict[chain_id]['pdb_info_lines'] = list()
            chain_dict[chain_id]['protein_dict'] = dict()
            chain_dict[chain_id]['ligand_dict'] = dict()
            chain_dict[chain_id]['mol_dict'] = dict()
            chain_dict[chain_id]['conect_dict'] = dict()

        for line in pdb_info_lines:
            if line[0:6] != 'SEQRES':
                for chain_id in chain_list:
                    chain_dict[chain_id]['pdb_info_lines'] += [line]
                continue
            chain_id = line[11]
            chain_dict[chain_id]['pdb_info_lines'] += [line]

        for atom_number in protein_dict:
            line = protein_dict[atom_number]
            chain_id = line[21]
            chain_dict[chain_id]['mol_dict'][atom_number] = line

        for atom_number in ligand_dict:
            line = ligand_dict[atom_number]
            chain_id = line[21]
            chain_dict[chain_id]['mol_dict'][atom_number] = line

        for chain_id in chain_list:
            chain = chain_dict[chain_id]['mol_dict']
            for atom_number in chain.keys():
                if atom_number not in conect_dict:
                    continue
                line = conect_dict[atom_number]
                chain_dict[chain_id]['conect_dict'][atom_number] = line

        for chain_id in chain_list:
            output_file = '%s/%s%s.pdb' % (out_dir, pdb_code, chain_id)
            model_dict = dict()
            model_dict[0] = (chain_dict[chain_id]['pdb_info_lines'],
                             chain_dict[chain_id]['mol_dict'],
                             chain_dict[chain_id]['conect_dict'])
#            print(model_dict[0][1])
            cls.write_model_pdb(model_dict, output_file)

    @classmethod
    def rot_pdb(cls, molecule_dict, U, t):

        molecule_dict_new = dict()
        molecule_atom_numbers = sorted(molecule_dict.keys())

        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
            coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            coor = np.array(coor)
            coor1 = np.dot(coor, U)
            coor2 = coor1 + t
            line_out = "%s%8.3f%8.3f%8.3f%s" % (
                line[0:30], coor2[0], coor2[1], coor2[2], line[54:])
            molecule_dict_new[atom_number] = line_out
        return molecule_dict_new

    @classmethod
    def align_3d(cls, input_file, protein_dict, ligand_dict, ref_file, tmp_dir):

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        matrix_file = tmp_dir+"/matrix.txt"
        run_line = "TMalign %s %s -m %s" % (input_file, ref_file, matrix_file)
        subprocess.check_output(run_line.split(),
                                stderr=subprocess.STDOUT,
                                timeout=10, universal_newlines=True)
        U = np.zeros([3, 3])
        t = np.zeros([3])
        fp = open(matrix_file)
        line = fp.readline()
        line = fp.readline()
        for i in range(3):
            line = fp.readline()
            lis = line.strip().split()
            t[i] = float(lis[1])
            for j in range(3):
                U[j][i] = float(lis[j+2])

        protein_dict_new = cls.rot_pdb(protein_dict, U, t)
        ligand_dict_new = cls.rot_pdb(ligand_dict, U, t)

        return protein_dict_new, ligand_dict_new

    @classmethod
    def split_ligand(cls, ligand_dict, conect_dict):
        ligand_mol_dict = dict()
        conect_mol_dict = dict()
        atom_number_list = ligand_dict.keys()
        for atom_number in atom_number_list:
            line = ligand_dict[atom_number]
            residue_name = line[17:20].strip()
#            residue_num = int(line[22:26])
            if residue_name not in ligand_mol_dict:
                ligand_mol_dict[residue_name] = dict()
                conect_mol_dict[residue_name] = dict()
            ligand_mol_dict[residue_name][atom_number] = line
            if atom_number not in conect_dict:
                continue
            conect_mol_dict[residue_name][atom_number] = conect_dict[atom_number]

        return ligand_mol_dict, conect_mol_dict

    @classmethod
    def fix_pdb(cls, protein_file, protein_fix_file, add_missing_residue=True, pH=7.4):
        fixer = pdbfixer.PDBFixer(filename=protein_file)

#        fixer.applyMutations(["ALA-57-GLY"], "A")
#        fixer.missingResidues = {}
        if add_missing_residue:
            a = fixer.findMissingResidues()
        else:
            fixer.missingResidues = {}
        missing_residues = fixer.missingResidues

        fixer.findNonstandardResidues()
        nonstandard_residues = fixer.nonstandardResidues
        fixer.replaceNonstandardResidues()
#        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        missing_atoms = fixer.missingAtoms
        missing_terminals = fixer.missingTerminals

        fixer.addMissingAtoms()
        if pH is not None:
            try:
                fixer.addMissingHydrogens(pH=pH)
            except Exception as e:
                print('Error:', e, protein_file)
        fp_out = open(protein_fix_file, 'w')

        app.PDBFile.writeFile(fixer.topology, fixer.positions, fp_out, True)
        fp_out.close()

        return (missing_residues, nonstandard_residues,
                missing_atoms, missing_terminals)

    @classmethod
    def fix_charge(cls, molecule_dict, conect_dict):
        molecule_dict_new = dict()
        conect_dict_new = copy.deepcopy(conect_dict)
        molecule_atom_numbers = sorted(molecule_dict.keys())
        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
            if atom_number not in conect_dict_new:
                continue
            atom_type = line[77]
            if atom_type == 'N':
                ans = conect_dict_new[atom_number]
                con_dict = dict()
                for an in ans:
                    if an in molecule_atom_numbers:
                        if an not in con_dict:
                            con_dict[an] = 0
                        con_dict[an] += 1
                con_dict2 = dict()
                ans2 = list()
                keys = sorted(con_dict.keys())
                check_o2 = False
                check_mod_o2 = False
                o2_atom_num = 0
                for key in keys:
                    bond_type = con_dict[key]
                    line2 = molecule_dict[key]
                    atom_type2 = line2[77]
                    if atom_type2 == 'O' and bond_type == 2 and not check_o2:
                        check_o2 = True
                    elif atom_type2 == 'O' and bond_type == 2 and check_o2:
                        bond_type = 1
                        o2_atom_num = key
                        check_mod_o2 = True
                    con_dict2[key] = bond_type
                    for i in range(bond_type):
                        ans2 += [key]
                conect_dict_new[atom_number] = ans2
                if check_mod_o2:
                    ans_o = conect_dict_new[o2_atom_num]
                    con_dict = dict()
                    for an in ans_o:
                        if an in molecule_atom_numbers:
                            if an not in con_dict:
                                con_dict[an] = 0
                            con_dict[an] += 1
                    ans_o2 = list()
                    keys = sorted(con_dict.keys())
                    for key in keys:
                        bond_type = con_dict[key]
                        if key == atom_number and bond_type == 2:
                            bond_type = 1
                        for i in range(bond_type):
                            ans_o2 += [key]
                    conect_dict_new[o2_atom_num] = ans_o2

        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
            atom_name = line[12:16]
#            residue_name2 = residue_name
#            if residue_name[-1].isdigit():
#                residue_name2 = residue_name[:-1]
#             if residue_name2 in cls.metal_atom_list + ['HEM'] :
            if atom_name in cls.metal_atom_list:
                molecule_dict_new[atom_number] = line

            if atom_number not in conect_dict_new:
                continue
            ans = conect_dict_new[atom_number]
            ans2 = list()
            for an in ans:
                if an in molecule_atom_numbers:
                    ans2 += [an]
            num_conect = len(ans2)
            atom_type = line[77]
            charge = 0
            if atom_type == 'N':
                charge = num_conect - 3
                if charge < 0:
                    charge = 0
            if atom_type == 'O':
                charge = num_conect - 2

            if charge > 0:
                line_out = '%s%1d+\n' % (line[0:78], charge)
            elif charge < 0:
                line_out = '%s%1d-\n' % (line[0:78], -charge)
            else:
                line_out = line
            molecule_dict_new[atom_number] = line_out

        return molecule_dict_new, conect_dict_new

    @classmethod
    def fix_pdb_atom_idx(cls, molecule_dict):
        molecule_dict_new = dict()
        molecule_atom_numbers = sorted(molecule_dict.keys())

        atom_dict = dict()
        check_change = False

        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
            atom = line[12:16]
            if atom.strip()[-1].isnumeric():
                line_out = line
            elif atom[2] != ' ':
                line_out = line
#            elif atom[0] == 'H' or atom[1] == 'H':
#                line_out = line
            else:
                at = atom[0:2]
                if at not in atom_dict:
                    atom_dict[at] = 0
                atom_dict[at] += 1
                idx = atom_dict[at]
                if idx >= 100:
                    return 'error: index >=100', False
                line_out = line[0:12] + '%s%-2d' % (at, idx) + line[16:]
                check_change = True

            molecule_dict_new[atom_number] = line_out
        return molecule_dict_new, check_change

    @classmethod
    def correct_ligand_pep(cls, molecule_dict):
        molecule_dict_new = dict()
        molecule_atom_numbers = sorted(molecule_dict.keys())

        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
            residue_name = line[17:20].strip()
#            residue_num = int(line[22:26])
#            residue_num2 = line[22:27]
            atom_name = line[12:16].strip()
            if residue_name != 'UNK':
                line_out = 'HETATM%sUNK%s%4d%s' % (
                    line[6:17], line[20:22], 1, line[26:])
                if atom_name == 'OXT':
                    line_out = '%s O  %s' % (line_out[:12], line_out[16:])
                if atom_name == 'CA':
                    line_out = '%s C  %s' % (line_out[:12], line_out[16:])

            else:
                line_out = line
            molecule_dict_new[atom_number] = line_out
        return molecule_dict_new

    @classmethod
    def gen_conect_protein(cls, molecule_dict, conect_dict):

        chain_dict = dict()
        molecule_atom_numbers = sorted(molecule_dict.keys())
        for atom_number in molecule_atom_numbers:
            line = molecule_dict[atom_number]
#            if line[0:6] != 'ATOM  ':
#                continue
            atom_idx = int(line[6:11])
            atom_name = line[12:16]
            chain = line[21]
            residue_name = line[17:20]
#            residue_idx = int(line[22:26])
            residue_idx2 = line[22:27]
#            atom_type = line[77]
#            atom_type = atom_name.strip()[0]

            charge = line[80:77:-1].strip()
            if charge == '':
                charge = 0
            else:
                charge = int(charge)
            if chain not in chain_dict:
                chain_dict[chain] = dict()
            if residue_idx2 not in chain_dict[chain]:
                chain_dict[chain][residue_idx2] = dict()
            if atom_name not in chain_dict[chain][residue_idx2]:
                chain_dict[chain][residue_idx2][atom_name] = dict()
                chain_dict[chain][residue_idx2][atom_name]['atom_idx'] = atom_idx
                if atom_name in conect_dict:
                    chain_dict[chain][residue_idx2][atom_name]['conect_list'] = copy.deepcopy(
                        conect_dict[atom_name])
                else:
                    chain_dict[chain][residue_idx2][atom_name]['conect_list'] = list()

                chain_dict[chain][residue_idx2][atom_name]['residue_name'] = residue_name

        chain_keys = sorted(chain_dict.keys())
#        print(chain_dict['A'].keys())
#        print(chain_dict['A']['501'])
#        sys.exit()
        for chain in chain_keys:
            for residue_idx2 in chain_dict[chain].keys():
                residue_dict = chain_dict[chain][residue_idx2]
                atom_name_keys = residue_dict.keys()
                for atom_name in atom_name_keys:
                    atom_dict = residue_dict[atom_name]
                    conect_list = atom_dict['conect_list']
                    atom_idx = atom_dict['atom_idx']
                    residue_name = atom_dict['residue_name']
                    if residue_name == 'HOH':
                        if atom_name == ' O  ':
                            continue
                        if atom_name[0:2] == ' H':
                            previous_atom_list = [[' O  ', 1]]
                            for previous_atom in previous_atom_list:
                                previous_atom_name = previous_atom[0]
                                previous_bond_type = previous_atom[1]
                                previous_atom_dict = residue_dict[previous_atom_name]
                                previous_atom_idx = previous_atom_dict['atom_idx']
                                previous_conect_list = previous_atom_dict['conect_list']
                                if previous_atom_idx not in conect_list:
                                    conect_list += [previous_atom_idx] * \
                                        previous_bond_type
                                if atom_idx not in previous_conect_list:
                                    previous_conect_list += [atom_idx] * \
                                        previous_bond_type

                    elif atom_name == ' N  ':
                        residue_idx2_previous = '%4d ' % (
                            int(residue_idx2[0:4])-1)
                        if residue_idx2_previous not in chain_dict[chain]:
                            continue
                        if ' C  ' not in chain_dict[chain][residue_idx2_previous]:
                            continue
                        previous_residue_c = chain_dict[chain][residue_idx2_previous][' C  ']
                        atom_idx_prc = previous_residue_c['atom_idx']
                        conect_list_prc = previous_residue_c['conect_list']
                        if atom_idx_prc not in conect_list:
                            conect_list += [atom_idx_prc]
                        if atom_idx not in conect_list_prc:
                            conect_list_prc += [atom_idx]

                    elif atom_name.strip()[0] != 'H':
                        check_atom = False
                        previous_atom_list = list()
                        if atom_name in cls.atom_conect_main_dict:
                            previous_atom_list = cls.atom_conect_main_dict[atom_name]
                            check_atom = True
                        elif residue_name in cls.atom_conect_heavy_dict:
                            atom_conect_dict = cls.atom_conect_heavy_dict[residue_name]
                            if atom_name in atom_conect_dict:
                                previous_atom_list = atom_conect_dict[atom_name]
                                check_atom = True
                            else:
                                print('atom_name', atom_name, 'is not in dict')
                        elif atom_name in cls.metal_atom_list:
                            check_atom = False
                        else:
                            print('residue_name', residue_name, 'is not in dict')
                        if not check_atom:
                            continue

                        for previous_atom in previous_atom_list:
                            previous_atom_name = previous_atom[0]
                            previous_bond_type = previous_atom[1]
                            previous_atom_dict = residue_dict[previous_atom_name]
                            previous_atom_idx = previous_atom_dict['atom_idx']
                            previous_conect_list = previous_atom_dict['conect_list']
                            if 1 < previous_bond_type < 2:
                                #                            if previous_bond_type == 1.5:
                                if previous_atom_name[0:2].strip() != 'N':
                                    print('strange:', chain, residue_idx2,
                                          atom_name, previous_atom_name)
                                    continue
                                ph_atom_name = ' H%s' % previous_atom_name[2:4]
                                if ph_atom_name in atom_name_keys:
                                    previous_bond_type = 1
                                else:
                                    previous_bond_type = 2

                            if previous_atom_idx not in conect_list:
                                conect_list += [previous_atom_idx] * \
                                    previous_bond_type
                            if atom_idx not in previous_conect_list:
                                previous_conect_list += [atom_idx] * \
                                    previous_bond_type

                    else:
                        check_atom = False
                        atom_name2 = atom_name.strip()[1:]
                        lenan = len(atom_name2)
                        previous_atom_list = list()
                        if lenan == 0:
                            previous_atom_list = [[' N  ', 1]]
                        else:
                            if atom_name2[0].isdecimal():
                                previous_atom_list = [[' N  ', 1]]
                            else:
                                for atom_name_pp in atom_name_keys:
                                    if atom_name_pp.strip()[0] == 'H':
                                        continue
                                    atom_name_pp2 = atom_name_pp[2:4].strip()
                                    if atom_name2 == atom_name_pp2:
                                        previous_atom_list = [
                                            [atom_name_pp, 1]]
                                    elif lenan >= 2:
                                        if atom_name2[: -1] == atom_name_pp2:
                                            previous_atom_list = [
                                                [atom_name_pp, 1]]
                        if len(previous_atom_list) > 0:
                            for previous_atom in previous_atom_list:
                                previous_atom_name = previous_atom[0]
                                previous_bond_type = previous_atom[1]
                                previous_atom_dict = residue_dict[previous_atom_name]
                                previous_atom_idx = previous_atom_dict['atom_idx']
                                previous_conect_list = previous_atom_dict['conect_list']
                                if previous_atom_idx not in conect_list:
                                    conect_list += [previous_atom_idx] * \
                                        previous_bond_type
                                if atom_idx not in previous_conect_list:
                                    previous_conect_list += [atom_idx] * \
                                        previous_bond_type
                        else:
                            print('conect search fail:', chain,
                                  residue_idx2, atom_name)

        conect_dict_new = dict()
        for chain in chain_dict.keys():
            for residue_idx2 in chain_dict[chain].keys():
                residue_dict = chain_dict[chain][residue_idx2]
                atom_name_keys = residue_dict.keys()
                for atom_name in atom_name_keys:
                    atom_dict = residue_dict[atom_name]
                    conect_list = atom_dict['conect_list']
                    atom_idx = atom_dict['atom_idx']
                    conect_dict_new[atom_idx] = sorted(conect_list)
        return conect_dict_new

    @classmethod
    def fix_protein(cls, input_file, output_file, add_missing_residue=True, pH=7.4):
        result = cls.read_pdb_protein(input_file, remain_remark=True,
                                      ligand_split=False)
#        (protein_chain_residues, pdb_info_lines, protein_dict,
#            ligand_dict, conect_dict) = result
        result_new = cls.initialize_residue_index(result, initialize_option=1)
        (protein_chain_residues, pdb_info_lines, protein_dict,
            ligand_dict, conect_dict) = result_new

        model_dict = dict()
        model_dict[0] = (pdb_info_lines, protein_dict, conect_dict)
        cls.write_model_pdb(model_dict, output_file)
#        cls.obabel_rewrite(output_file, output_file, option=None)
        cls.fix_pdb(output_file, output_file, add_missing_residue=True, pH=7.4)
#        cls.obabel_rewrite(output_file, output_file)

        result = cls.read_pdb_protein(
            output_file, remain_remark=True, ligand_split=False)
        (protein_chain_residues, pdb_info_lines, protein_dict,
         ligand_dict, conect_dict) = result
        conect_dict_new = cls.gen_conect_protein(protein_dict, conect_dict)
#        model_dict = dict()
#        model_dict[0] = (pdb_info_lines, protein_dict, conect_dict_new)
#        cls.write_model_pdb(model_dict, output_file)

        result_new = cls.fix_charge(protein_dict, conect_dict_new)
        (protein_dict_new, conect_dict_new) = result_new
        model_dict = dict()
        model_dict[0] = (pdb_info_lines, protein_dict_new, conect_dict_new)
        cls.write_model_pdb(model_dict, output_file)
#        cls.obabel_rewrite(output_file, output_file)

    @ classmethod
    def read_cif(cls, input_file):
        bond_name_dict = {'SING': 1, 'DOUB': 2, 'TRIP': 3}
        fp = open(input_file)
        lines = fp.readlines()
        fp.close()

        t_line = False
        is_atom = False
        is_bond = False
#        mol_title = list()
#        atom_title = list()
#        bond_title = list()
#        atom_list = list()
        bond_list = list()
        for line in lines:
            if line.startswith('data'):
                continue
            if line.startswith('#'):
                is_atom = False
                is_bond = False
            if line.startswith('loop_'):
                t_line = True
                continue
            if t_line:
                if line.startswith('_'):
                    if line.startswith('_chem_comp_atom'):
                        is_atom = True
#                        atom_title += [line.strip().split('.')[1]]
                        continue
                    if line.startswith('_chem_comp_bond'):
                        is_bond = True
#                        bond_title += [line.strip().split('.')[1]]
                        continue
                else:
                    t_line = False
            if is_atom:
                #                lis = line.strip().split()
                #                atom_name = lis[1].strip('"')
                continue
            if is_bond:
                lis = line.strip().split()
                atom1 = lis[1].strip('"')
                atom2 = lis[2].strip('"')
                bond_type = bond_name_dict[lis[3]]
                aromatic = False
                if lis[4] == 'Y':
                    aromatic = True
                bond_list += [(atom1, atom2, bond_type, aromatic)]
                continue
        return bond_list

    @ classmethod
    def fix_conect(cls, ligand_dict, conect_dict, ref_bond_list):
        conect_dict_new = dict()
        keys = ligand_dict.keys()
        atom_num_name_dict = dict()
        for key in keys:
            line = ligand_dict[key]
            atom_number = int(line[6:11])
            atom_name = line[12:16].strip()
            atom_num_name_dict[atom_name] = atom_number

        for ref_bond in ref_bond_list:
            atom1_name = ref_bond[0]
            atom2_name = ref_bond[1]
            bond_type = ref_bond[2]
            if atom1_name in atom_num_name_dict:
                atom1_num = atom_num_name_dict[atom1_name]
            else:
                continue
            if atom2_name in atom_num_name_dict:
                atom2_num = atom_num_name_dict[atom2_name]
            else:
                continue
            if atom1_num in conect_dict:
                clist = conect_dict[atom1_num]
                if atom2_num not in clist:
                    print('strange:', atom2_num, 'not in conect_dict')
                    continue
            else:
                print('strange:', atom1_num, 'not in conect_dict')
                continue
            if atom2_num in conect_dict:
                clist = conect_dict[atom2_num]
                if atom1_num not in clist:
                    print('strange:', atom1_num, 'not in conect_dict')
                    continue
            else:
                print('strange:', atom2_num, 'not in conect_dict')
                continue
            if atom1_num not in conect_dict_new:
                conect_dict_new[atom1_num] = list()
            if atom2_num not in conect_dict_new:
                conect_dict_new[atom2_num] = list()
            conect_dict_new[atom1_num] += [atom2_num] * bond_type
            conect_dict_new[atom2_num] += [atom1_num] * bond_type

        return conect_dict_new

    @ classmethod
    def fix_ligand_ref(cls, input_file, ref_file, output_file, neutralize=False, pH=None, is_peptide=False, is_fix_atom_idx=False):
        ref_bond_list = cls.read_cif(ref_file)
        model_dict = cls.read_pdb_ligand(input_file)
        model_dict_new = dict()
        keys = model_dict.keys()
        for model_id in keys:
            (pdb_info_lines, ligand_dict, conect_dict) = model_dict[model_id]
            conect_dict_new = cls.fix_conect(
                ligand_dict, conect_dict, ref_bond_list)
            model_dict_new[model_id] = (
                pdb_info_lines, ligand_dict, conect_dict_new)
        cls.write_model_pdb(model_dict_new, output_file)

        option_h = '-h'
        cls.obabel_rewrite(output_file, output_file, option=option_h)

        if neutralize:
            option1 = '--neutralize'
            cls.obabel_rewrite(output_file, output_file, option=option1)
        if pH is not None:
            option2 = '-p %.1f' % (pH)
            cls.obabel_rewrite(output_file, output_file, option=option2)
            if not is_peptide:
                model_dict = cls.read_pdb_ligand(output_file)
                model_dict_new = dict()
                keys = model_dict.keys()
                for model_id in keys:
                    (pdb_info_lines, ligand_dict,
                     conect_dict) = model_dict[model_id]
                    ligand_dict_new = cls.correct_ligand_pep(ligand_dict)
                    model_dict_new[model_id] = (
                        pdb_info_lines, ligand_dict_new, conect_dict)
                cls.write_model_pdb(model_dict_new, output_file)

        if is_fix_atom_idx and not is_peptide:
            model_dict = cls.read_pdb_ligand(output_file)
            model_dict_new = dict()
            keys = model_dict.keys()
            check_change = False
            for model_id in keys:
                (pdb_info_lines, ligand_dict,
                 conect_dict) = model_dict[model_id]
                ligand_dict_new, check_change0 = cls.fix_pdb_atom_idx(
                    ligand_dict)
                model_dict_new[model_id] = (
                    pdb_info_lines, ligand_dict_new, conect_dict)
                if check_change0:
                    check_change = True
            if check_change:
                cls.write_model_pdb(model_dict_new, output_file)

    @ classmethod
    def ligand_to_pdbqt(cls, pdb_file, pdbqt_file):
        run_line = 'prepare_ligand4.py -l %s -o %s' % (pdb_file, pdbqt_file)
        run_line += ' -U nphs_lps'
        e = None
        try:
            subprocess.check_output(run_line.split(),
                                    stderr=subprocess.STDOUT,
                                    universal_newlines=True)
        except Exception as e:
            return e
        return e

    @ classmethod
    def protein_to_pdbqt(cls, pdb_file, pdbqt_file):
        run_line = 'prepare_receptor4.py -r %s -o %s' % (pdb_file, pdbqt_file)
#        run_line += ' -U nphs_lps_waters_nonstdres'
        run_line += ' -U nphs_lps_nonstdres'

        e = None
        try:
            subprocess.check_output(run_line.split(),
                                    stderr=subprocess.STDOUT,
                                    universal_newlines=True)
        except Exception as e:
            return e
        return e


def main():

    import argparse
    title_line = 'Fixer for protein pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdb', help='output protein pdb file')

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

    out_format = output_file.split('.')[-1]
    if out_format == 'pdbqt':
        tmp_file = '.'.join(output_file.split('.')[: -1]) + '.pdb'
    else:
        tmp_file = output_file
    PDBtools.fix_protein(input_file, tmp_file)
    if out_format == 'pdbqt':
        e = PDBtools.protein_to_pdbqt(tmp_file, output_file)
        if e is not None:
            print(e)


if __name__ == "__main__":
    main()
