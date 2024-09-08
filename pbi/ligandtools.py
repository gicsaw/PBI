#!/usr/bin/env python
import numpy as np
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess


def ligand_preparation(smi, neutralize, pH):
    """
        input smi
            -> neutralize
            -> protonate
        output protonated smiles
    """
    neutralize = neutralize
    pH = pH
    if neutralize:
        run_line = 'obabel -:%s -osmi --neutralize' % (smi)
        result = subprocess.check_output(run_line.split(),
                                         stderr=subprocess.STDOUT,
                                         universal_newlines=True)
        for line in result.split('\n'):
            if line.startswith('1 molecule converted'):
                continue
            if line.startswith('0 molecule converted'):
                smi0 = smi
                break
            if len(line.strip()) < 1:
                continue
            smi0 = line.strip()
    else:
        smi0 = smi
    if pH is None:
        smi_p = smi0
    else:
        try:
            m = pybel.readstring("smi", smi0)
            m.OBMol.AddHydrogens(False, True, pH)
            smi_p = m.write("smi").strip()
        except Exception:
            smi_p = smi0
    return smi_p


def check_gen3d_rd(m3):
    conformer = m3.GetConformer()
    positions = conformer.GetPositions()
    count_0 = 0
    for i in range(0, positions.shap[0]):
        coor = positions[i]
        x, y, z = coor
        if x == 0.0 and y == 0.0 and z == 0.0:
            count_0 += 1
        if count_0 >= 2:
            check_error = True
            break

    return check_error


def fix_ligand_atom_idx(line_list):

    atom_dict = dict()
    total_line_out = str()

    for line in line_list:
        line = line.rstrip('\n')
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            atom = line[12:16]
            at = atom[0:2]
            chain = line[21]
            line = 'HETATM%s%s   UNK %s   1 %s' % (line[6:12], at,
                                                   chain, line[27:])
        if line[0:6] == 'HETATM':
            atom = line[12:16]
            at = atom[0:2]
            if at not in atom_dict:
                atom_dict[at] = 0
#            chain = line[21]
            atom_dict[at] += 1
            idx = atom_dict[at]
            line_out = line[0:12] + '%s%-2d' % (at, idx) + line[16:]
        else:
            line_out = line
        total_line_out += line_out + '\n'

    return total_line_out


def fix_ligand(input_file, output_file, neutralize=False, pH=None,
               add_hydrogen=True, is_fix_atom_idx=True):
    tmp_file = output_file
    option_h = ''
    if neutralize:
        option_h += ' --neutralize'
    if pH is None and add_hydrogen:
        option_h += ' -h'
    obabel_rewrite(input_file, tmp_file, option=option_h)
#    if option_h != '':
#        obabel_rewrite(input_file, tmp_file, option=option_h)
#    else:
#        tmp_file = input_file

    option_h = ''
    if pH is not None:
        option_h += ' -p %.1f' % (pH)
        obabel_rewrite(tmp_file, output_file, option=option_h)
        tmp_file = output_file

    if is_fix_atom_idx:
        fp = open(tmp_file)
        line_list = fp.readlines()
        fp.close()
        total_line_out = fix_ligand_atom_idx(line_list)
        fp = open(output_file, 'w')
        fp.write(total_line_out)
        fp.close()


def add_mol_id(result, mol_id):

    line_list = result.rstrip('\n').split('\n')
    total_line_out = str()
    check_mol_st = True
    for line in line_list:
        if check_mol_st is True:
            if line.strip() == '':
                line_out = mol_id + '\n'
            else:
                line_out = line + '\n'
            check_mol_st = False
        elif line[0:4] == '$$$$':
            check_mol_st = True
            line_out = line + '\n'
        else:
            line_out = line + '\n'

        total_line_out += line_out

    return total_line_out


def gen_3d(smi, ligand_file, mol_id=None, file_format=None, timeout=10):

    if file_format is None:
        file_format = ligand_file.strip().split('.')[-1]
    if file_format not in ['pdb', 'mol', 'sdf']:
        e = 'error: file_format not in pdb, mol, sdf'
        return e

    m = Chem.MolFromSmiles(smi)
    m3 = Chem.AddHs(m)
    check_error = AllChem.EmbedMolecule(m3)

#    check_error = check_gen3d_rd(m3)
    if check_error:
        e = 'error: gen 3d, two or more (0,0,0)'
        return e

    if file_format == 'pdb':

        #        Chem.MolToPDBFile(m3, ligand_file, flavor=4)
        result_data = Chem.MolToPDBBlock(m3, flavor=4)
        if mol_id is not None:
            total_line_out = result_data.replace('UNL', mol_id)
        else:
            total_line_out = result_data
        fp = open(ligand_file, 'w')
        fp.write(total_line_out)
        fp.close()

    if file_format == 'mol' or file_format == 'sdf':
        result_data = Chem.MolToMolBlock(m3)
        if mol_id is not None:
            total_line_out = add_mol_id(result_data, mol_id)
        else:
            total_line_out = result_data
        fp = open(ligand_file, 'w')
        fp.write(total_line_out)
        fp.close()


def obabel_rewrite(input_file, output_file, option=None):
    #    ms = pybel.readfile(mol_format, input_file)
    #    m = list(ms)[0]
    #    m.write('pdb', output_file, overwrite=True)
    run_line = 'obabel %s -O %s' % (input_file, output_file)
    if option is not None:
        run_line += ' %s' % (option)
    e = subprocess.check_output(run_line.split(),
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)
    return e


def pdb_to_pdbqt(pdb_file, pdbqt_file):

    run_line = 'prepare_ligand4.py -l %s -o %s' % (pdb_file, pdbqt_file)
    run_line += ' -U nphs_lps'
    e = None
    try:
        subprocess.check_output(run_line.split(), stderr=subprocess.STDOUT,
                                universal_newlines=True)
    except Exception as e:
        return e
    return e


def read_ref_pdb_ligand(pdb_file):
    fp = open(pdb_file)
    lines = fp.readlines()
    fp.close()
    atom_dict = dict()
    conect_dict = dict()
    for line in lines:
        if line[0:6] == 'HETATM':
            atom_num = int(line[6:11])
#            atom_name = line[12:16]
            atom_dict[atom_num] = line
        if line[0:6] == 'CONECT':
            conect_list = []
            for i in range(0, 8):
                ini = i * 5 + 6
                fin = (i + 1) * 5 + 6
                atom_num = line[ini:fin].strip()
                if len(atom_num) > 0:
                    conect_list += [int(atom_num)]
            conect_idx = conect_list[0]
            if conect_idx not in conect_dict:
                conect_dict[conect_idx] = conect_list[1:]
            else:
                conect_dict[conect_idx] = conect_dict[conect_idx] + \
                    conect_list[1:]

    return atom_dict, conect_dict


def read_pdbqt_file(pdbqt_file):
    model_dict = dict()
    model_num = 0
    fp = open(pdbqt_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        if line[0:6] == 'MODEL ':
            model_num = int(line[6:].strip())
        if model_num not in model_dict:
            model_dict[model_num] = dict()
            model_dict[model_num]['REMARK'] = list()
            model_dict[model_num]['HETATM'] = dict()

        if line[0:6] == 'REMARK':
            model_dict[model_num]['REMARK'] += [line]
        if line[0:6] == 'HETATM' or line[0:6] == 'ATOM  ':
            atom_name = line[12:16]
            pos = line[30:54]
            model_dict[model_num]['HETATM'][atom_name] = pos

    return model_dict


def write_pdb_one_ref(model, ref_atom_dict, ref_conect_dict):

    total_line_out = ''
    remark_list = model['REMARK']
    for line in remark_list:
        total_line_out += line
    coor_dict = model['HETATM']

    total_atom_list = list()
    keys = ref_atom_dict.keys()
    for atom_num in keys:
        atom_line = ref_atom_dict[atom_num]
        atom_name = atom_line[12:16]
        if atom_name in coor_dict:

            total_atom_list += [atom_num]
            line_out = '%s%s%s' % (
                atom_line[:30], coor_dict[atom_name], atom_line[54:])
            total_line_out += line_out

    keys = ref_conect_dict.keys()
    for atom_num in keys:
        if atom_num not in total_atom_list:
            continue
        ans = ref_conect_dict[atom_num]
        ans2 = list()
        for an in ans:
            if an in total_atom_list:
                ans2 += [an]
        num_conect = len(ans2)
        line_out = ''
        for i_con in range(num_conect):
            if i_con % 4 == 0:
                line_out += 'CONECT%5d' % (atom_num)
            line_out += '%5d' % (ans2[i_con])
            if i_con % 4 == 3:
                line_out += '\n'
        if len(line_out.strip()) < 1:
            continue
        if line_out[-1] != '\n':
            line_out += '\n'
        total_line_out += line_out
    return total_line_out


def pdbqt_to_pdb_ref(input_pdbqt_file, output_pdb_file, ref_pdb_file):
    ref_atom_dict, ref_conect_dict = read_ref_pdb_ligand(ref_pdb_file)
    model_dict = read_pdbqt_file(input_pdbqt_file)
    model_list = model_dict.keys()
    num_model = len(model_list)
    fp_out = open(output_pdb_file, 'w')
    for model_id in model_list:
        total_line_out = write_pdb_one_ref(
            model_dict[model_id], ref_atom_dict, ref_conect_dict)

        if num_model > 1:
            line_out = 'MODEL %8d\n' % model_id
            fp_out.write(line_out)
        fp_out.write(total_line_out)
        if num_model > 1:
            line_out = 'ENDMDL\n'
            fp_out.write(line_out)
    line_out = 'END\n'
    fp_out.write(line_out)
    fp_out.close()

    run_line = 'obabel %s -h -O %s' % (output_pdb_file, output_pdb_file)
    subprocess.check_output(run_line.split(), stderr=subprocess.STDOUT,
                            universal_newlines=True)
    return


def read_coor_pdb(input_file, exclude_Hs=True):
    fp = open(input_file)
    lines = fp.readlines()
    fp.close()
    model_dict = dict()
    ligand_dict = dict()
    for line in lines:
        if line[0: 6] == 'MODEL ':
            model_id = int(line[6:].strip())
            ligand_dict = dict()

        if line[0: 6] == 'ATOM  ' or line[0: 6] == 'HETATM':
            residue_num = int(line[22: 26])
#            residue_name = line[17:20].strip()
#            residue_num2 = line[22:27]
#            atom_number = int(line[6:11])
            atom_name = line[12:16].strip()
            if atom_name.startswith('H') and exclude_Hs:
                continue
            coor = [float(line[30:38]), float(
                line[38:46]), float(line[46:54])]
            coor = np.array(coor)
            if residue_num not in ligand_dict:
                ligand_dict[residue_num] = list()
            ligand_dict[residue_num] += [coor]

        if line[0: 6] == 'ENDMDL':
            model_dict[model_id] = ligand_dict
    if len(model_dict.keys()) == 0:
        model_dict[1] = ligand_dict

    return model_dict


def gen_conf_ref(m_new_m):
    m_new = Chem.RemoveAllHs(m_new_m)
    m_new_h = Chem.AddHs(m_new)
    num_atoms_mcs = m_new_h.GetNumAtoms()
    coordmap = dict()
    conf0 = m_new_h.GetConformer()

    for i in range(num_atoms_mcs):
        pos = conf0.GetAtomPosition(i)
        if np.linalg.norm((pos.x, pos.y, pos.z)) == 0:
            continue
        coordmap[i] = pos
    cids = AllChem.EmbedMultipleConfs(m_new_h, numConfs=1, numThreads=1,
                                      clearConfs=True, useRandomCoords=True,
                                      coordMap=coordmap)
    return m_new_h


def gen_conf_ref_mol(m, m_ref):
    m_new = Chem.RemoveAllHs(m)
    m_h = Chem.AddHs(m_new)
    check_error = AllChem.EmbedMolecule(m_h)
    if check_error:
        return 'error'

    match = m_h.GetSubstructMatch(m_ref)
    atommap = [[x, i] for i, x in enumerate(match)]
    rmsd = Chem.rdMolAlign.AlignMol(m_h, m_ref, atomMap=atommap)

    return m_h, rmsd


def find_root_atom_idx(m_new, m_ref_com, match_idx=0, atom_idx=0):
    match_atoms_list = m_new.GetSubstructMatches(m_ref_com, uniquify=False)
    root_atom_idx = match_atoms_list[match_idx][atom_idx]
    return root_atom_idx


def pdbqt_to_flex(pdbqt0, flex_pdbqt):
    fp = open(pdbqt0)
    lines = fp.readlines()
    fp.close()

    fp_out = open(flex_pdbqt, 'w')
    residue_name = 'UNL'
    chain_id = ' '
    res_idx = '   1'
    line_begin = 'BEGIN_RES %s %s %s\n' % (residue_name, chain_id, res_idx)
    fp_out.write(line_begin)
    for line in lines:
        if line[0:7] == 'TORSDOF':
            continue
        fp_out.write(line)
    line_end = 'END_RES %s %s %s\n' % (residue_name, chain_id, res_idx)
    fp_out.write(line_end)
    fp_out.close()


def split_pdbqt_dict(pdbqt):
    fp = open(pdbqt)
    lines = fp.readlines()
    fp.close()

    model_dict = dict()

    for line in lines:
        if line[0:5] == 'MODEL':
            model_id = int(line[6:].strip())
            model_dict[model_id] = list()
            continue
        elif line[0:6] == 'ENDMDL':
            continue
        model_dict[model_id].append(line)
    return model_dict


def cal_ligand_size(ligand):
    coor_list = list()
    for atom_coor in ligand:
        coor_list.append(atom_coor)
    coor_list = np.array(coor_list)
    cmin = coor_list.min(axis=0)
    cmax = coor_list.max(axis=0)
    return cmin, cmax


def cal_box(ligand_file_list, exclude_Hs=True):
    cmins = list()
    cmaxs = list()
    for ligand_file in ligand_file_list:
        ligand_model_dict = read_coor_pdb(
            ligand_file, exclude_Hs=exclude_Hs)
        for model_idx in ligand_model_dict.keys():
            ligand_dict = ligand_model_dict[model_idx]
            for ligand_num in ligand_dict.keys():
                ligand_coor = ligand_dict[ligand_num]
                cmin0, cmax0 = cal_ligand_size(ligand_coor)
                cmins.append(cmin0)
                cmaxs.append(cmax0)
    cmins = np.array(cmins)
    cmaxs = np.array(cmaxs)
    cmin = cmins.min(axis=0)
    cmax = cmaxs.max(axis=0)
    return cmin, cmax


def cal_box_size(ligand_file_list, margin=4.0, use_hydrogen=False):
    """
        cal box size from ligands
        input:
            ligand file list
            margin: addtional box-size to ligand size, default 3.0
            use_hydrogen: include hydrogen atom position, defalut Flase
        output:
            box_center : tuple (x, y, z)
            box_size : tuple (wx, wy, wz)
    """
    cmins = list()
    cmaxs = list()
    for ligand_file in ligand_file_list:
        file_format = ligand_file.split(".")[-1]
        ms = list(pybel.readfile(file_format, ligand_file))
        m = ms[0]
        if not use_hydrogen:
            m.removeh()
        atoms = m.atoms
        coor_list = list()
        for atom in atoms:
            coor_list.append(atom.coords)
        coor = np.array(coor_list)

        cmin0 = coor.min(axis=0)
        cmax0 = coor.max(axis=0)
        cmins.append(cmin0)
        cmaxs.append(cmax0)
    cmins = np.array(cmins)
    cmaxs = np.array(cmaxs)
    cmin = cmins.min(axis=0)
    cmax = cmaxs.max(axis=0)

    box_center = tuple((cmax+cmin)/2.0)
    box_size = tuple((cmax-cmin) + margin*2)

    return box_center, box_size


def main():

    import argparse
    title_line = 'convert pdbqt to pdb using reference pdb file'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', required=True,
                        help='input ligand pdbqt file')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output ligand pdb file')
    parser.add_argument('-r', '--ref_file', required=True,
                        help='reference ligand pdb file')

    args = parser.parse_args()
    ligand_input_file = args.input_file
    ligand_output_file = args.output_file
    ref_file = args.ref_file

    e = pdbqt_to_pdb_ref(ligand_input_file, ligand_output_file, ref_file)
    if e is not None:
        print(e)


if __name__ == "__main__":
    main()
