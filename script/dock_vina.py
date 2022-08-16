#!/usr/bin/env python
import sys
import os
import subprocess
from pbi import ligandtools
from pbi.pdbtools import PDBtools
from openbabel import pybel


def add_hydrogen(smi0, pH=7.4):
    try:
        m = pybel.readstring("smi", smi0)
        m.OBMol.AddHydrogens(False, True, pH)
        smi_p = m.write("smi").strip()
    except Exception:
        smi_p = smi0
    return smi_p


def docking(vina, config_file, ligand_file, output_file):

    run_line = '%s --config %s --ligand %s --out %s' % (
        vina, config_file, ligand_file, output_file)
    e = None
    try:
        subprocess.check_output(run_line.split(),
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)
    except Exception as e:
        return e
    return e


def main():

    if len(sys.argv) < 5:
        print('dock_vina.py vina_file config_file smi_list_file out_dir')
        sys.exit()

    vina = sys.argv[1]
    config_file = sys.argv[2]
    list_file = sys.argv[3]
    out_dir = sys.argv[4]
    pH = 7.4

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        lis = line.strip().split()
        mol_id = lis[0]
        smi = lis[1]
        print(mol_id)
        pdb_file = out_dir + '/' + mol_id + '.pdb'
        pdbqt_file = out_dir + '/' + mol_id + '.pdbqt'
        dock_pdbqt_file = out_dir + '/' + '/dock_' + mol_id + '.pdbqt'
        dock_pdb_file = out_dir + '/' + '/dock_' + mol_id + '.pdb'

        if pH is not None:
            smi_p = add_hydrogen(smi, pH=pH)
        else:
            smi_p = smi
        ligandtools.gen_3d(smi_p, pdb_file, mol_id=mol_id, file_format='pdb',
                           timeout=20)
        PDBtools.ligand_to_pdbqt(pdb_file, pdbqt_file)
        docking(vina, config_file, pdbqt_file, dock_pdbqt_file)
        ligandtools.pdbqt_to_pdb_ref(dock_pdbqt_file, dock_pdb_file, pdb_file)


if __name__ == "__main__":
    main()
