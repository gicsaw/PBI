#!/usr/bin/env python
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS


def main():
    if len(sys.argv) < 3:
        print('dock_rmsd.py mol ref_mol')
        sys.exit()

    m_file = sys.argv[1]
    ref_file = sys.argv[2]
    file_format = m_file.strip().split('.')[-1].lower()
    if file_format == 'pdb':
        m = Chem.MolFromPDBFile(m_file)
    elif file_format == 'mol' or file_format == 'sdf':
        m = Chem.MolFromMolFile(m_file)
    elif file_format == 'mol2':
        m = Chem.MolFromMol2File(m_file)
    if m is None:
        print(m_file, 'is strange')

    ref_format = ref_file.strip().split('.')[-1].lower()
    if ref_format == 'pdb':
        m_ref = Chem.MolFromPDBFile(ref_file)
    if ref_format == 'mol' or ref_format == 'sdf':
        m_ref = Chem.MolFromMOLFile(ref_file)
    elif ref_format == 'mol2':
        m_ref = Chem.MolFromMol2File(ref_file)
    if m_ref is None:
        print(ref_file, 'is strange')

    m = Chem.RemoveHs(m)
    m_ref = Chem.RemoveHs(m_ref)

    smarts = Chem.MolToSmarts(m_ref)
#    z = rdFMCS.FindMCS([m, m_ref])
#    print(z.numAtoms, m.GetNumAtoms(), m_ref.GetNumAtoms())
#    smarts = z.smartsString

    patt = Chem.MolFromSmarts(smarts)
    cs_list = m.GetSubstructMatches(patt, uniquify=False)
#    print(cs_list)
    cs_ref = m_ref.GetSubstructMatch(patt)

    num_atoms = len(cs_ref)
    conformers = m.GetConformers()
    conformer_ref = m_ref.GetConformer()
    num_conf = len(conformers)
    num_hatoms = m.GetNumHeavyAtoms()
    for i_conf in range(num_conf):
        conformer = conformers[i_conf]
        rmsd_min = 99.999
        for cs in cs_list:
            tmp = 0.0
            for i in range(num_atoms):
                idx = cs[i]
                idx_ref = cs_ref[i]
                atom = m.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'H':
                    continue
                p = conformer.GetAtomPosition(idx)
                p_ref = conformer_ref.GetAtomPosition(idx_ref)
                pos = np.array([p.x, p.y, p.z])
                pos_ref = np.array([p_ref.x, p_ref.y, p_ref.z])
                tmp += ((pos-pos_ref)**2).sum()
            rmsd = np.sqrt(tmp/num_hatoms)
            if rmsd < rmsd_min:
                rmsd_min = rmsd
        print(i_conf, rmsd_min)


if __name__ == '__main__':
    main()
