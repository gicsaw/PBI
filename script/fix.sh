#!/bin/bash

list_file=$1
ref_file=$2
mkdir fix
while read line
do
    array=($line)
#    ligand_file=${array[0]}
    code=${array[0]}
    ligand=${array[1]}
#    echo $ligand_file
#    fix_protein.py -i select/"$code"_receptor.pdb -o fix/"$code"_receptor.pdb
    fix_ligand_ref.py -i select/"$code"_$ligand.pdb -r cif/$ligand.cif -o fix/"$code"_$ligand"_p.pdb" -n -p 7.4 --fix_atom_idx
    fix_ligand_ref.py -i select/"$code"_$ligand.pdb -r cif/$ligand.cif -o fix/"$code"_$ligand.pdb --fix_atom_idx

#    fix_ligand.py -i select/"$code"_$ligand.pdb -o fix/"$code"_$ligand.pdb -n -p 7.4 --fix_atom_idx

#    obabel fix/"$code"_$ligand.pdb -O fix/"$code"_$ligand"_p.pdb"
#    pdb2pdbqt.py -i fix/"$code"_$ligand"_p.pdb" -o fix/"$code"_$ligand"_p.pdbqt" -l

done < $list_file


