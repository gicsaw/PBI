#!/bin/bash

list_file=$1
dock_dir=dock
mkdir $dock_dir
while read line
do
    array=($line)
    mol_id=${array[0]}
    smi=${array[1]}
    smi_p=`obabel -:$smi -p7.4 -osmi`
    echo $smi_p
    gen3d.py -i $smi_p -o $dock_dir/$mol_id.pdb
    pdb2pdbqt.py -i $dock_dir/$mol_id.pdb -o $dock_dir/$mol_id.pdbqt -l
    qvina02 --config config.txt --ligand $dock_dir/$mol_id.pdbqt --out $dock_dir/dock_$mol_id.pdbqt
#    smina --config config.txt --ligand $dock_dir/$mol_id.pdbqt --out $dock_dir/dock_$mol_id.pdbqt

    pdbqt2pdb_ref.py -i $dock_dir/dock_$mol_id.pdbqt -o $dock_dir/dock_$mol_id.pdb -r $dock_dir/$mol_id.pdb

done < $list_file

