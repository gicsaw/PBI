#!/bin/bash

list_file=$1
ref_file=$2

while read line
do
    array=($line)
    ligand_file=${array[0]}
#    code=${array[0]}
#    ligand=${array[1]}
#    echo $ligand_file
    dist_ligand.py -i "$ligand_file" -r $ref_file -c 3
#    dist_ligand.py -i "$code"_$ligand.pdb -r 2P2IA_608.pdb

done < $list_file


