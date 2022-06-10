#!/bin/bash

list_file=$1
#ref_file=$2
while read line
do
    array=($line)
    #protein=${array[0]}
#    chain_file=${array[0]}
    pdb_code=${array[0]}
    chain_id=${array[1]}
    chain_file=chain/$pdb_code$chain_id"_rotate.pdb"
    echo $chain_file
    split_ligand.py -i $chain_file -d chain/
#    split_ligand.py -i $chain_file -r $ref_file -d chain/


done < $list_file
#ls ????A_???.pdb > list_ligand.txt

