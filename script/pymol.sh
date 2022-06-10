#!/bin/bash

list_file=$1
line_pymol="pymol"
while read line
do
    array=($line)
    ligand_file=${array[0]}
#    code=${array[0]}
#    ligand=${array[1]}
#    echo $ligand_file
    line_pymol+=" $ligand_file"

done < $list_file
echo $line_pymol



