#!/bin/bash

line_pymol="pymol "
while read line
do
    array=($line)
    ligand_file=${array[0]}
    code=${array[0]}
    ligand=${array[1]}
#    echo $ligand_file
    line_pymol="$line_pymol $code"_"$ligand.pdb"
#    line_pymol="$line_pymol $ligand_file"

done < list_new.txt
echo $line_pymol



