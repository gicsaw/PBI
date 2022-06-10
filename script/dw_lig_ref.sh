#!/bin/bash

list_file=$1
mkdir cif
line_cif="wget"

while read line
do
    array=($line)
    code=${array[0]}
    ligand=${array[1]}
    echo $code $ligand
    line_cif+=" https://files.rcsb.org/ligands/view/$ligand.cif"

done < $list_file
$line_cif -P cif



