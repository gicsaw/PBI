#!/bin/bash

list_file=$1
mkdir select
while read line
do
    code=${line:6:5}
    ligand=${line:12:3}
    echo $code $ligand
    check_binding_water.py -r chain/$code"_receptor.pdb" -l chain/$code"_"$ligand.pdb -w chain/$code"_HOH.pdb" -o select/$code"_HOH.pdb"

done < $list_file


