#!/bin/bash

list_file=$1
mkdir select
while read line
do
    code=${line:6:5}
    ligand=${line:12:3}
    echo $code $ligand
    cp chain/"$code"_receptor.pdb chain/"$code"_$ligand.pdb select

done < $list_file


