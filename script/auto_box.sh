#!/bin/bash

list_file=$1
autobox_file='autobox.txt'
echo "" > $autobox_file
while read line
do
    array=($line)
    code=${array[0]}
    ligand=${array[1]}
    echo "--autobox_file_list=fix/"$code"_$ligand.pdb" >> $autobox_file

done < $list_file
cal_box.py --arg_file $autobox_file -m 4


