#!/bin/bash

list_file=$1
ref_file=$2
while read line
do
    array=($line)
    #protein=${array[0]}
    code=${array[0]}
    chain_id=${array[1]}
    echo $code $chain_id
    code2=$code$chain_id
    align_3d.py -i chain/$code2.pdb -r $ref_file -o chain/"$code2"_rotate.pdb
done < $list_file

