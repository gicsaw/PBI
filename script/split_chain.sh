#!/bin/bash

list_file=$1
mkdir chain
lines=`awk -F";" '{print $1}' $list_file`
for code in $lines
do
    echo $code
    split_chain.py -i pdb/$code.pdb -d chain/

done 

