#!/bin/bash

### see help: uniprot_dict.py 
# echo "db_dir=$UniPort_db_dir" > db_dir.txt
### db initialize
### uniprot_dict.py make_pickle


### find target protein at uniprot https://www.uniprot.org 
### uniprot_dict.py pdb UniProt_ID
### example:
uniprot_dict.py pdb EGFR_HUMAN

### filtering pdb list
### option: resolution (--resolution value), include NMR (--include_nmr), sequence range (-r ini-fin)...
### filter_pdb_list.py --resolution value -r ini-fin list_file > output_file
### example:
filter_pdb_list.py --resolution 10 P00533_pdb_list.txt > P00533_pdb_list_filter.txt

### download pdb files to pdb_dir (default: pdb)
### dw_pdb.py $list_file $pdb_dir
### example:
dw_pdb.py P00533_pdb_list_filter.txt pdb

### check mutation and ligand 
### find_mutation_pdb.py $list_file $fasta_file $pdb_dir > list_mutation.txt
### example:
find_mutation_pdb.py P00533_pdb_list_filter.txt P00533.fasta pdb > list_mutation.txt

### split chain, result is chain folder 
### ./split_chain.sh $list_file
### example:
./split_chain.sh list_mutation.txt

### generate chain list file from pdb list file 
### option: --use_all_chain: use all chain
### python gen_chain_list.py -i $input_list_file -o $output_list_file
### example:
python gen_chain_list.py -i list_mutation.txt -o list_chain.txt

### align pdb by reference pdb 
### ./rot.sh $list_chain_file $ref_pdb_file
### example:
./rot.sh list_chain.txt chain/2ITYA.pdb

### split receptor and ligand 
### split_ligand.sh $list_file
### example:
./split_ligand.sh list_chain.txt

#gather ligand_list 
### if you want only chain A, ls chain/????A_???.pdb > list_ligand_A.txt
### or grep A_ list_ligand.txt > list_ligand_A.txt
### example:
ls chain/?????_???.pdb |grep -v HOH > list_ligand.txt

### find pocket ligand using reference ligand
###./dist.sh $ligand_list_file $ref_lig.pdb > list_ligand_new.txt
### example
./dist.sh list_ligand.txt chain/2ITYA_IRE.pdb > list_ligand_new.txt
cp list_ligand_new.txt list_ligand_select.txt

### show ligand 
./pymol.sh list_ligand_new.txt
### run_output script and check ligands
### fix list_ligand_select.txt  manualy 
### vi list_ligand_select.txt

### copy pdb files to select folder
./cp_pdb.sh list_ligand_select.txt > list_final.txt

### check bound waters
### example: 
./check_water.sh list_ligand_select.txt
ls select/*_HOH.pdb > water_file_list.txt

### check consensus water : I don't recommand... 
consensus_water_coor.py list_ligand_select.txt > consensus_water_list.txt
#select water in consensus_water_list and remove other waters
find_consensus_water.py select/2ITYA_HOH.pdb consensus_water_list.txt

### copy and paste HOH to receptor file 

### see select directory

### download ligand cif file (reference)
./dw_lig_ref.sh list_final.txt
### fix ligand using cif reference file
./fix.sh list_final.txt

### fix_protein.py -i $input_receptor_pdb -o $output_receptor_pdb
### example:
fix_protein.py -i select/2ITYA_receptor.pdb -o fix/2ITYA_receptor.pdb

### pdb2pdbqt
### option -r: receptor, -l: ligand
### pdb2pdbqt -i $input_pdb -o $output_pdbqt -r (or -l)
### example
pdb2pdbqt.py -i fix/2ITYA_receptor.pdb -o fix/2ITYA_receptor.pdbqt -r
pdb2pdbqt.py -i fix/2ITYA_IRE.pdb -o fix/2ITYA_IRE.pdbqt -l

### generate docking box parameters for Vina
./auto_box.sh list_final.txt > config.txt

#cd fix
#qvina02 --config config.txt --ligand 2ITYA_IRE.pdbqt --out a.pdbqt
#pdbqt2pdb_ref.py -i a.pdbqt -o a.pdb -r 2ITYA_IRE.pdb
