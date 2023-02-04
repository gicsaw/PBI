from setuptools import setup, find_packages

setup(name='pbi',
        version='0.1',
        packages=['pbi'],
#        packages=find_packages(),
        url='https://github.com/gicsaw/PBI',
        license='MIT LICENSE',
        author='Seung Hwan Hong',
        author_email='gicsaw0@gmail.com',
        description='',
        scripts=['bin/align_3d.py', 'bin/cal_box.py',
                'bin/check_binding_water.py', 'bin/dist_ligand.py',
                'bin/dock_rmsd.py', 'bin/dw_pdb.py', 'bin/filter_pdb_list.py',
                'bin/find_mutation_pdb.py', 'bin/fix_ligand_ref.py',
                'bin/fix_protein.py', 'bin/gen3d.py', 'bin/nwalign.py',
                'bin/pdb2pdbqt.py', 'bin/pdbqt2pdb_ref.py',
                'bin/split_chain.py', 'bin/split_ligand.py',
                'bin/uniprot_dict.py'
               ]
)


