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
                'bin/dist_ligand.py', 'bin/fix_ligand_ref.py',
                'bin/fix_protein.py', 'bin/gen3d.py', 'bin/pdb2pdbqt.py',
                'bin/pdbqt2pdb_ref.py', 'bin/nwalign.py',
                'bin/split_chain.py', 'bin/split_ligand.py',
                'bin/fix_ligand.py', 'bin/dock_rmsd.py'
               ]
)


