# cMD MC simultion
#python main.py -p True -dyn cMD -int MC -rep 1 -nstp 2000 -T 300 -freq 1
#python main.py -p True -dyn cMD -int MC -rep 1 -nstp 5 -T 300 -freq 1 -rst True
#python main.py -p True -dyn cMD -int MC -ml LKR -rep 1 -nstp 5 -T 300 -freq 1
#python main.py -p True -dyn cMD -int MC -ml N2P2 -rep 1 -nstp 10 -T 300 -freq 1
#python main.py -p True -dyn cMD -int MC -bsnld XTB -rep 1 -nstp 10 -T 300 -freq 1
# cMD VVL simulation
#python main.py -p True -dyn cMD -int VVL -rep 1 -T 300 -freq 1 -nstp 100 -lgv 0.01
#python main.py -p True -dyn cMD -int VVL -rep 1 -T 300 -freq 1 -nstp 100 -rst True
#python main.py -p True -dyn cMD -int VVL -ml N2P2 -rep 1 -T 300 -freq 1 -nstp 100 -lgv 0.01
#python main.py -p True -dyn cMD -int VVL -bsnld XTB -rep 1 -T 300 -freq 1 -nstp 100 -lgv 0.01
# cMD VV simulation
#python main.py -p True -dyn cMD -int VV -rep 1 -ts 1 -T 300 -freq 10 -nstp 100 -rseed 1897
python main.py -p True -dyn cMD -int VV -rep 1 -ts 1 -T 300 -freq 10 -nstp 100 -rseed 1897 -plm True
#python main.py -p True -dyn cMD -int VV -rep 1 -T 300 -freq 100 -nstp 100 -rst True
#python main.py -p True -dyn cMD -int VV -ml N2P2 -rep 1 -T 300 -freq 1 -nstp 100 
#python main.py -p True -dyn cMD -int VV -bsnld XTB -rep 1 -T 300 -freq 1 -nstp 1000
# cMD VV_MTS simulation
#python main.py -p True -dyn cMD -int VV_MTS -bsnld XTB -ml N2P2 -rep 1 -T 300 -freq 10 -nstp 1000

# hRES MC simulation
#python main.py -p True -dyn hRES -int MC -T 300 -freq 1 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 20 --verbose True
#python main.py -p True -dyn hRES -int MC -T 300 -freq 3 -nstp 3 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 10
#python main.py -p True -dyn hRES -int MC -ml LKR -T 300 -freq 3 -nstp 3 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50
#python main.py -p True -dyn hRES -int MC -ml N2P2 -T 300 -freq 3 -nstp 3 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50
#python main.py -p True -dyn hRES -int MC -bsnld XTB -T 300 -freq 3 -nstp 3 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50
# hRES VVL simulation
#python main.py -p True -dyn hRES -int VVL -T 300 -rep 4 -ts 0.5 -nstp 20 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 20
#python main.py -p True -dyn hRES -int VVL -ml N2P2 -T 300 -rep 4 -ts 0.5 -nstp 3 -freq 3 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 20
# hRES VV simulation
#python main.py -p True -dyn hRES -int VV -T 300 -rep 4 -ts 0.5 -nstp 3 -freq 3 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50
#python main.py -p True -dyn hRES -int VV -ml N2P2 -T 300 -rep 4 -ts 0.5 -nstp 3 -freq 3 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50
#python main.py -p True -dyn hRES -int VV -bsnld XTB -T 300 -rep 4 -ts 0.5 -nstp 3 -freq 3 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50

# USEFUL ANNEXED COMMANDS
# => Activate the deepmd environment: 
#		conda activate deepmd
# => for dftb: import dftb path manually:
		export PATH=$PATH:/home/celerse/1_software/dftbplus/_install/bin/
# => for qml: already installed within anaconda 2019
# => for n2p2, we need to export the python path: 
#		export PYTHONPATH=$PYTHONPATH:/home/celerse/1_software/n2p2/lib/
# => for plumed, we still need to export the python path and import INTELMPI libraries:
#		export PLUMED_KERNEL=/home/celerse/miniconda3/envs/deepmd/lib/libplumedKernel.so
		module load intelmpi/17.0.8
