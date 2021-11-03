NEW MORESIM 2.0
===============
**Modular Replica Exchange Simulator
VERSION 2.0**

About NEW MORESIM 2.0
---------------------
MORESIM is a made home package that allows to perform standard MD simulations
at a low level of theory (baselined) and correct the obtained PES using some 
ML potentials !
	--> The different MD approaches:
		* Conventional Molecular Dynamics (cMD)
		* Hamiltonian Reservoir Replica Exchanges Molecular Dynamics (hRES)
	--> The different available integrators:
		* Monte Carlo
		* Velocity Verlet Langevin Modified
		* Velocity Verlet
	--> The different available baselined and ML corrections:
		* XTB, DFTB
		* KRR, N2P2, DeepMD

Requirements
------------
Mandatory:
	* Python 3.9 (should also work with 3.6/3.7/3.8 but not tested)
	* Python basic packages:
		- os
		- sys
		- numpy 1.20.3
		- argparse 1.1
		- psutil 5.8.0
		- time 
		- cpuinfo
		- scipy 1.6.3 (deprecated, is useless)
		- dill 0.3.4
		- concurrent
	* Miniconda 3
	* Atomic Simulation Enviroment (ASE) version 3.22.0

Optional (depends on your target):
	1/ For baselined:
		* DFTB: https://dftbplus.org/download			
			- NOTE: DFTB should be compiled with python and DFTD3 options in preference
		* XTB: https://xtb-docs.readthedocs.io/en/latest/setup.html
	2/ For ML corrections:
		* KRR: qml python package https://github.com/qmlcode/qml/blob/master/README.md
		* N2P2: pynnp python package https://github.com/CompPhysVienna/n2p2
		* DeepMD: deepmd python package https://github.com/deepmodeling/deepmd-kit
		
Standard Installation 
---------------------
The best way to install the package is to create our own python environment. 
For this you can use the miniconda approach by typing in your terminal:

.. code:: bash

	~$ conda create --name mymoresim python=3.9 os sys numpy=1.20.3 argparse=1.1 psutil=5.8.0 time cpuinfo scipy=1.6.3 dill=0.3.4 concurrent ase=3.22.0

You can then activate your 
new environment:

.. code:: bash

	~$ conda activate mymoresim

And decide to deactivate it 
once finished:

.. code:: bash

	~$ conda deactivate mymoresim
	
Concerning now the different baselined and ML corrections, you should install them independently, depending on
what type fo simulations you would like to perform:
	* DFTB: ensure that the dftb+ executable is in your bin
	* XTB: ensure that the xtb executable is in your bin
	* qml: install the qml package within your python environment
	* pynnp: compile N2P2 with the pynnp module, and export your library in your PATH
	* deepmd: install it using conda, or compile deepmd with shared libraries and export it in your PATH

You can then download the MORESIM 2.0 code directly in your directory
using the followning command:

.. code:: bash

	~$ git clone https://github.com/lcmd-epfl/MORESIM.gitExample 

The directory where the code was downloaded is thus ready to be used 
for some simulations !

Specific installation in case of DeepMD
---------------------------------------
If DeepMD is the ML you would like target, you have to know that DeepMD-kit
allow you to directly install it using conda. Therefore, you can create your
own conda environment directly with deepmd:

.. code:: bash

	~$ conda create -n deepmd deepmd-kit=*=*cpu libdeepmd=*=*cpu lammps-dp -c https://conda.deepmodeling.org
	
and then add manually each python packages using the pip install command.

Possibility for GPUs plateform is also allowed:

.. code:: bash

	~$ conda create -n deepmd deepmd-kit=*=*gpu libdeepmd=*=*gpu lammps-dp cudatoolkit=11.3 horovod -c https://conda.deepmodeling.org

Fast and small tutorial for hurry users
---------------------------------------
We list here the most important things that a user has to know
in order to correctly use the code.

Examples
--------
Future good tutorials are in current statement !
launch.sh lists some basic commands to launch simulations !
We list here some examples of possible commands. 
Note that it is not the whole possible commands but just use here to show how a computation is basically launched.

**cMD simulations**

- DFTB/DeepMD // Monte Carlo: 

.. code:: bash

	~$ python main.py -p True -dyn cMD -int MC -rep 1 -nstp 2000 -T 300 -freq 1

- DFTB/DeepMD // Restart // Monte Carlo: 
		
.. code:: bash

	~$ python main.py -p True -dyn cMD -int MC -rep 1 -nstp 2000 -T 300 -freq 1 -rst True

- DFTB/KRR // Monte Carlo: 

.. code:: bash

	~$ python main.py -p True -dyn cMD -int MC -ml LKR -rep 1 -nstp 5 -T 300 -freq 1

- DFTB/DeepMD // Velocity Verlet Langevin Modified: 

.. code:: bash

	~$ python main.py -p True -dyn cMD -int VVL -rep 1 -T 300 -freq 1 -nstp 100 -lgv 0.01

- DFTB/DeepMD // Velocity Verlet:

.. code:: bash

	~$ python main.py -p True -dyn cMD -int VV -rep 1 -ts 1 -T 300 -freq 100 -nstp 100 -rseed 1897

- DFTB/N2P2 // Velocity Verlet:

.. code:: bash

	~$ python main.py -p True -dyn cMD -int VV -ml N2P2 -rep 1 -T 300 -freq 1 -nstp 100 

- XTB/N2P2 // Velocity Verlet: 

.. code:: bash

	~$ python main.py -p True -dyn cMD -int VV -bsnld XTB -rep 1 -T 300 -freq 1 -nstp 1000

**hRES simulations**

- DFTB/DeepMD // Monte Carlo:

.. code:: bash

	~$ python main.py -p True -dyn hRES -int MC -T 300 -freq 1 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 20

- DFTB/DeepMD // Restart // Monte Carlo:

.. code:: bash

	~$ python main.py -p True -dyn hRES -int MC -T 300 -freq 1 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 20 -rst True

- XTB/N2P2 // Monte Carlo:

.. code:: bash

	~$ python main.py -p True -dyn hRES -int MC -bslnd XTB -ml N2P2 -T 300 -freq 3 -nstp 3 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50

- XTB/N2P2 // Velocity Verlet:	

.. code:: bash

	~$ python main.py -p True -dyn hRES -int VV -bslnd XTB -ml N2P2 -T 300 -freq 20 -nstp 20 -rep 4 -rsv /home/celerse/ASE-lammps/pool_dithiacyclophene/new_reservoir/ -exc 50

Authors
-------
	* Raimon Fabregat: raimon.fabregat@epfl.ch
	* Frederic Celerse: frederic.celerse@epfl.ch
	* Alberto Fabrizio: alberto.fabrizio@epfl.ch
	* Veronika Juraskova: veronika.juraskova@epfl.ch
	* Benjamin Meyer: benjamin.meyer@epfl.ch
	* Theo Jaffrelot Inizant: theo.jaffrelot-inizant@sorbonne-universite.fr
	* Daniel Hollas: daniel.hollas@epfl.ch
	* Clemence Corminboeuf: clemence.corminboeuf@epfl.ch
