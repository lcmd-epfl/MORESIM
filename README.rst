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
		
Installation 
------------
The best way to install the package is to create our own python environment. 
For this you can use the miniconda approach by typing in your terminal:
	.. line-block::
		conda create --name mymoresim python=3.9 os sys numpy=1.20.3 argparse=1.1 psutil=5.8.0 time cpuinfo scipy=1.6.3 dill=0.3.4 concurrent ase=3.22.0
You can then activate your 
new environment:
	.. line-block::
		conda activate mymoresim
And decide to deactivate it 
once finished:
	.. line-block::
		conda deactivate mymoresim	
Concerning now the different baselined and ML corrections, you should install them independently, depending on
what type fo simulations you would like to perform:
	* DFTB: ensure that the dftb+ executable is in your bin
	* XTB: ensure that the xtb executable is in your bin
	* qml: install the qml package within your python environment
	* pynnp: compile N2P2 with the pynnp module, and export your library in your PATH
	* deepmd: install it using conda, or compile deepmd with shared libraries and export it in your PATH
You can then download the MORESIM 2.0 code directly in your directory
using the followning command:
	.. line-block::
		git clone 
Example
-------
See launch.sh for some trivial first examples
Future good tutorials are in current statement !

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
