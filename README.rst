NEW MORESIM 2.0
===============
**Modular Replica Exchange Simulator
VERSION 2.0**

About NEW MORESIM 2.0
---------------------
| MORESIM is a made home package that allows to perform standard MD simulations
at a low level of theory (baselined) and correct the obtained PES using some 
ML potentials !
| --> The different MD approaches:
	* Conventional Molecular Dynamics (cMD)
	* Hamiltonian Reservoir Replica Exchanges Molecular Dynamics (hRES)
| --> The different available integrators:
	* Monte Carlo
	* Velocity Verlet Langevin Modified
	* Velocity Verlet
| --> The different available baselined and ML corrections:
	* XTB, DFTB
	* KRR, N2P2, DeepMD

Requirements
------------
Mandatory:
* Python_ 3.9 (should also work with 3.6/3.7/3.8 but not tested)
* Python_ packages: os, 
* Anaconda_ 3
* Atomic Simulation Enviroment (ASE_)

Optional (depends on your target):
1/ For baselined:
	* DFTB: dftb+ software
	* XTB: xtb software
2/ For ML corrections:
	* KRR: qml python package
	* N2P2: pynnp python package
	* DeepMD: deepmd python package

Installation 
------------
TOTO IS TOTO :)

Example
-------
See launch.sh for some trivial first examples
Future good tutorials are in current statement !

Authors
-------
Raimon Fabregat: raimon.fabregat@epfl.ch
Frederic Celerse: frederic.celerse@epfl.ch
Alberto Fabrizio: alberto.fabrizio@epfl.ch
Veronika Juraskova: veronika.juraskova@epfl.ch
Benjamin Meyer: benjamin.meyer@epfl.ch
Theo Jaffrelot Inizant: theo.jaffrelot-inizant@sorbonne-universite.fr
Daniel Hollas: daniel.hollas@epfl.ch
Clemence Corminboeuf: clemence.corminboeuf@epfl.ch
