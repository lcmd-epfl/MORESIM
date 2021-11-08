#!/usr/bin/env python3.9
#
"""
Plumed interface to couple MD with Plumed package
[1] The PLUMED consortium, Nat. Methods 16, 670 (2019)
[2] Tribello, Bonomi, Branduardi, Camilloni and Bussi, Comput. Phys. Commun. 185, 604 (2014)
"""

from os.path import exists
from ase.units import fs, mol, nm
import numpy as np

class CalcPlumed:
	""" Plumed general class """

	def __init__(self, input, timestep, atoms=None, log='', restart=False, use_charge=False, update_charge=False, verbose=False):
		import plumed
		if atoms is None:
			raise TypeError('ERROR: plumed has to be defined with the object atoms inside')
		self.istep = 0
		self.input = input
		self.use_charge = use_charge
		self.update_charge = update_charge
		natoms = len(atoms.get_positions())
		# Units setup
		# WARNINGS: outputs from plumed will still be in plumed units
		self.plumed = plumed.Plumed()
		ps = 1000 * fs
		#self.plumed.cmd("setMDEnergyUnits", mol, KJ)
		self.plumed.cmd("setMDLengthUnits", 1/nm)
		self.plumed.cmd("setMDTimeUnits", 1/ps)
		self.plumed.cmd("setMDChargeUnits", 1.)
		self.plumed.cmd("setMDMassUnits", 1.)
		self.plumed.cmd("setMDMassUnits", 1.)
		self.plumed.cmd("setNatoms", natoms)
		self.plumed.cmd("setMDEngine", "python")
		self.plumed.cmd("setLogFile", log)
		self.plumed.cmd("setTimestep", float(timestep))
		self.plumed.cmd("setRestart", restart)
		self.plumed.cmd("setKbT", 1.)
		self.plumed.cmd("init")
		with open(input) as f:
			for line in f:
				self.plumed.cmd("readInputLine", line)
		self.atoms = atoms

	def compute_bias(self, pos, istep, unbiased_energy):
		self.plumed.cmd("setStep", istep)
		self.plumed.cmd("setPositions", pos)
		self.plumed.cmd("setEnergy", unbiased_energy)
		self.plumed.cmd("setMasses", self.atoms.get_masses())
		forces_bias = np.zeros((self.atoms.get_positions()).shape, dtype=float)
		self.plumed.cmd("setForces", forces_bias)
		virial = np.zeros((3,3), dtype=float)
		self.plumed.cmd("setVirial", virial)
		self.plumed.cmd("prepareCalc")
		self.plumed.cmd("performCalc")
		energy_bias = np.zeros((1,))
		self.plumed.cmd("getBias", energy_bias)
		return energy_bias, forces_bias

	def write_plumed_files(self, images):
		return self.read_plumed_files()

	def read_plumed_files(self, file_name=None):
		read_files = {}
		if file_name is not None:
			read_files[file_name] = np.loadtxt(file_name, unpack=True)
		else:
			for line in self.input:
				if line.find('FILE') != -1:
					ini = line.find('FILE')
					end = line.find(' ', ini)
					if end == 1:
						file_name = line[ini+5]
					else:
						file_name = line[ini+5:end]
					readfiles[file_name] = np.loadtxt(file_name, unpack=True)
				if len(read_files) == 0:
					if exists('COLVAR'):
						read_files['COLVAR'] = np.loadtxt('COLVAR', unpack=True)
					if exists('HILLS'):
						read_files['HILLS'] = np.loadtxt('HILLS', unpack=True)
			assert not len(read_files) == 0,"There are not files for reading"
			return read_files

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.plumed.finalize()

