#!/usr/bin/env python3.9
#
"""
Reservoir module containing:
	=> class Reservoir
"""

import numpy as np
import ase
from lib import Trajectory as Trj
from lib import Params as prm
from lib import Printing as Prt
import ase.io as aio
import time

class Reservoir:
	"""A constructor for Reservoirs"""

	def __init__(self, structures_folder, energies, temperature, energy_func, id_rep):
		self.structures_folder = structures_folder
		self.energies = energies
		self.size = len(energies)
		self.temperature = temperature
		self.beta = (prm.kb * self.temperature) ** - 1
		self.energy_func = energy_func
		self.id_rep = id_rep

	def simulation_details(self):
		details = {'temperature': self.temperature}
		return details

	def flush(self):
		pass

	def run(self, init_state, *args):
		np.random.seed()
		idx = np.random.choice(np.arange(self.size))
		mol = aio.read(self.structures_folder + '/{}.xyz'.format(idx+1))
		calc = init_state.calc
		mol.calc = calc
		# Check if periodic
		try:
			f = open('pbc.dat','r')
			a = f.readline()
			b = a.split()
			mol.pbc = (True, True, True)
			mol.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
			f.close()
		except:
			pass
		#
		ener = self.energy_func(mol,dp=1)
		mol.energy = ener
		traj = Trj.Trajectory(state_traj=[mol],generation_details=self.simulation_details())
		time = 0
		velocity = []
		num_at = mol.get_atomic_numbers()
		for i in range(len(num_at)):
			velocity.append([0.0, 0.0, 0.0])
		name = 'ms.rep'+str(self.id_rep)
		Trj.flush(mol,name)
		Prt.print_restart_hres(self.id_rep,mol,i+1,integration='VV',v=velocity)
		return [traj, mol, time, velocity]

