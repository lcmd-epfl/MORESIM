#!/usr/bin/env python3.9
#
"""
Monte Carlo module containing:
	=> class RandomParticleMove
	=> class RandGen
	=> class GaussianVar
	=> class Simulation
	=> function mc_prob
"""

import numpy as np
import time
import ase
import os
import psutil
from lib import Trajectory as Trj
from lib import Printing as Prt
from lib import Params as prm
os.environ['KMP_WARNINGS'] = '0'

class RandomParticleMove:
	"""Constructor of random particle MC moves"""

	def __init__(self, rand_num_gen, part_indexs='all'):
		self.rand_num_gen = rand_num_gen
		self.part_indexs = part_indexs

	def move(self, old_position):
		if self.part_indexs == 'all':
			new_pos = old_position + self.rand_num_gen.generate(old_position.shape)
		return new_pos

class RandGen:
	"""Required form of the generator for this library."""

	def __init__(self):
		pass

	def generate(self, size):
		return None

class GaussianVar(RandGen):
	"""Normally distributed random numbers generator."""

	def __init__(self, loc, var):
		self.loc = loc
		self.var = var

	def generate(self, size):
		return np.random.normal(self.loc, self.var, size)

class Simulation:
	"""A Monte Carlo simulation constructor"""

	def __init__(self, atoms, energy_func, temperature, move_list, move_weight_list=None, dyn_mode='MC', verbose=False):
		self.temperature = temperature
		self.beta = (prm.kb * self.temperature) ** - 1
		self.energy_func = energy_func
		self.move_list = move_list
		self.move_weight_list = move_weight_list
		self.moves_used = []
		self.moves_accepted = []
		self.dyn_mode = dyn_mode
		self.atoms = atoms
		self.ase_molecule = ase.Atoms(atoms)
		self.charges = self.ase_molecule.get_atomic_numbers() 
		self.masses = self.ase_molecule.get_masses()
		self.chemsymb = self.ase_molecule.get_chemical_symbols()
		self.verbose = verbose

	def mc_prob(self, weight_new, weight_old):
		prob = min([1, np.exp(- weight_new + weight_old)])
		return prob

	def simulation_details(self):
		details = {'temperature': self.temperature,'moves_used': self.moves_used,'moves_accepted': self.moves_accepted}
		return details

	def _advance(self, old_state, step, weights, lkr=True, n2p2=False, id_rep=0, p=1):
		move_idx = np.random.choice(list(range(len(self.move_list))), p=self.move_weight_list)
		move = self.move_list[move_idx]
		old_pos = old_state.positions
		old_ener = old_state.energy
		new_pos = move.move(old_pos)
		new_struc = ase.Atoms(self.atoms, new_pos)
		# Check the periodicity
		if os.path.isfile('pbc.dat'):
			f = open('pbc.dat','r')
			a = f.readline()
			b = a.split()
			new_struc.pbc = (True, True, True)
			new_struc.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
			f.close()
		#
		if lkr:
			new_ener, t1, t2 = self.energy_func(new_struc, weights, verbose=self.verbose)
		if n2p2:
			new_ener, t1, t2 = self.energy_func(new_struc, weights, id_rep, p, verbose=self.verbose)
		new_weight = self.beta * new_ener
		old_weight = self.beta * old_ener
		prob = self.mc_prob(weight_new=new_weight, weight_old=old_weight)
		accepted = np.random.rand() < prob
		if not accepted:
			print_accepted = 'False'
			new_pos = old_pos
			new_ener = old_ener
		else:
			print_accepted = 'True'
		new_state = ase.Atoms(self.atoms, positions=new_pos)
		new_state.energy = new_ener
		if os.path.isfile('pbc.dat'):
			new_state.pbc = (True, True, True)
			new_state.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
		if weights == 1.0:
			isprint=True
		else:
			isprint=False
		if isprint:
			print_ener = float(new_ener)
			Prt.print_mc_step(step+1,print_ener,print_accepted)
		return new_state, move_idx, accepted, t1, t2

	def _advance_dp(self, old_state, step, weights, dp):
		move_idx = np.random.choice(list(range(len(self.move_list))), p=self.move_weight_list)
		move = self.move_list[move_idx]
		old_pos = old_state.positions
		old_ener = old_state.energy
		new_pos = move.move(old_pos)
		new_struc = old_state
		new_struc.set_positions(new_pos)
		#print(old_state.get_positions(), new_struc.get_positions())
		#new_struc = ase.Atoms(self.atoms, new_pos)
		# Check the periodicity
		if os.path.isfile('pbc.dat'):
			f = open('pbc.dat','r')
			a = f.readline()
			b = a.split()
			new_struc.pbc = (True, True, True) 
			new_struc.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
			f.close()
		#
		new_ener, t1, t2 = self.energy_func(new_struc,dp,weights,verbose=self.verbose)
		new_weight = self.beta * new_ener
		old_weight = self.beta * old_ener
		prob = self.mc_prob(weight_new=new_weight, weight_old=old_weight)
		accepted = np.random.rand() < prob
		if not accepted:
			print_accepted = 'False'
			new_pos = old_pos
			new_ener = old_ener
		else:
			print_accepted = 'True'
		new_state = new_struc
		#new_state = ase.Atoms(self.atoms, positions=new_pos)
		new_state.energy = new_ener
		if os.path.isfile('pbc.dat'):
			new_state.pbc = (True, True, True)
			new_state.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
		if weights == 1.0:
			isprint=True
		else:
			isprint=False
		if isprint:
			print_ener = float(new_ener)
			Prt.print_mc_step(step+1,print_ener,print_accepted)
		return new_state, move_idx, accepted, t1, t2

	def run(self, init_state, steps, weights=1.0, stride=10, return_last=False, deepmd=True, lkr=False, n2p2=False, id_rep=0, pn2p2=1, *args):
		if weights == 1.0:
			isprint=True
		else:
			isprint=False
		# Check the periodicity
		try:
			f = open('pbc.dat','r')
			a = f.readline()
			b = a.split()
			init_state.pbc = (True, True, True)
			init_state.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])),(float(b[6]), float(b[7]), float(b[8]))])
			f.close()
		except:
			pass
		#
		self.atoms = init_state.get_chemical_symbols()
		np.random.seed()
		if deepmd:
			if self.verbose:
				if isprint:
					print('VERBOSE: DeepMD ML potential has been well detected')
			from deepmd.infer import DeepPot
			dp = DeepPot('graph.pb')
			if 'energy' in dir(init_state):
				ener = init_state.energy
			else:
				ener, tt1, tt2 = self.energy_func(init_state,dp,weights,verbose=self.verbose)
				print_ener = float(ener)
				init_state.energy = print_ener
				if isprint:
					Prt.print_mc_step(0,print_ener,'False')
		if lkr:
			if self.verbose:
				if isprint:
					print('VERBOSE: LKR ML potential has been detected')
			if 'energy' in dir(init_state):
				ener = init_state.energy
			else:
				ener, tt1, tt2 = self.energy_func(init_state,weights,verbose=self.verbose)
				print_ener = float(ener)
				init_state.energy = print_ener
				if isprint:
					Prt.print_mc_step(0,print_ener,'False')
		if n2p2:
			if self.verbose:
				if isprint:
					print('VERBOSE: N2P2 ML potential has been detected')
			if pn2p2 == 1:
				import pynnp
				pn2p2 = pynnp.Prediction()
				pn2p2.log.writeToStdout = False
				pn2p2.setup()
			if 'energy' in dir(init_state):
				ener = init_state.energy
			else:
				ener, tt1, tt2 = self.energy_func(init_state,weights,id_rep,pn2p2,verbose=self.verbose)
				print_ener = float(ener)
				init_state.energy = print_ener
				if isprint:
					Prt.print_mc_step(0,print_ener,'False')
		traj = []
		moves_used = []
		moves_accepted = []
		times = []
		state = init_state
		if self.verbose:
			if isprint:
				print('VERBOSE: Starting MC simulation ...')
		for i in range(steps):
			t1 = time.time()
			if deepmd:
				state, move_idx, accepted, tt1, tt2 = self._advance_dp(state,i,weights,dp)
			if lkr:
				state, move_idx, accepted, tt1, tt2 = self._advance(state,i,weights,lkr=True,n2p2=False,p=1)
			if n2p2:
				state, move_idx, accepted, tt1, tt2 = self._advance(state,i,weights,lkr=False,n2p2=True,id_rep=id_rep,p=pn2p2)
			moves_used.append(move_idx)
			moves_accepted.append(accepted)
			t2 = time.time()
			if (i+1) % stride == 0:
				traj.append(state)
				# Flushing trajectory in case of cMD simulations
				if self.dyn_mode == 'cMD':
					Trj.flush(state,'traj',verbose=self.verbose)
					velocity = []
					acceleration = []
					new_position = state.get_positions()
					Prt.print_restart_cmd(self.charges,self.chemsymb,new_position,velocity,i+1,integration='MC',verbose=self.verbose)
				if self.dyn_mode == 'hRES':
					name = 'ms.rep'+str(id_rep)
					Trj.flush(state,name,verbose=self.verbose)
					# Updating restart file with positions and velocities
					velocity = []
					for j in range(len(self.masses)):
						velocity.append([0.0, 0.0, 0.0]) 
					Prt.print_restart_hres(id_rep,state,i+1,integration='VV',v=velocity,verbose=self.verbose)
				if isprint:
					Prt.print_timing(t1,t2,tt1,tt2,steps,integration='MC',verbose=self.verbose)
			if self.verbose:
				if isprint:
					print('VERBOSE: step {}    Timing (in sec) = {}'.format(i+1,t2-t1))
			if isprint:
				print('STEP {}'.format(i+1))
		self.moves_used = moves_used
		self.moves_accepted = moves_accepted
		if self.verbose:
			if isprint:
				print('VERBOSE: Finishing MC simulation, traj will be printed ...')
		traj = Trj.Trajectory(traj, generation_details=self.simulation_details(),verbose=self.verbose)
		if return_last is True:
			return [traj, state, np.mean(times)]
		else:
			return traj

def mc_prob(weight_new, weight_old):
	prob = min([1, np.exp(- weight_new + weight_old)])
	return prob
