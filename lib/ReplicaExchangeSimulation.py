#!/usr/bin/env python3.9
#
"""
Replica Exchange classes to define RE simulation
	=> REXC class
	=> _smap function
"""

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import dill
import os
import time
import ase
import psutil
from lib import MC
from lib import VelocityVerlet_Langevin as VVL
from lib import Trajectory as Trj
from lib import Printing as Prt

class REXC:
	"""A general constructor for Replica Exchange Simulation"""
	
	def __init__(self,num_reps,simulations,init_states,weights,stride=20,num_processors=1,rep_steps=20,integration='MC',parallel=True,xtb=False,deepmd=True,lkr=False,n2p2=False,verbose=False):
		directory = '.'
		self.num_reps = num_reps
		self.num_processors = num_processors
		self.temperatures = [sim.temperature for sim in simulations]
		self.energy_funcs = [sim.energy_func for sim in simulations]
		self.simulations = simulations
		self.in_rep_states = init_states
		self.integration = integration
		self.stride=stride
		self.weights = weights
		#self.par_exec = ProcessPoolExecutor(max_workers=num_processors)
		self.par_exec = ThreadPoolExecutor(max_workers=num_processors)
		self.xtb = xtb
		self.deepmd = deepmd
		self.lkr = lkr
		self.n2p2 = n2p2
		self.verbose = verbose
		# toto is used here as passing method for NN calculator for each 
		# replica as it is not pickabled !
		toto = []
		for i in range(num_reps):
			toto.append(1)
		self.toto = toto
		# We have to redefine the weights in order to let them pickable
		# => baselined will still be 1.0, so we declare a 1D array with
		#    the different values from the self.weights
		weights_ml = np.zeros(num_reps, dtype=float)
		for i in range(len(weights_ml)):
			weights_ml[i] = self.weights[i][1]
		self.weights_ml = weights_ml
		if self.verbose:
			a = dill.pickles(self.energy_funcs)
			b = dill.pickles(self.in_rep_states)
			c = dill.pickles(weights_ml)
			print('VERBOSE: is energy_funcs pickabled  ?   => {}'.format(a))
			print('VERBOSE: is in_rep_states pickabled ?   => {}'.format(b))
			print('VERBOSE: is new_weights pickabled   ?   => {}'.format(c))
			if a and b and c:
				print('VERBOSE: all is pickable, it can be thus distributed on {} processors !'.format(num_processors))
			else:
				print('VERBOSE: one object is not pickable ! It can not be distributed on {} processors !'.format(num_processors))
		# Define the number of each replicas 
		id_of_reps = []  
		for i in range(num_reps):
			id_of_reps.append(i) 
		self.id_of_reps = id_of_reps
		if sum(['energy' in dir(state) for state in init_states]) == self.num_reps:
			pass
		else:
			if deepmd:
				self.in_rep_eners = list(self.par_exec.map(_smap, self.energy_funcs, self.in_rep_states, self.toto, self.weights_ml))
			if lkr:
				self.in_rep_eners = list(self.par_exec.map(_smap, self.energy_funcs, self.in_rep_states, self.weights_ml))
			if n2p2:
				self.in_rep_eners = list(self.par_exec.map(_smap, self.energy_funcs, self.in_rep_states, self.weights_ml, self.id_of_reps))
		self.rep_index = np.arange(self.num_reps)
		self.even_sims = self.rep_index[::2]
		self.odd_sims = self.rep_index[::2]
		self.accepted_exchanges = {(i, (i + 1) % self.num_reps):[] for i in range(self.num_reps)}
		self.strides = [stride for i in range(num_reps)]
		self.rep_steps = rep_steps
		for stride in self.strides:
			if self.rep_steps % stride != 0:
				raise ValueError('Rep_steps must be multiple of stride')
		self.rep_stepss = [rep_steps for i in range(self.num_reps)]
		self.directory = directory

	def run(self, num_exchanges, return_last=True, dftbs=1):
		trajectories = [Trj.Trajectory(generation_details=sim.simulation_details()) for sim in self.simulations]
		mem = psutil.Process().memory_info().rss * 1e-9
		print("Memory allocated at the moment 'RUNNING REPLICA EXCHANGES': {} GB".format(mem))
		if self.n2p2:
			import pynnp
			p = pynnp.Prediction()
			p.log.writeToStdout = False
			p.setup()
			pn2p2 = [p for i in range(self.num_reps)]
		for i in range(num_exchanges):
			t1 = time.time()
			return_last = [True for l in range(self.num_reps)]
			if self.deepmd:
				deepmd = [True for l in range(self.num_reps)]
			else:
				deepmd = [False for l in range(self.num_reps)]
			if self.lkr:
				lkr = [True for l in range(self.num_reps)]
			else:
				lkr = [False for l in range(self.num_reps)]
			if self.n2p2:
				n2p2 = [True for l in range(self.num_reps)]
			else:
				n2p2 = [False for l in range(self.num_reps)]
			mem = psutil.Process().memory_info().rss * 1e-9
			print("Memory allocated at the moment 'BEFORE ALLOCATION': {} GB".format(mem))
			# We launch MD within each replicas
			if self.integration == 'MC':
				if self.deepmd or self.lkr:
					simulation_results = list(self.par_exec.map(_run_simulation, self.simulations, self.in_rep_states, self.rep_stepss, self.weights_ml, self.strides, return_last, deepmd, lkr, n2p2, self.id_of_reps))
				if self.n2p2:
					simulation_results = list(self.par_exec.map(_run_simulation, self.simulations, self.in_rep_states, self.rep_stepss, self.weights_ml, self.strides, return_last, deepmd, lkr, n2p2, self.id_of_reps, pn2p2))
				# SEQUENTIAL
				#simulation_results = []
				#for j in range(4):
				#	simulation_results.append(_run_simulation(self.simulations[j], self.in_rep_states[j], self.rep_stepss[j], self.weights_ml[j], self.strides[j], return_last[j], deepmd[j], lkr[j], n2p2[j], self.id_of_reps[j], pp[j]))
			else:
				rseed = np.ones(len(self.simulations), dtype=float)
				if self.deepmd or self.lkr:
					simulation_results = list(self.par_exec.map(_run_simulation, self.simulations, self.in_rep_states, self.rep_stepss, self.strides, rseed, return_last, self.weights_ml, self.id_of_reps, deepmd, n2p2))
				if self.n2p2:
					simulation_results = list(self.par_exec.map(_run_simulation, self.simulations, self.in_rep_states, self.rep_stepss, self.strides, rseed, return_last, self.weights_ml, self.id_of_reps, deepmd, n2p2, pn2p2))
				# SEQUENTIAL
				#simulation_results = []
				#for j in range(4):
				#	simulation_results.append(_run_simulation(self.simulations[j], self.in_rep_states[j], self.rep_stepss[j], self.strides[j], rseed[j], return_last[j], self.weights_ml[j]))
			# From the list we now treat the exchange part
			mem = psutil.Process().memory_info().rss * 1e-9
			print("Memory allocated at the moment 'AFTER ALLOCATION': {} GB".format(mem))
			rep_trajs = [res[0] for res in simulation_results]
			exchange_states = [res[1] for res in simulation_results]
			sim_times = [res[2] for res in simulation_results]
			if self.integration != 'MC':
				velocities = [res[3] for res in simulation_results]
			else:
				velocities = []
			if self.integration == 'MC':
				velocities = []
			Prt.print_hres_iniexc(i)
			print('Exchange {} starting ...'.format(i+1))
			if self.lkr or self.deepmd:
				self.in_rep_states = self.replica_exchange(exchange_states, velocities, dftbs, pn2p2=1)
			if self.n2p2:
				self.in_rep_states = self.replica_exchange(exchange_states, velocities, dftbs, pn2p2)
			print('Exchange {} is done !'.format(i+1))
			mem = psutil.Process().memory_info().rss * 1e-9
			print("Memory allocated at the moment 'AFTER EXCHANGE': {} GB".format(mem))
			self.exchange_probabilities = {key: (0.001 + sum(val)) / (len(val) + 0.001) for key, val in self.accepted_exchanges.items()}
			Prt.print_hres_probaexc(self.exchange_probabilities)

	def replica_exchange(self, exchange_states, velocities, dftbs, pn2p2=1):
		# Initializing parameters for exchange
		shift = np.random.choice([1, -1])
		rep_index = np.arange(self.num_reps)
		group1 = rep_index[::2]
		group2 = rep_index[1::2]
		exchange_structs = [xx.positions for xx in exchange_states]
		exchange_eners = [xx.energy for xx in exchange_states]
		if self.integration != 'MC':
			exchange_vels = velocities
		if shift == 1:
			ex_index = np.vstack((group2, group1)).flatten(order='F')
		else:
			ex_index = np.roll(np.vstack((group1, np.roll(group2, 1))).flatten(order='F'), -1)
		pairs = list(zip(group1, ex_index[::2]))
		old_structs = exchange_structs
		old_energies = exchange_eners
		new_structs = [old_structs[i] for i in ex_index]
		if self.integration != 'MC':
			old_velocities = np.array(exchange_vels)
			new_velocities = [old_velocities[i] * np.sqrt(self.temperatures[j] / self.temperatures[i])for j, i in enumerate(ex_index)]
		new_states = [ase.Atoms(symbols=exchange_states[0].get_chemical_symbols(), positions=new_structs[i]) for i in range(self.num_reps)]
		if self.xtb:
			for i in range(self.num_reps):
				new_states[i].set_calculator(dftbs[i])
		# Check if periodicity
		for i in range(self.num_reps):
			try:
				f = open('pbc.dat','r')
				periodic=True
				a = f.readline()
				b = a.split()
				new_states[i].pbc = (True, True, True)
				new_states[i].set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
				f.close() 
				periodic1 = (True, True, True)
				periodic2 = [(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))]
			except:
				periodic=False
		#
		mem = psutil.Process().memory_info().rss * 1e-9
		print("Memory allocated at the moment 'BEFORE ENER COMPUT': {} GB".format(mem))
		if self.deepmd:
			new_energiess = list(self.par_exec.map(_smap, self.energy_funcs, new_states, self.toto, self.weights_ml))
			# SEQUENTIAL
			#new_energiess = []
			#for j in range(4):
			#	new_energiess.append(self.energy_funcs[j](new_states[j], self.toto[j], self.weights_ml[j]))
		if self.lkr:
			new_energiess = list(self.par_exec.map(_smap, self.energy_funcs, new_states, self.weights_ml))
		if self.n2p2:
			new_energiess = list(self.par_exec.map(_smap, self.energy_funcs, new_states, self.weights_ml, self.id_of_reps, pn2p2))
		# WARNING: In this case, new_energiess will return as a tuple due to the t1 and t2 variables returned in the energy function !
		# => In this case, we redefine manually new_energiess in order to match with old_energies !
		mem = psutil.Process().memory_info().rss * 1e-9
		print("Memory allocated at the moment 'AFTER ENER COMPUT': {} GB".format(mem))
		new_energies = []
		for i in range(len(new_energiess)):
			a = new_energiess[i]
			# In case of Hamiltonian => tuple // in case of reservoir => float
			# => So we have to try here ...
			try:
				b = a[0]
				new_energies.append(b)
			except:
				new_energies.append(a)
		# Printing energetics informations before the exchange
		Prt.print_hres_exc1(old_energies,new_energies)
		# Performing the exchange part of RE
		for pair in pairs:
			rep0 = self.simulations[pair[0]]
			rep1 = self.simulations[pair[1]]
			old_e0 = old_energies[pair[0]]
			old_e1 = old_energies[pair[1]]
			new_e0 = new_energies[pair[0]]
			new_e1 = new_energies[pair[1]]
			old_weight = rep0.beta * old_e0 + rep1.beta * old_e1
			new_weight = rep0.beta * new_e0 + rep1.beta * new_e1
			prob = MC.mc_prob(weight_new=new_weight, weight_old=old_weight)
			accepted = np.random.rand() < prob
			if shift == 1:
				self.accepted_exchanges[(pair[0], pair[1])].append(accepted)
			else:
				self.accepted_exchanges[(pair[1], pair[0])].append(accepted)
			if accepted:
				pass
			else:
				new_structs[pair[0]] = old_structs[pair[0]]
				new_structs[pair[1]] = old_structs[pair[1]]
				new_energies[pair[0]] = old_energies[pair[0]]
				new_energies[pair[1]] = old_energies[pair[1]]
				try:
					new_velocities[pair[0]] = old_velocities[pair[0]]
					new_velocities[pair[1]] = old_velocities[pair[1]]
				except:
					pass
			# Printing exchange informations in the control file
			Prt.print_hres_exc2(pair,prob,accepted)
		# Initializing replicas for next loop
		new_states = [ase.Atoms(symbols=exchange_states[0].get_chemical_symbols()) for i in range(self.num_reps)]
		if self.xtb:
			for i in range(self.num_reps):
				new_states[i].set_calculator(dftbs[i])
		for i, state in enumerate(new_states):
			state.positions = new_structs[i]
			state.energy = new_energies[i]
			if self.integration != 'MC':
				state.velocities = new_velocities[i]
			if periodic:
				state.pbc = periodic1
				state.set_cell(periodic2)
		return new_states

def _smap(f, *args):
	return f(*args)

def _run_simulation(simulation, *args):
	return simulation.run(*args)
