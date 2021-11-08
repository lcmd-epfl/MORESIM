#!/usr/bin/env python3.9
#
"""
Langevin dynamics modified Velocity Verlet algorithm
from Martin Kroger
https://www.complexfluids.ethz.ch/Langevin_Velocity_Verlet.pdf
"""

import numpy as np
from lib import Maxwell_boltzmann_distribution as MBD
import time
import os
import ase
import sys
from deepmd.infer import DeepPot
from lib import Trajectory as Trj
from lib import Printing as Prt
from lib import Params as prm
os.environ['KMP_WARNINGS'] = '0'

class Simulation:
	""" Class that propagate the atomic positions by a Velocity Verlet
	algorithm for any representation class """

	def __init__(self, time_step, atoms, energy_func, langevin_friction_coeff=0.01, temperature=300, dyn_mode='cMD', use_plumed=False, verbose=False):
		"""
		This is small note about the significance and respective units of each variables declared here
		in the Simulation class:
			=> time_step 			= step used to propagate the dynamics 		(in fs)
			=> atoms			= atom types in the dynamics
			=> energy_force_calculator	= function (in EnergyCalculator.py --> MixedPotentials) 
							used to compute potential energy and forces 
			=> langevin_friction_coeff	= langevin damping coefficient			(in fs-1)
			=> temperature			= targetted temperature through the dynamics	(in K)
			=> dyn_mode			= dynamics mode (useful for the flush option)
			=> use_plumed			= plumed coupling
			=> verbose			= debugging mode
		"""
		self.atoms = atoms
		self.time_step = time_step
		self.ase_molecule = ase.Atoms(atoms)
		self.masses = self.ase_molecule.get_masses()
		self.charges = self.ase_molecule.get_atomic_numbers()
		self.chemsymb = self.ase_molecule.get_chemical_symbols()
		self.energy_func = energy_func
		self.langevin_friction_coeff = langevin_friction_coeff
		self.temperature = temperature
		self.dyn_mode = dyn_mode
		self.use_plumed = use_plumed
		self.verbose = verbose
		self.beta = (prm.kb * self.temperature) ** - 1

	def simulation_details(self):
		return vars(self)

	def simulation_type(self):
		return Trj.Trajectory(generation_details=self.simulation_details())

	def run(self, init_state, steps, stride=10, rseed=1, return_last=True, weights=1.0, id_rep=0, deepmd=True, n2p2=False, pn2p2=1, *args):
		np.random.seed()
		"""
		This is small note about how this attribute work.
		It starts from a structures passing in argument and propagate the Velocity Verlet
		Algorithm following a modified Langevin Version. Therefore, it does not need
		any use of thermostat !
		The routine works like that:
			1/ Initialize velocities from
				a) already existed structures
				b) nothing => Maxwell Boltzmann distribution is thus used
			2/ Compute initial forces and accelerations => update new velocities
			3/ Start the propagation along the desired number of steps:
				a) Compute Langevin coefficients (there are constant, so it 
				could be done out of the loop
				b) Choose eta and compute v1/2 => update new positions then
				c) Compute new forces and accelerations => update velocities
				d) compute kinetic energy and T inst
				e) Append stuffs and print stuffs !
		WARNING: an important note about the units which are used here:
			=> These are Standard Units in most of the cases !
			=> velocity     	= m.s-1
			=> acceleration 	= m.s-2
			=> potential energy     = kcal/mol
			=> forces       	= J.m-1
			=> coeff1       	= --
			=> coeff2       	= m.s-1
			=> coeff3       	= fs
			=> kinetic energy	= J
			=> Temperature		= K
			###
			=> potential energy is done in kcal/mol for hRES reason ! so for 
			printing we decided to convert kinetic energy in kcal/mol !
			=> However, to compute T_inst we need to keep it in J !
		"""
		##
		if weights == 1.0:
			isprint=True
		else:
			isprint=False
		##
		if self.verbose:
			print('VERBOSE: Starting Velocity Verlet Langevin modified dynamics !')
			print('VERBOSE: setupping initial velocities')
		if 'velocities' in dir(init_state):
			if self.verbose:
				print('VERBOSE: Velocities detected from a previous initialization')
			input_velocity = init_state.velocities
			# Check if velocity is coming from the resrevoir (cf null)
			if input_velocity[0][0] == 0.0:
				if self.verbose:
					print('VERBOSE: first vx/vy/vz equal to 0 => should come from the reservoir')
					print('VERBOSE: => will be reinitialized using the Maxwell Boltzmann distribution !')
				input_velocity = MBD.maxwell_boltzmann_distribution(self.atoms, self.temperature, rseed, isprint)
		else:
			if self.verbose:
				print('VERBOSE: No initial velocities have been detected')
				print('VERBOSE: => will be reinitialized using the Maxwell Boltzmann distribution !')
			input_velocity = MBD.maxwell_boltzmann_distribution(self.atoms, self.temperature, rseed, isprint)
		# STEP 0 => Initialization of the dynamics
		langevin_friction_coeff = self.langevin_friction_coeff
		temperature = self.temperature
		numb_iterations = steps
		positions = []
		velocity = []
		T = []
		times = []
		accelerations = []
		potential_energies = []
		kinetic_energies = []
		total_energies = []
		forces = []
		traj = []
		positions.append(init_state.get_positions())
		### Compute initial energy/forces --> initial accelerations --> new velocities ###
		if deepmd:
			dp = DeepPot('graph.pb')
			potential, force, tt1, tt2 = self.energy_func(init_state,dp,weights,verbose=self.verbose)
		if n2p2 and pn2p2 == 1:
			import pynnp
			pn2p2 = pynnp.Prediction()
			pn2p2.log.writeToStdout = False
			pn2p2.setup()
		if n2p2:
			potential, force, tt1, tt2 = self.energy_func(init_state,weights,id_rep,pn2p2,verbose=self.verbose)
		masses_kg = self.masses * prm.amu2kg
		inverse_masses = np.array(1 / masses_kg)
		acceleration = force * inverse_masses[:, np.newaxis]
		accelerations.append(acceleration)
		velocity.append(input_velocity + accelerations[0] * self.time_step * 1e-15)
		#
		# Initialize plumed if used
		if self.use_plumed:
			plumed_calculator = PLM.CalcPlumed(input='plumed.in',timestep=self.time_step,atoms=self.atoms,log='PLUMED.dat')
		#
		# 1/ Compute Langevin constants
		coeff1 = (2 - langevin_friction_coeff * self.time_step) / (2 + langevin_friction_coeff * self.time_step)
		coeff2 = np.sqrt(prm.kB * temperature * self.time_step * (0.5 * langevin_friction_coeff / masses_kg))
		coeff3 = 2 * self.time_step / (2 + langevin_friction_coeff * self.time_step)
		### VVL scheme for one step: ###
		for i in range(steps):
			t1 = time.time()
			if weights == 1.0:
				isprint=True
			else:
				isprint=False
			# 2/ Random eta for Langevin noise
			eta = np.random.normal(0, 1, (len(self.atoms), 3))
			# 3/ Update of velocities at half step
			vel_half_step = velocity[i] + (accelerations[i] * 0.5 * self.time_step * 1e-15) + (coeff2[:, np.newaxis] * eta)
			# 4/ Update of positions at step
			new_position = positions[i] + (coeff3 * vel_half_step * 1e-5)
			new_state = init_state
			new_state.set_positions(new_position)
			positions.append(new_position)
			# 5/ Compute new energy/forces --> new accelerations
			if deepmd:
				potential, force, tt1, tt2 = self.energy_func(new_state,dp,weights,verbose=self.verbose)
			if n2p2:
				potential, force, tt1, tt2 = self.energy_func(new_state,weights,id_rep,pn2p2,verbose=self.verbose)
			# 6/ Constraints are placed here  
			if self.use_plumed:
				if self.verbose:
					print("VERBOSE: PLUMED computation !") 
				curr_pos = new_position
				curr_ener = potential * (prm.hartree2joule ** -1)
				energy_bias, forces_bias = plumed_calculator.compute_bias(curr_pos, i, curr_ener)
				potential += energy_bias[0]  
				potential_energies.append(potential)
				forces_bias *= 1e-10
				force += forces_bias
				potential_energies.append(potential)
				forces.append(force)
			else:
				forces.append(force)
				potential_energies.append(potential)
			acceleration = force * inverse_masses[:, np.newaxis]
			accelerations.append(acceleration)
			# 7/ Update of velocities at step
			velocity.append(coeff1 * vel_half_step + (coeff2[:, np.newaxis] * eta) + (0.5 * accelerations[i + 1] * self.time_step * 1e-15))
			times.append(self.time_step + self.time_step * i)
			# 8/ Compute kinetic energy and virial
			kinetic_energy = 0.5 * np.vdot(velocity[i + 1] * np.array(masses_kg)[:, np.newaxis], velocity[i + 1])
			kinetic_energy_kcal = kinetic_energy * prm.joule2kcal * prm.Na
			kinetic_energies.append(kinetic_energy_kcal)
			total_energy = potential + kinetic_energy_kcal
			total_energies.append(total_energy)
			number_degrees_freedom = 3 * len(self.masses)
			T_inst = (2 * kinetic_energy * prm.Na) / (number_degrees_freedom * prm.gasconst)
			T.append(T_inst)
			t2 = time.time()
			if isprint:
				print_potential = potential
				print_kinetic_energy = float(kinetic_energy_kcal)
				print_total_energy = float(total_energy)
				Prt.print_vvl_step(i+1,print_potential,print_kinetic_energy,print_total_energy,T_inst)
				if self.verbose:
					print('VERBOSE: VVL step {}        Potential energy = {}, Kinetic energy = {}, Total energy = {}, Temperature = {}'.format(i+1,potential,kinetic_energy,total_energy,T_inst))
					print('VERBOSE: Timing (in sec) = {}'.format(t2-t1))
			# 9/ Built new ASE object if it has to be stored within traj
			if (i+1) % stride == 0:
				state = init_state
				state.set_positions(new_position)
				state.energy = float(potential)
				traj.append(state)
				# Flushing trajectory in case of cMD simulations
				if self.dyn_mode == 'cMD':
					Trj.flush(state,'traj',verbose=self.verbose)
					# Updating restart file with positions, velocities and accelerations
					Prt.print_restart_cmd(self.charges,self.chemsymb,new_position,velocity[i + 1],i+1,integration='VVL',verbose=self.verbose)
				if self.dyn_mode == 'hRES':
					name = 'ms.rep'+str(id_rep)
					Trj.flush(state,name,verbose=self.verbose) 
					# Updating restart file with positions and velocities 
					Prt.print_restart_hres(id_rep,state,i+1,integration='VV',v=velocity[i + 1],verbose=self.verbose)
				if isprint and i != 0:
					Prt.print_timing(t1,t2,tt1,tt2,steps,timestep=self.time_step,integration='VVL',verbose=self.verbose)
			if isprint:
				print('STEP {}'.format(i+1))
			t2 = time.time()
		potential_energies = np.array(potential_energies)
		if self.verbose:
			print('VERBOSE: Finishing Velocity Verlet Langevin modified simulation, traj will be printed ...')
		traj = Trj.Trajectory(traj,energy_traj=potential_energies,generation_details=self.simulation_details())
		ener = potential_energies[-1]
		if return_last is True:
			return [traj, state, np.mean(times), velocity[-1]]
		else:
			return traj

