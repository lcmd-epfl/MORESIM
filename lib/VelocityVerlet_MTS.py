#!/usr/bin/env python3.9
#
"""
Conventional Velocity Verlet algorithm
Verlet, Loup (1967). "Computer "Experiments" on Classical Fluids. I. Thermodynamical Properties of Lennard-Jones Molecules". Physical Review. 159 (1): 98-103.
Swope, William C.; H. C. Andersen; P. H. Berens; K. R. Wilson (1 January 1982). "A computer simulation method for the calculation of equilibrium constants for the formation of physical clusters of molecules: Application to small water clusters". The Journal of Chemical Physics. 76 (1): 648
"""

import numpy as np
from lib import Maxwell_boltzmann_distribution as MBD
import time
import os
import shutil
import ase
import sys
import psutil
from lib import Trajectory as Trj
from lib import Printing as Prt
from lib import Params as prm
from lib import Thermostat as Thr
os.environ['KMP_WARNINGS'] = '0'

class Simulation:
	""" Class that propagate the atomic positions by the conventional
	Velocity Verlet althorithm for any representation class """

	def __init__(self, time_step, atoms, energy_func, mode='NVT', thermostat='Bussi', temperature=300, dyn_mode='cMD', mts=[1,6], verbose=False):
		"""
		This is small note about the significance and respective units of each variables declared here
		in the Simulation class:
			=> time_step                    = step used to propagate the dynamics           (in fs)
			=> atoms                        = atom types in the dynamics  
			=> energy_force_calculator      = function (in EnergyCalculator.py --> MixedPotentials)
							used to compute potential energy and forces
			=> mode				= Thermodynamics ensemble
			=> temperature                  = targetted temperature through the dynamics    (in K)
			=> dyn_mode                     = dynamics mode (useful for the flush option)
			=> mts				= Multiple Timesteping mode 
			=> verbose                      = debugging mode
		"""
		self.atoms = atoms
		self.time_step = time_step
		self.ase_molecule = ase.Atoms(atoms)
		self.masses = self.ase_molecule.get_masses()
		self.charges = self.ase_molecule.get_atomic_numbers()
		self.chemsymb = self.ase_molecule.get_chemical_symbols()
		self.energy_func = energy_func
		self.mode = mode
		self.temperature = temperature
		self.dyn_mode = dyn_mode
		self.mts = mts
		self.verbose = verbose
		self.beta = (prm.kb * self.temperature) ** - 1
		if thermostat == 'Bussi':
			self.thermostat_name = 'Bussi'
			self.thermostat = Thr.Bussi(time_step,temperature,verbose=verbose)
		elif thermostat == 'Nose-Hoover':
			self.thermostat_name = 'Nose-Hoover'
			self.tautemp = 0.2
			self.thermostat = Thr.NoseHoover(time_step,temperature,verbose=verbose)
		else:
			raise('Wrong value of thermostat ! Should be Bussi or Nose-Hoover ! I kill you')
		if mode == 'NVT':
			if verbose:
				print('VERBOSE: NVT mode detected !')
		if mode == 'NPT':
			if verbose:
				print('VERBOSE: NPT mode detected !')

	def simulation_details(self):
		return vars(self)

	def simulation_type(self):
		return Trj.Trajectory(generation_details=self.simulation_details())

	def run(self, init_state, steps, stride=10, rseed=1, return_last=True, weights=1.0, id_rep=0, deepmd=True, n2p2=False, pn2p2=1, *args):
		np.random.seed()
		"""
		This is small note about how this attribute work.
		It starts from a structures passing in argument and propagate the Velocity Verlet
		Algorithm following the conventional algorithm. Thus, the Bussi thermostat is
		used to ensure the velocity rescaling !
		The routine works like that:
			1/ Initialize velocities from
				a) already existed structures
				b) nothing => Maxwell Boltzmann distribution is thus used
			2/ Compute initial forces and accelerations => update new velocities
			3/ Start the propagation along the desired number of steps:
				a) Compute v1/2 => update new positions then
				b) Compute new forces and accelerations => update velocities
				c) Compute kinetic energy and T inst
				d) Apply the thermostat to rescale velocities according to the
				targetted temperature
				e) Append stuffs and print stuffs !
		WARNING: an important note about the units which are used here:
			=> These are Standard Units in most of the cases !
			=> velocity             = m.s-1
			=> acceleration         = m.s-2
			=> potential energy     = kcal/mol
			=> forces               = J.m-1
			=> kinetic energy       = J
			=> Temperature          = K
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
		if self.verbose:
			print('VERBOSE: Starting Conventional Velocity Verlet dynamics !')
			print('VERBOSE: setupping initial velocities')
		input_velocity = init_state.get_velocities()
		# Check if velocity is coming from the resrevoir (cf null)
		if input_velocity[0][0] == 0.0:
			if self.verbose:
				print('VERBOSE: No initial velocities have been detected')
				if dyn_mode == 'hRES':
					print('VERBOSE: first vx/vy/vz equal to 0 => maybe it could  come from the reservoir')
				print('VERBOSE: => will be reinitialized using the Maxwell Boltzmann distribution !')
			input_velocity = MBD.maxwell_boltzmann_distribution(self.atoms, self.temperature, rseed, isprint)
		else:
			if self.verbose:
				print('VERBOSE: already existed velocities have been detected')
			input_velocity = MBD.maxwell_boltzmann_distribution(self.atoms, self.temperature, rseed, isprint)
		# STEP 0 => Initialization of the dynamics
		thermostat = self.thermostat
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
			from deepmd.infer import DeepPot
			dp = DeepPot('graph.pb') 
			potential, force, tt1, tt2 = self.energy_func(init_state,dp,weights,verbose=self.verbose)
		if n2p2 and pn2p2 == 1:
			import pynnp
			# 1/ Direct N2P2 object
			try:
				os.chdir('DIRECT/')
				pn2p2_d = pynnp.Prediction()
				pn2p2_d.log.writeToStdout = False
				pn2p2_d.setup()
				os.chdir('../')
			except:
				raise ('ERROR: error occurs during the DIRECT N2P2 object creation')
			# 2/ Baselined N2P2 object
			try:
				os.chdir('BASELINED/')
				pn2p2_b = pynnp.Prediction()
				pn2p2_b.log.writeToStdout = False
				pn2p2_b.setup()
				os.chdir('../')
			except:
				raise ('ERROR: error occurs during the BASELINED N2P2 object creation')
		if n2p2:
			potential, force, tt1, tt2 = self.energy_func(init_state,weights,id_rep,pn2p2_b,verbose=self.verbose)
		masses_kg = self.masses * prm.amu2kg
		inverse_masses = np.array(1 / masses_kg)
		acceleration = force * inverse_masses[:, np.newaxis]
		accelerations.append(acceleration)
		velocity.append(input_velocity + accelerations[0] * self.time_step * 1e-15)
		#
		# Initialize arrows for Nose-Hoover case
		if self.thermostat_name == 'Nose-Hoover':
			vnh = [0.0, 0.0, 0.0, 0.0]
			gnh = [0.0, 0.0, 0.0, 0.0]
			qnh = [0.0, 0.0, 0.0, 0.0]
			# set masses for Nose-Hoover thermostat
			ekt = 1.987204259*1e-3 * self.temperature
			qterm = ekt * self.tautemp * self.tautemp
			for i in range(len(qnh)):
				if qnh[i] == 0.0:
					qnh[i] = qterm
			number_degrees_freedom = (3 * len(self.masses)) - 6
			qnh[1] = float(number_degrees_freedom) * qnh[1]
		#
		# MTS SCHEME !
		# => default is defined by self.mts
		# => mts count will be used to switch from direct to baselined
		mts_count = 0
		### VV scheme for one step: ###
		for i in range(steps):
			t1 = time.time()
			# 1/ Find half-step velocities and full-step positions via Verlet recursion
			vel_half_step = velocity[i] + (accelerations[i] * 0.5 * self.time_step * 1e-15)
			new_position = positions[i] + (self.time_step * vel_half_step * 1e-5)
			new_state = init_state
			new_state.set_positions(new_position)
			positions.append(new_position)
			# 2/ Get constraint-corrected positions and half-step velocities
			# Place here future constraints !!
			# 3/ Get the potential energy and atomic forces
			if deepmd:
				if mts_count % self.mts[1] == 0:
					if self.verbose:
						print('VERBOSE: MTS is baselined')
					potential, force, tt1, tt2 = self.energy_func(new_state,dp,weights,verbose=self.verbose)
				else:
					if self.verbose:
						print('VERBOSE: MTS is direct')
					potential, force, tt1, tt2 = self.energy_func(new_state,dp,weights,verbose=self.verbose)
			if n2p2:
				if mts_count % self.mts[1] == 0:
					if self.verbose:
						 print('VERBOSE: MTS is baselined')
					potential, force, tt1, tt2 = self.energy_func(new_state,weights,id_rep,pn2p2_b,direct=False,baselined=True,verbose=self.verbose)
				else:
					if self.verbose:
						print('VERBOSE: MTS is direct')
					potential, force, tt1, tt2 = self.energy_func(new_state,weights,id_rep,pn2p2_d,direct=True,baselined=False,verbose=self.verbose)
			forces.append(force)
			potential_energies.append(potential)
			# 4/ Make half-step temperature and pressure corrections
			# => ONLY FOR NOSE-HOOVER THERMOSTAT
			if self.thermostat_name == 'Nose-Hoover':
				kinetic_energy = 0.5 * np.vdot(vel_half_step * np.array(masses_kg)[:, np.newaxis], vel_half_step)
				number_degrees_freedom = (3 * len(self.masses)) - 6
				T_inst = (2 * kinetic_energy * prm.Na * 0.000239) / (number_degrees_freedom * 1.987204259*1e-3)
				if self.verbose and isprint:
					print('VERBOSE: Half step current kinetic energy : {}'.format(kinetic_energy))
					print('VERBOSE: Half step current temperature : {}'.format(T_inst))
				self.thermostat.run(vel_half_step, number_degrees_freedom, 1, kinetic_energy, vnh, gnh, qnh)
				kinetic_energy = 0.5 * np.vdot(vel_half_step * np.array(masses_kg)[:, np.newaxis], vel_half_step)
				T_inst = (2 * kinetic_energy * prm.Na * 0.000239) / (number_degrees_freedom * 1.987204259*1e-3)
				if self.verbose and isprint:
					print('VERBOSE: Half step new kinetic energy : {}'.format(kinetic_energy))
					print('VERBOSE: Half step new temperature : {}'.format(T_inst))
			# 5/ Get the next accelerations and find the full-step velocities using the Verlet recursion
			acceleration = force * inverse_masses[:, np.newaxis]
			accelerations.append(acceleration)
			vel = vel_half_step + (0.5 * (accelerations[i + 1]) * self.time_step * 1e-15)
			times.append(self.time_step + self.time_step * i)
			# 6/ Find the constraint-corrected full-step velocities
			# Place here future constraints !! 
			# 7/ Compute kinetic energy and make full-step temperature and pressure corrections
			kinetic_energy = 0.5 * np.vdot(vel * np.array(masses_kg)[:, np.newaxis], vel)
			number_degrees_freedom = (3 * len(self.masses)) - 6
			T_inst = (2 * kinetic_energy * prm.Na * 0.000239) / (number_degrees_freedom * 1.987204259*1e-3)
			#########################
			# Thermostat at full step here !
			if self.thermostat_name == 'Bussi':
				new_vel, new_kinetic_energy, new_T_inst = self.thermostat.run(vel, number_degrees_freedom, masses_kg, T_inst)
				T_inst = new_T_inst
				kinetic_energy = new_kinetic_energy
				kinetic_energy_kcal = kinetic_energy * prm.joule2kcal * prm.Na
				kinetic_energies.append(kinetic_energy_kcal)
				total_energy = potential + kinetic_energy_kcal
				total_energies.append(total_energy)
				# Adjust new velocities
				Thr.mdrest(masses_kg, new_vel, new_position)
				velocity.append(new_vel)
			if self.thermostat_name == 'Nose-Hoover':
				self.thermostat.run(vel, number_degrees_freedom, 1, kinetic_energy, vnh, gnh, qnh)
				kinetic_energy = 0.5 * np.vdot(vel * np.array(masses_kg)[:, np.newaxis], vel)
				T_inst = (2 * kinetic_energy * prm.Na * 0.000239) / (number_degrees_freedom * 1.987204259*1e-3)
				kinetic_energy_kcal = kinetic_energy * prm.joule2kcal * prm.Na
				kinetic_energies.append(kinetic_energy_kcal) 
				total_energy = potential + kinetic_energy_kcal
				total_energies.append(total_energy)
				# Adjust new velocities
				Thr.mdrest(masses_kg, vel, new_position)
				velocity.append(vel)
			#########################
			t2 = time.time() 
			if isprint:
				if deepmd:
					print_potential = potential
					print_kinetic_energy = float(kinetic_energy_kcal)
					print_total_energy = float(total_energy) 
					Prt.print_vvl_step(i+1,print_potential,print_kinetic_energy,print_total_energy,T_inst)
				if n2p2:
					print_potential = potential
					print_kinetic_energy = float(kinetic_energy_kcal)
					print_total_energy = float(total_energy)
					Prt.print_vvl_step(i+1,potential,print_kinetic_energy,print_total_energy,T_inst)
				if self.verbose:
					print('VERBOSE: VVL step {}        Potential energy = {}, Kinetic energy = {}, Total energy = {}'.format(i+1,potential,kinetic_energy,total_energy)) 
					print('VERBOSE: Timing (in sec) = {}'.format(t2-t1)) 
			# 8/ Built new ASE object if it has to be stored withinn traj
			if (i+1) % stride == 0:
				state = init_state
				init_state.set_positions(new_position)
				state.energy = float(potential)
				traj.append(state)
				# Flushing trajectory in case of cMD simulations 
				if self.dyn_mode == 'cMD':
					Trj.flush(state,'traj',verbose=self.verbose) 
					# Updating restart file with positions, velocities and accelerations
					Prt.print_restart_cmd(self.charges,self.chemsymb,new_position,velocity[i + 1],i+1,integration='VV',verbose=self.verbose)
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
			# 9/ Increment MTS counter
			mts_count += 1
		potential_energies = np.array(potential_energies)
		if self.verbose:
			print('VERBOSE: Finishing Conventional Velocity Verlet simulation, traj will be printed ...')
		traj2 = Trj.Trajectory(traj,energy_traj=potential_energies,generation_details=self.simulation_details()) 
		del traj
		ener = potential_energies[-1] 
		if return_last is True: 
			return [traj2, state, np.mean(times), velocity[-1]]  
		else:
			return traj2 

