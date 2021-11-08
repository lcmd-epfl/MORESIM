#!/usr/bin/env python3.9
#
import sys
import os
import numpy as np
import argparse
import psutil
from ase.io import read
from lib import MC, EnergyCalculators
from lib import Reservoir as RSV
from lib import Printing as Prt
from lib import Restart as Rst
from lib import Trajectory as Trj
from lib import VelocityVerlet_Langevin as VVL
from lib import VelocityVerlet as VV
from lib import VelocityVerlet_MTS as VVMTS
from lib import ReplicaExchangeSimulation as RXC

def main(structure, periodic, charge, dynamics, integration, mts, replicas, verbose, baselined, gfn, T, mlcorrection, rsv, szrsv, nexc, nsteps, timestep, stride, langevin, restart, thermostat, plumed, rseed):
	#####################################################################################################################
	# MEMORY USAGE !
	mem = psutil.Process().memory_info().rss * 1e-9
	print("Initial memory allocated at the moment: {} GB".format(mem))
	# Initial parametrization
	try:
		init_state = read(structure)
	except:
		raise('An error occurs during the XYZ reading ! I kill you')
	try:
		atoms = []
		with open('type.dat','r') as f:
			for line in f:
				line = line.strip()
				atoms.append(line)
		atnums = init_state.get_atomic_numbers()
		charges = np.zeros(len(atnums))
		charges[0] = charge
		init_state.set_initial_charges(charges)
	except:
		raise('An error occurs during the atom type reading ! I kill you')
	if periodic:
		try:
			f = open('pbc.dat','r')
			a = f.readline()
			b = a.split()
			init_state.pbc = (True, True, True)
			init_state.set_cell([(float(b[0]), float(b[1]), float(b[2])), (float(b[3]), float(b[4]), float(b[5])), (float(b[6]), float(b[7]), float(b[8]))])
			f.close()
		except:
			raise('An error occurs during the PBC reading ! I kill you')
	if dynamics == 'cMD':
		weights=[1,1]
	if dynamics == 'hRES':
		weights = [np.array([1, 1 - i / (replicas - 1)]) for i in range(replicas)]
	#####################################################################################################################

	#####################################################################################################################
	# Checking procedures for the initialization
	if verbose:
		print('VERBOSE: verbose option has been activated !')
	if dynamics != 'cMD' and dynamics != 'hRES':
		raise('ERROR: dynamics type has to be either cMD or hRES in the python arguments ! I kill you')
	if integration != 'MC' and integration != 'VVL' and integration != 'VV' and integration != 'VV_MTS':
		raise('ERROR: integration mode has to be either MC, VV or VVL ! I kill you')

	if integration == "MC":
		rand_num_gen = MC.GaussianVar(loc=0, var=0.0010)
		mc_moves = [MC.RandomParticleMove(rand_num_gen=rand_num_gen)]

	if baselined != 'DFTB' and baselined != 'XTB':
		raise('ERROR: baselined has to be DFTB ! I kill you')
	if mlcorrection != 'N2P2' and mlcorrection != 'DeepMD' and mlcorrection != 'LKR':
		raise('ERROR: ML correction has to be either n2p2, DeepMD or LKR ! I kill you !')
	if mlcorrection == 'LKR' and integration != 'MC':
		raise('ERROR: ML correction LKR is ONLY compatible with MC integration ! I kill you !')

	if dynamics == 'cMD':
		replicas = 1
		if timestep <= 0:
			raise('Error: timestep has to be positive ! I kill you')
	if dynamics == 'hRES':
		if replicas % 2 != 0:
			raise('ERROR: Number of replicas must be pair ! I kill you')
		if len(weights) != replicas:
			raise('ERROR: Wrong number of temperatures ! I kill you')
	if integration == 'VVL':
		if verbose:
			print('VERBOSE: Velocity Verlet Langevin modified has been detected !')
		if langevin <= 0:
			raise('ERROR: langevin friction coefficient has to be positive ! I kill you !')
	if integration == 'VV':
		if verbose:
			print('VERBOSE: Conventional Velocity Verlet has been detected !')
		if timestep <= 0:
			raise('ERROR: timestep has to be positive ! I kill you')
		if thermostat != 'Bussi' and thermostat != 'Nose-Hoover':
			raise('ERROR: wrong thermostat, has to be either Bussi or Nose-Hoover ! I kill you !')
	if integration == 'VV_MTS':
		if verbose:
			print('VERBOSE: MTS Velocity Verlet has been detected !')
		if timestep <= 0:
			raise('ERROR: timestep has to be positive ! I kill you')
		if thermostat != 'Bussi' and thermostat != 'Nose-Hoover':
			raise('ERROR: wrong thermostat, has to be either Bussi or Nose-Hoover ! I kill you !')
	#####################################################################################################################

	#####################################################################################################################
	# Updating the control file
	Prt.print_init(structure,init_state,periodic,dynamics,integration,replicas,weights,verbose)
	if integration != 'MC':
		Prt.print_dynparams(dynamics,timestep,nsteps,integration,nexc,langevin,stride,thermostat,verbose=verbose)
	#####################################################################################################################

	#####################################################################################################################
	# DFTB baselined parametrization
	# NOTE: This is the only part where the user has to check if all is good for him !
	if baselined == 'DFTB':
		try:
			dftbs = [EnergyCalculators.DftbEnergy(
				atoms=atoms,
				directory='./dftb_comp_{}'.format(i),
				Hamiltonian_ReadInitialCharges='No',
				Hamiltonian_MaxSCCIterations=500,
				Hamiltonian_SCC='Yes',
				Hamiltonian_ThirdOrderFull='Yes',
				Hamiltonian_SCCTolerance=1.0E-008,
				Hamiltonian_Filling='Fermi {\nTemperature = 0.0009500425602573001 \n}',
				Hamiltonian_SlaterKosterFiles_Prefix = '/home/celerse/ML-photoswitchables/3_Dynamics_on_frags/2_Simone/1_PSW-orgcat/dftb-3ob-3-1_files/',
				Hamiltonian_Dispersion_='SlaterKirkwood',
				Hamiltonian_Dispersion_PolarRadiusCharge_='HybridDependentPol',
				Hamiltonian_Dispersion_PolarRadiusCharge_H='{\n CovalentRadius [Angstrom] = 0.4 \n HybridPolarisations [Angstrom^3,Angstrom,] = {\n 0.386 0.396 0.400 0.410 0.410 0.410 3.5 3.5 3.5 3.5 3.5 3.5 0.8 \n }\n }',
				Hamiltonian_Dispersion_PolarRadiusCharge_C='{\n CovalentRadius [Angstrom] = 0.76 \n HybridPolarisations [Angstrom^3,Angstrom,] = {\n 1.382 1.382 1.382 1.064 1.064 1.064 3.8 3.8 3.8 3.8 3.8 3.8 2.50\n }  \n}',
				Hamiltonian_Dispersion_PolarRadiusCharge_S='{\n CovalentRadius [Angstrom] = 1.02 \n HybridPolarisations [Angstrom^3,Angstrom,] = {\n 3.000 3.000 3.000 3.000 3.000 3.000 4.7 4.7 4.7 4.7 4.7 4.7 4.80\n }  \n}',
				Hamiltonian_MaxAngularMomentum_='',
				Hamiltonian_MaxAngularMomentum_H='"s"',
				Hamiltonian_MaxAngularMomentum_C='"p"',
				Hamiltonian_MaxAngularMomentum_S='"d"',
				Hamiltonian_HubbardDerivs_='',
				Hamiltonian_HubbardDerivs_C='-0.1492',
				Hamiltonian_HubbardDerivs_H='-0.1857',
				Hamiltonian_HubbardDervis_S='-0.11',
				Hamiltonian_Charge=0,
				Hamiltonian_KPointsAndWeights='{\n 0 0 0 1.0\n}',
				Analysis_='',
				Analysis_CalculateForces='Yes',
				) for i in range(replicas)]
		except:
			raise('An error occurs during the reading of DFTB parameters ! I kill you')
	# xtb baselined parametrization
	if baselined == 'XTB':
		try:
			dftbs = [EnergyCalculators.XTB(gfn,i,init_state) for i in range(replicas)]
			init_state.set_calculator(dftbs[0])
		except:
			raise('An error occurs during the reading of XTB parameters ! I kill you')

	#####################################################################################################################

	#####################################################################################################################
	# ML correction parametrization
	try:
		if mlcorrection == "LKR":
			import qml
			training_points, alphas, sigma = np.load('./model_lkr.npy', allow_pickle=True)
			if verbose:
				print('VERBOSE: qml library has been imported for the use of SLATM')
			atoms_symb = init_state.get_chemical_symbols()
			slatm = EnergyCalculators.SLATMGenerator(atoms_symb)
			ml = EnergyCalculators.LKR(representation_generator=slatm, training_representations=training_points, alpha_values=alphas, sigma=sigma)
		if mlcorrection == "N2P2":
			ml = EnergyCalculators.N2P2(periodic)
		if mlcorrection == "DeepMD":
			at_types = []
			a = init_state.get_chemical_symbols()
			for i in range(len(a)):
				for j in range(len(atoms)):
					if a[i] == atoms[j]:
						at_types.append(j)
			ml = EnergyCalculators.DeepMD(periodic,at_types)
	except:
		raise('An error occurs during the ML parametrization ! I kill you')
	#####################################################################################################################

	#####################################################################################################################
	# Updating the control file
	Prt.print_check_ml(baselined,mlcorrection,verbose)
	#####################################################################################################################

	#####################################################################################################################
	# Defining mixed potentials for the simulation
	if dynamics == "cMD":
		# 1/ Monte Carlo cMD case
		if integration == 'MC':
			mixed_potentials = [EnergyCalculators.MixedPotential([dftbs[0].energy, ml.energy],forces=False,weights=[1,1])]
		# 2/ Conventional Velocity Verlet and Langevin modified cMD case
		else:
			mixed_potentials = [EnergyCalculators.MixedPotential([dftbs[0].energy_forces, ml.energy_forces],forces=True,weights=[1,1])]
	if dynamics == "hRES":
		# 1/ Monte Carlo hRES case
		if integration == 'MC':
			mixed_potentials = [EnergyCalculators.MixedPotential([dftbs[i].energy, ml.energy],weights=weights) for i in range(replicas)]
		# 2/ Conventional Velocity Verlet and Langevin modified hRES case
		else:
			mixed_potentials = [EnergyCalculators.MixedPotential([dftbs[i].energy_forces, ml.energy_forces],weights=weights) for i in range(replicas)]
	#####################################################################################################################

	#####################################################################################################################
	# check for a plumed activation
	if plumed:
		if verbose:
			print('VERBOSE: Plumed has been activated !')
		if integration == 'MC':
			raise('ERROR: Plumed is not available with MC algorithm ! I kill you !')
		if dynamics != 'cMD':
			raise('ERROR: Plumed is only available with the cMD mode ! I kill you !')
		if os.path.isfile('plumed.in') == False:
			raise('ERROR: plumed.in file does not exist ! I kill you !')
		Prt.print_plumed_init()
	#####################################################################################################################

	#####################################################################################################################
	# check if a restart procedure is activated for cMD ! 
	if restart:
		if verbose:
			print('VERBOSE: Restart procedure has been activated !')
		if dynamics == 'cMD':
			new_init_state = Rst.Restart_cmd(atoms, init_state, integration, periodic=periodic, verbose=verbose)
			init_state = new_init_state
			# WARNING: XTB calculator has to be reload in this restart case
			if baselined == "XTB":
				init_state.calc = dftbs[0]
	#####################################################################################################################

	#####################################################################################################################
	# Defining the dynamics main setup
	# 1/ cMD case
	if dynamics == "cMD":
		# 1/ Monte Carlo cMD case
		if integration == "MC":
			if mlcorrection == "DeepMD":
				cmd_simulation = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[0].energy_dp,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
			if mlcorrection == "LKR":
				cmd_simulation = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[0].energy_lkr,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
			if mlcorrection == "N2P2":
				cmd_simulation = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[0].energy_n2p2,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
		# 2/ Velocity Verlet Langevin modified cMD case
		if integration == "VVL":
			if mlcorrection == "DeepMD":
				cmd_simulation = VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[0].energy_forces_dp,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,use_plumed=plumed,verbose=verbose)
			if mlcorrection == "N2P2":
				cmd_simulation = VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[0].energy_forces_n2p2,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,use_plumed=plumed,verbose=verbose)
		# 3/ Conventional Velocity Verlet cMD case
		if integration == "VV":
			if mlcorrection == "DeepMD":
				cmd_simulation = VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[0].energy_forces_dp,temperature=T,dyn_mode=dynamics,thermostat=thermostat,use_plumed=plumed,verbose=verbose)
			if mlcorrection == "N2P2":
				cmd_simulation = VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[0].energy_forces_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,use_plumed=plumed,verbose=verbose)
		# 4/ MTS Velocity Verlet cMD case
		if integration == "VV_MTS":
			if mlcorrection == "N2P2":
				cmd_simulation = VVMTS.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[0].energy_forces_mts_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,use_plumed=plumed,verbose=verbose)

	# 2/ hRES case
	if dynamics == "hRES":
		# 1/ Monte Carlo hRES case
		if integration == "MC":
			if mlcorrection == "DeepMD":
				rep_ex_simulations = [MC.Simulation(atoms=init_state,energy_func=mixed_potentials[i].energy_hres_dp,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose) for i in range(replicas - 1)]
			if mlcorrection == "LKR":
				rep_ex_simulations = [MC.Simulation(atoms=init_state,energy_func=mixed_potentials[i].energy_lkr,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose) for i in range(replicas - 1)]
			if mlcorrection == "N2P2":
				rep_ex_simulations = [MC.Simulation(atoms=init_state,energy_func=mixed_potentials[i].energy_hres_n2p2,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose) for i in range(replicas - 1)]
		# 2/ Velocity Verlet Langevin modified hRES case
		if integration == 'VVL':
			if mlcorrection == "DeepMD":
				rep_ex_simulations = [VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[i].energy_forces_hres_dp,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,verbose=verbose) for i in range(replicas - 1)]
			if mlcorrection == "N2P2":
				rep_ex_simulations = [VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[i].energy_forces_hres_n2p2,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,verbose=verbose) for i in range(replicas - 1)]
		# 3/ Conventional Velocity Verlet hRES case
		if integration == 'VV':
			if mlcorrection == "DeepMD":
				rep_ex_simulations = [VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[i].energy_forces_hres_dp,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose) for i in range(replicas - 1)]
			if mlcorrection == "N2P2":
				rep_ex_simulations = [VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[i].energy_forces_hres_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose) for i in range(replicas - 1)]
		# 4/ MTS Velocity Verlet cMD case
		if integration == "VV_MTS":
			if mlcorrection == "N2P2":
				rep_ex_simulations = [VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[i].energy_forces_mts_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose) for i in range(replicas - 1)]
	#####################################################################################################################

	#####################################################################################################################
	# Defining the reservoir for hRES simulations
	if dynamics == "hRES":
		try:
			reserv = RSV.Reservoir(structures_folder=rsv,size=szrsv,temperature=T,energy_func=mixed_potentials[-1].energy_reservoir,id_rep=replicas-1)
			rep_ex_simulations.append(reserv)
			if verbose:
				print('VERBOSE: Be aware that reservoir structures coordinates have to be in angstrom and energies in kcal/mol !')
			# Updating the control file
			Prt.print_check_reservoir(rsv,szrsv)
		except:
			print("WARNING: no reservoir was detected ! Thus the last windows will be a full baselined level0")
			if integration == "MC":
				if mlcorrection == "DeepMD":
					reserv = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[-1].energy_hres_dp,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
				if mlcorrection == "LKR":
					reserv = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[-1].energy_lkr,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
				if mlcorrection == "N2P2":
					reserv = MC.Simulation(atoms=init_state,energy_func=mixed_potentials[-1].energy_hres_n2p2,temperature=T,move_list=mc_moves,dyn_mode=dynamics,verbose=verbose)
			if integration == 'VVL':
				if mlcorrection == "DeepMD":
					reserv = VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[-1].energy_forces_hres_dp,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,verbose=verbose)
				if mlcorrection == "N2P2":
					reserv = VVL.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[-1].energy_forces_hres_n2p2,temperature=T,langevin_friction_coeff=langevin,dyn_mode=dynamics,verbose=verbose)
			if integration == 'VV':
				if mlcorrection == "DeepMD":
					reserv = VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[-1].energy_forces_hres_dp,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose)
				if mlcorrection == "N2P2":
					reserv = VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[-1].energy_forces_hres_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose)
			if integration == 'VV_MTS':
				if mlcorrection == "N2P2":
					reserv = VV.Simulation(time_step=timestep,atoms=init_state,energy_func=mixed_potentials[-1].energy_forces_mts_n2p2,temperature=T,dyn_mode=dynamics,thermostat=thermostat,verbose=verbose)
			rep_ex_simulations.append(reserv)

		# Defining Replica Exchange itself and check if restart is activated for hres
		init_states = []
		init_states = [init_state.copy() for i in range(replicas)]
		if restart:
			new_init_states = Rst.Restart_hres(replicas, init_states, integration, periodic=periodic, verbose=verbose)
			init_states = new_init_states
		# WARNING: XTB calculator has to be reload in these cases
		if baselined == "XTB":
			for i in range(replicas):
				init_states[i].set_calculator(dftbs[i])
		if verbose:
			print('VERBOSE: Initializing the {} replicas from the REXC class'.format(replicas))
		if baselined == 'XTB':
			xtb = True
		else:
			xtb = False
		if baselined == 'XTB':
			for i in range(replicas):
				init_states[i].calc = dftbs[i]
		if mlcorrection == 'DeepMD':
			rep_ex = RXC.REXC(num_reps=replicas, simulations=rep_ex_simulations, init_states=init_states, weights=weights, stride=stride, num_processors=replicas, rep_steps=nsteps, integration=integration, xtb=xtb, verbose=verbose)
		if mlcorrection == 'LKR':
			rep_ex = RXC.REXC(num_reps=replicas, simulations=rep_ex_simulations, init_states=init_states, weights=weights, stride=stride, num_processors=replicas, rep_steps=nsteps, integration=integration, xtb=xtb, deepmd=False, lkr=True, n2p2=False, verbose=verbose)
		if mlcorrection == 'N2P2':
			rep_ex = RXC.REXC(num_reps=replicas, simulations=rep_ex_simulations, init_states=init_states, weights=weights, stride=stride, num_processors=replicas, rep_steps=nsteps, integration=integration, xtb=xtb, deepmd=False, lkr=False, n2p2=True, verbose=verbose)
	#####################################################################################################################

	#####################################################################################################################
	# Launching the dynamics
	mem = psutil.Process().memory_info().rss * 1e-9
	print("Memory allocated at the moment 'LAUNCHING THE DYNAMICS': {} GB".format(mem))
	if dynamics == "cMD":
		if integration == "MC":
			Prt.dyn_start(verbose)
			Prt.dyn_mc_label(verbose)
			if mlcorrection == 'DeepMD':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, return_last=False)
			if mlcorrection == 'LKR':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, deepmd=False, return_last=False)
			if mlcorrection == 'N2P2':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, deepmd=False, lkr=False, n2p2=True, return_last=False)
		if integration == "VVL":
			Prt.dyn_start(verbose)
			Prt.dyn_vvl_label(verbose)
			if mlcorrection == 'DeepMD':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, rseed=rseed, return_last=False)
			if mlcorrection == 'N2P2':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, rseed=rseed, deepmd=False, n2p2=True, return_last=False)
		if integration == "VV":
			Prt.dyn_start(verbose)
			Prt.dyn_vvl_label(verbose)
			if mlcorrection == 'DeepMD':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, rseed=rseed, return_last=False)
			if mlcorrection == 'N2P2':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, rseed=rseed, deepmd=False, n2p2=True, return_last=False)
		if integration == "VV_MTS":
			Prt.dyn_start(verbose)
			Prt.dyn_vvl_label(verbose)
			if mlcorrection == 'N2P2':
				traj = cmd_simulation.run(init_state=init_state, steps=nsteps, stride=stride, rseed=rseed, deepmd=False, n2p2=True, return_last=False)
	if dynamics == "hRES":
		if integration == "MC":
			Prt.dyn_start(verbose)
			Prt.dyn_mc_label(verbose)
			traj = rep_ex.run(num_exchanges=nexc, return_last=True, dftbs=dftbs)
		if integration == "VVL":
			Prt.dyn_start(verbose)
			Prt.dyn_vvl_label(verbose)
			traj = rep_ex.run(num_exchanges=nexc, return_last=True, dftbs=dftbs)
		if integration == "VV":
			Prt.dyn_start(verbose) 
			Prt.dyn_vvl_label(verbose) 
			traj = rep_ex.run(num_exchanges=nexc, return_last=True, dftbs=dftbs)
		if integration == "VV_MTS":
			Prt.dyn_start(verbose)
			Prt.dyn_vvl_label(verbose)
			traj = rep_ex.run(num_exchanges=nexc, return_last=True, dftbs=dftbs)
	#####################################################################################################################

	#####################################################################################################################
	#####################################################################################################################
	#####################################################################################################################

#####################################################################################################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=None)
	parser.add_argument('-s', '--structure', type = str, default = 'system.xyz', help = 'The path with the corresponded xyz filename. Default is system.xyz')
	parser.add_argument('-p', '--periodic', type = bool, default = False, help = 'Boolean to specify if the dynamics is periodic or not. Default is False')
	parser.add_argument('-chg', '--charge', type = int, default = 0, help = 'Total charge of the system. Should be an integer. Default is 0')
	parser.add_argument('-dyn', '--dynamics', type = str, default = 'hRES', help = 'The dynamics modes (cMD or hRES). Default is hRES')
	parser.add_argument('-int', '--integration', type = str, default = 'MC', help = 'The integration mode (MC, VV, VVL or VV_MTS). Default is MC')
	parser.add_argument('-mts', '--multitimestep', type = list, default = [1,6], help = 'MTS parametrization. Default is 1:6')
	parser.add_argument('-rep', '--replicas', type = int, default = 4, help = 'The number of replicas in the dynamics. Default is 4')
	parser.add_argument('-verb', '--verbose', type = bool, default = False, help = 'Activate or not Verbose option. Default is False')
	parser.add_argument('-bsnld', '--baselined', type = str, default = 'DFTB', help = "The baselined level. Default is DFTB")
	parser.add_argument('-gfn', '--gfn', type = int, default = '0', help = "GFN functional for XTB baselined (0, 1, 2). Default is 0")
	parser.add_argument('-T', '--temperature', type = float, default = '298.15', help = "The temperature of the simulation in K. Default is 298.15 K")
	parser.add_argument('-ml', '--mlcorrection', type = str, default = 'DeepMD', help = "The ML correction level (LKR, N2P2, DeepMD). Default is DeepMD")
	parser.add_argument('-rsv', '--reservoir', type = str, default = './', help = "The path for the reservoir structures and energies. Default is ./")
	parser.add_argument('-szrsv', '--size_reservoir', type = int, default = 0, help = "The size of the reservoir. Default is 0")
	parser.add_argument('-exc', '--exchange', type = int, default = 100, help = "The number of exchange betweeen windows. Default is 100")
	parser.add_argument('-nstp','--nsteps', type = int, default = 20, help = "The number of steps within onne window. Default is 20")
	parser.add_argument('-ts', '--timestep', type = float, default = 0.5, help = "The timestep in case of VVL. Default is 0.5")
	parser.add_argument('-freq', '--stride', type = int, default = 10, help = "The stride for trajectory file. Default is 10")
	parser.add_argument('-lgv', '--langevin', type = float, default = 0.1, help = "The Langevin factor for VVL dynamics. Default is 10")
	parser.add_argument('-rst', '--restart', type = bool, default = False, help = "Boolean to specify if we retart or not the dynamics. Default is False")
	parser.add_argument('-thr', '--thermostat', type = str, default = 'Nose-Hoover', help = "Thermostat for Velocity Verlet simulation. Default is Nose-Hoover")
	parser.add_argument('-plm', '--plumed', type = str, default = 'False', help = "Plumed activation. Default is False")
	parser.add_argument('-rseed', '--randomseed', type = int, default = 1, help = "Randomseed")
	args = parser.parse_args()
	sys.exit(main(args.structure, args.periodic, args.charge, args.dynamics, args.integration, args.multitimestep, args.replicas, args.verbose, args.baselined, args.gfn, args.temperature, args.mlcorrection, args.reservoir, args.size_reservoir, args.exchange, args.nsteps, args.timestep, args.stride, args.langevin, args.restart, args.thermostat, args.plumed, args.randomseed))

#####################################################################################################################
