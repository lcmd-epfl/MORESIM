#!/usr/bin/env python3.9
#
"""
Functions to ensure proper writing during the processes
"""

import os
import ase
from lib import Params as prm

def print_init(structure,init_state,periodic,dynamics,integration,replicas,weights,verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for checking structure')
	control_file = open('control_file.dat','w')
	with open('lib/MORESIM_LOGO.dat','r') as f:
		for line in f:
			control_file.write(line)
	control_file.write("NOTE: This a control file to check what you're doing on ! \n\n")
	control_file.write("################################################### \n")
	control_file.write("   directory of the structure : {} \n".format(os.getcwd()))
	control_file.write("   name of the file : {} \n".format(structure))
	control_file.write("   number of atoms : {} \n".format(len(init_state.get_positions())))
	control_file.write("   PBC : {} \n".format(periodic))
	if periodic:
		control_file.write("   PBC cell : {} \n".format(init_state.get_cell()))
	control_file.write("################################################### \n")
	control_file.write("\n")
	control_file.write("################################################### \n")
	control_file.write("Checking of the dynamics mode : \n")
	control_file.write("   dynamics mode : {} \n".format(dynamics))
	control_file.write("   integration mode : {} \n".format(integration))
	if integration == 'MC':
		control_file.write("      => MONTE-CARLO integration algorithm ! \n")
	if integration == 'VVL':
		control_file.write("      => VELOCITY VERLET LANGEVIN MODIFIED algorithm \n")
	if integration == 'VV':
		control_file.write("      => VELOCITY VERLET algorithm \n")
	if dynamics == "hRES":
		control_file.write("   number of replicas : {} \n".format(replicas))
		control_file.write("   weights values : {} \n".format(weights))
	control_file.write("################################################### \n")
	if verbose:
		print('VERBOSE: closing of the control_file.dat for checking structure')
	control_file.close()

def print_dynparams(dynamics,timestep,num_steps,integration,exchanges,langevin,stride,thermostat,verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for checking dynamics parameters')
	control_file = open('control_file.dat', 'a')
	control_file.write("\n################################################### \n")
	control_file.write("Dynamics parameters have been detected ! \n")
	if integration == "VVL":
		integrator = "Velocity Verlet Langevin modified"
	if integration == "VV":
		integrator = "Conventional Velocity Verlet"
	if integration == "VV_MTS":
		integrator = "MTS Velocity Verlet"
	control_file.write("   Integrator : {} \n".format(integrator))
	control_file.write("   Timestep   : {} fs \n".format(timestep))
	control_file.write("   Stride     : each {} steps \n".format(stride))
	if integration == "VVL":
		control_file.write('   Langevin friction coefficient : {} fs-1 \n'.format(langevin))
	if integration == "VV":
		control_file.write('   Thermostat : {} \n'.format(thermostat))
	if dynamics == "cMD":
		control_file.write("   Number of steps : {} \n".format(num_steps))
	if dynamics == "hRES":
		control_file.write("   Number of steps between each exchanges : {} \n".format(num_steps))  
		control_file.write("   Number of exchanges : {} \n".format(exchanges))
	control_file.write("################################################### \n")
	control_file.close() 
	if verbose:
		print('VERBOSE: closing of the control_file.dat for checking dynamics parameters')

def print_check_ml(baselined,mlcorrection,verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for checking ML parametrization') 
	control_file = open('control_file.dat','a') 
	control_file.write("\n################################################### \n")  
	control_file.write("Checking of the baselined and the ML correction \n")
	control_file.write("selected : \n")
	control_file.write("   selected baselined : {} \n".format(baselined))   
	control_file.write("   selected ML correction : {} \n".format(mlcorrection))
	control_file.write("################################################### \n") 
	if verbose: 
		print('VERBOSE: closing of the control_file.dat for checking ML parametrization') 
	control_file.close() 

def print_check_reservoir(rsv,szrsv,verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for reservoir update')  
	control_file = open('control_file.dat','a') 
	control_file.write("\n################################################### \n")
	control_file.write("Checking of the reservoir : \n") 
	control_file.write("   path of the energy file : {} \n".format(rsv)) 
	control_file.write("   number of structures and energies : {} \n".format(szrsv))
	control_file.write("################################################### \n")   
	if verbose:
		print('VERBOSE: closing of the control_file.dat for reservoir update')
	control_file.close()  

def print_init_restart_cmd(verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for restarting')
	control_file = open('control_file.dat','a')  
	control_file.write("\n################################################### \n") 
	control_file.write("Restarting procedure has been activated for cMD ! \n")
	control_file.write("################################################### \n")
	if verbose:
		print('VERBOSE: closing of the control_file.dat for restarting')
	control_file.close()  

def print_init_restart_hres(verbose=False):
	if verbose:
		print('VERBOSE: opening of the control_file.dat for restarting')
	control_file = open('control_file.dat','a')
	control_file.write("\n################################################### \n")
	control_file.write("Restarting procedure has been activated for hRES ! \n")
	control_file.write("################################################### \n")
	if verbose:
		print('VERBOSE: closing of the control_file.dat for restarting')
	control_file.close()

def print_randomseed(rseed):
	control_file = open('control_file.dat','a')
	control_file.write("\nRandomseed : {} \n".format(rseed))
	control_file.close()

def dyn_start(verbose=False):
	if verbose: 
		print('VERBOSE: entering in the MC.simulation.run attribut to launch cMD simulation at the MC level')
		print('VERBOSE: opening of the control_file.dat for dynamics starting')
	control_file = open('control_file.dat','a')  
	control_file.write("\n \n \n################################################### \n")
	control_file.write("################################################### \n") 
	control_file.write("############### DYNAMICS STARTS NOW ############### \n") 
	control_file.write("################################################### \n")   
	control_file.write("################################################### \n") 
	control_file.write("\n")  
	control_file.close() 
	if verbose: 
		print('VERBOSE: closing of the control_file.dat for dynamics starting')

def dyn_mc_label(verbose=False):
	if verbose:
		print('VERBOSE: printing of the MC label for each arrow')
	control_file = open('control_file.dat','a')
	control_file.write("STEP \t MC ENERGY (kcal/mol) \t Accepted \n \n".expandtabs(30))
	control_file.close()

def dyn_vvl_label(verbose=False):
	if verbose:
		print('VERBOSE: printing of the VVL label for each arrow')
	control_file = open('control_file.dat','a')
	control_file.write(("STEP \t POTENTIAL ENERGY (kcal/mol) \t KINETIC ENERGY (kcal/mol) \t TOTAL ENERGY (kcal/mol) \t CURRENT TEMPERATURE \n\n").expandtabs(30))
	control_file.close()

def print_mc_step(i,ener,accepted,verbose=False,decomp=False):
	if verbose:
		print('VERBOSE: printing MC step and energy')
	if verbose and decomp:
		print('VERBOSE: energy will be decomposed into the baselined and ML corrected contributions')
	control_file = open('control_file.dat','a')
	control_file.write(("{0:6d} \t {1:12.5f} \t {2:s} \n".format(i, ener, accepted)).expandtabs(30))
	control_file.close()

def print_vvl_step(i,pot,kin,tot,T_inst,verbose=False,decomp=False):
	if verbose:
		print('VERBOSE: printing VVL step and energies')
	if verbose and decomp:
		print('VERBOSE: energy will be decomposed into the baselined and ML corrected contributions')
	control_file = open('control_file.dat','a')
	control_file.write(("{0:6d} \t  {1:12.5f}  \t {2:12.5f} \t {3:12.5f} \t {4:12.5f} \n".format(i, pot, kin, tot, T_inst)).expandtabs(30))
	control_file.close()

def print_timing(t1,t2,tt1,tt2,steps,timestep=0.5,integration='MC',verbose=False):
	if verbose:
		print('VERBOSE: printing timing')
	control_file = open('control_file.dat','a')
	sec_in_day = float(3600*24)
	num_mc_steps = int(sec_in_day/(t2-t1))
	num_step = int(sec_in_day/(t2-t1))*timestep/1000.0
	control_file.write('TIMING: Time per step             (sec): {} \n'.format(t2-t1))
	control_file.write('TIMING: Time for baselined        (sec): {} \n'.format(tt1))
	control_file.write('TIMING: Time for ML corr          (sec): {} \n'.format(tt2))
	if integration == 'MC':
		control_file.write('TIMING: Estimated # of MC steps per day: {} \n'.format(num_mc_steps))
	if integration == 'VVL':
		control_file.write('TIMING: Estimated ps/day: {} \n'.format(num_step))
	control_file.close()

def print_hres_iniexc(i):
	control_file = open('control_file.dat','a')
	control_file.write('\n')
	control_file.write('###### Exchange {} ##### \n'.format(i+1))
	control_file.close()

def print_hres_exc1(old_energies,new_energies,verbose=False):
	if verbose:
		print('VERBOSE: printing hRES exchange informations')
	control_file = open('control_file.dat','a')
	old_energies_bis = []
	new_energies_bis = []
	for i in range(len(old_energies)):
		old_energies_bis.append(old_energies[i])
		new_energies_bis.append(new_energies[i])
	control_file.write('Older energies       : {} \n'.format(old_energies_bis))
	control_file.write('New current energies : {} \n'.format(new_energies_bis))
	control_file.write('\n')
	control_file.close()

def print_hres_exc2(pair,prob,accepted,verbose=False):
	if verbose:
		print('VERBOSE: printing hRES exchange informations')
	control_file = open('control_file.dat','a')
	if accepted:
		acc = '=> Exchange accepeted !'
	else:
		acc = '=> Exchange refused !'
	control_file.write('{} \t {} (ref: 1.0) \t {} \n'.format(pair,prob,acc).expandtabs(5))
	control_file.close()

def print_hres_probaexc(proba,verbose=False):
	if verbose:  
		print('VERBOSE: printing hRES global exchange probabilities')
	control_file = open('control_file.dat','a')
	control_file.write("hRES global probability exchanges: \n")
	[control_file.write('{0}: {1:.3}\n'.format(x, y)) for x, y in proba.items()]
	control_file.write("\n")
	control_file.close()  

def print_mc_newsteps():
	control_file = open('control_file.dat','a')
	control_file.write('\nSTEP \t MC ENERGY (kcal/mol) \t Accepted \n \n'.expandtabs(30))
	control_file.close()

def print_restart_cmd(a_type, symb, p, v, step, integration, verbose=False):
	if verbose:
		print('VERBOSE: printing restart file')
	restart_file = open('restart.dat','w')
	restart_file.write("Restart file from MORESIM containing positions, velocities and accelerations from structure at step {}\n".format(step))
	# 1/ Positions
	restart_file.write("\n")
	restart_file.write("NUMBER OF ATOMS = {} \n".format(len(p)))
	restart_file.write("\n")
	restart_file.write("POSITIONS FOR EACH ATOM (A) \n")
	for i in range(len(p)):
		restart_file.write("{0:5d} \t {1:3d} \t {2:s} \t {3:12.5f} \t {4:12.5f} \t {5:12.5f} \n".format(i+1,a_type[i],symb[i],p[i][0],p[i][1],p[i][2]).expandtabs(5))
	# 2/ Velocities
	if integration != 'MC':
		restart_file.write("\n")
		restart_file.write("VELOCITIES PER ATOM (m.s-1) \n")
		for i in range(len(p)):
			restart_file.write("{0:12.5f} \t {1:12.5f} \t {2:12.5f} \n".format(v[i][0],v[i][1],v[i][2]).expandtabs(5))
	restart_file.close()

def print_restart_hres(rep, struct, step, integration, v, verbose=False):
	if verbose:
		print('VERBOSE: printing restart file for replica {}'.format(rep))
	restart_file = open('restart_{}.dat'.format(rep),'w')
	restart_file.write("Restart file from MORESIM specific to hRES containing positions, velocities and accelerations from structure at exchange {}\n".format(step))
	restart_file.write("REPLICA {} \n".format(rep))
	# 1/ Positions
	a_type = struct.get_atomic_numbers()
	symb = struct.get_chemical_symbols()
	p = struct.get_positions()
	restart_file.write("\n")
	restart_file.write("NUMBER OF ATOMS = {} \n".format(len(p)))
	restart_file.write("\n")
	restart_file.write("POSITIONS FOR EACH ATOM (A) \n")
	for i in range(len(p)):
		restart_file.write("{0:5d} \t {1:3d} \t {2:s} \t {3:12.5f} \t {4:12.5f} \t {5:12.5f} \n".format(i+1,a_type[i],symb[i],p[i][0],p[i][1],p[i][2]).expandtabs(5))
	# 2/ Velocities
	restart_file.write("\n")
	if integration != 'MC':
		restart_file.write("VELOCITIES PER ATOM (m.s-1) \n")
		for i in range(len(p)):
			restart_file.write("{0:12.5f} \t {1:12.5f} \t {2:12.5f} \n".format(v[i][0],v[i][1],v[i][2]).expandtabs(5))
	restart_file.close()

def print_n2p2data(struc, id_rep, periodic, verbose=False):
	if verbose:
		print('VERBOSE: printing input.data file for the n2p2 procedure')
	coords = struc.get_positions()
	if periodic:
		cell = struc.get_cell()
	file_name = "input.data."+str(id_rep)
	f = open(file_name, 'w')
	# Write input.data file in the n2p2 style => RuNNer format style !
	f.write("begin \n")
	f.write("comment \n")
	if periodic:
		f.write("lattice {} {} {} \n".format(cell[0][0]*(1.0/prm.bohr2angst), cell[0][1]*(1.0/prm.bohr2angst), cell[0][2]*(1.0/prm.bohr2angst)))
		f.write("lattice {} {} {} \n".format(cell[1][0]*(1.0/prm.bohr2angst), cell[1][1]*(1.0/prm.bohr2angst), cell[1][2]*(1.0/prm.bohr2angst)))
		f.write("lattice {} {} {} \n".format(cell[2][0]*(1.0/prm.bohr2angst), cell[2][1]*(1.0/prm.bohr2angst), cell[2][2]*(1.0/prm.bohr2angst)))
	at_symb = struc.get_chemical_symbols()
	for i in range(len(coords)):
		f.write("atom {} {} {} {} 0 0 0.0 0.0 0.0 \n".format(coords[i][0]*(1.0/prm.bohr2angst),coords[i][1]*(1.0/prm.bohr2angst),coords[i][2]*(1.0/prm.bohr2angst),at_symb[i]))
	f.write("energy 0.0 \n")
	f.write("charge 0.0 \n")
	f.write("end \n")
	f.close()

def print_end(verbose=False):
	if verbose:
		print('VERBOSE: printing ending of the simulation')
	control_file = open('control_file.dat','a')
	control_file.write("\n\n")
	control_file.write("NOTE: End of the simulation \n")
	control_file.write("NOTE: Simulation ended normally \n")
	control_file.write("NOTE: See you soon !")
	control_file.close()
