#!/usr/bin/env python3.9
#
"""
Restart initialization at the beginning of the dynamics !
	=> Restart_cmd
	=> Restart_hres
"""

import numpy as np
import ase
from lib import Printing as Prt

def Restart_cmd(atoms, init_state, integration, periodic=False, verbose=False):
	try:
		f = open('restart.dat','r')
		for i in range(2):
			f.readline()
		a = f.readline()
		b = a.split()
		num_at = int(b[4])
		for i in range(2):
			f.readline()
		new_p = np.zeros((num_at,3), dtype=float)
		for i in range(num_at):
			a = f.readline() 
			b = a.split()
			for j in range(3):
				new_p[i][j] = float(b[j+3])
		if integration != 'MC':
			for i in range(2):
				f.readline()
			vel = np.zeros((num_at,3), dtype=float)
			for i in range(num_at):
				a = f.readline() 
				b = a.split()  
				for j in range(3): 
					vel[i][j] = float(b[j]) 
		f.close()
	except:
		raise("An error occurs during the reading of restart.dat file ! I kill you !")
	symb = init_state.get_chemical_symbols()
	pbcc = init_state.get_pbc()
	celll = init_state.get_cell()
	new_state = ase.Atoms(symb, new_p)
	if periodic:
		new_state.pbc = pbcc
		new_state.cell = celll
	if integration != 'MC':
		new_state.set_velocities(vel)
	Prt.print_init_restart_cmd()
	return new_state

def Restart_hres(num_reps, init_states, integration, periodic=False, verbose=False):
	try:
		states = []
		for i in range(num_reps):
			f = open('restart_{}.dat'.format(i),'r')
			for i in range(3):
				f.readline()
			a = f.readline()
			b = a.split()
			num_at = int(b[4])
			for j in range(2):
				f.readline()
			new_p = np.zeros((num_at,3), dtype=float)
			for j in range(num_at):
				a = f.readline()
				b = a.split()
				for k in range(3):
					new_p[j][k] = float(b[k+3])
			f.close()
			states.append(new_p)
	except:
		raise("An error occurs during the reading of restart.dat file ! I kill you !")
	symb = init_states[0].get_chemical_symbols()
	pbcc = init_states[0].get_pbc()
	celll = init_states[0].get_cell()
	new_states = []
	for i in range(num_reps):
		new_state = ase.Atoms(symb, states[i])
		if periodic:
			new_state.pbc = pbcc
			new_state.cell = celll
		new_states.append(new_state)
	Prt.print_init_restart_hres()
	return new_states

