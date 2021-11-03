#!/usr/bin/env python3.9
#
"""
Maxwell boltzmann distribution for velocities
	=> are in m.s-1 (international units)
"""

import numpy as np
import ase
from lib import Params as prm
from lib import Printing as Prt
from scipy.stats import maxwell

def maxwell_boltzmann_distribution(atoms, temperature, rseed, isprint=False):
	input_velocities = np.zeros(shape=[len(atoms), 3])
	T = temperature
	mass = [1e-3 * M / prm.Na for M in ase.Atoms(atoms).get_masses()]
	standard_deviation = [np.sqrt((prm.kB * T) / m) for m in mass]
	for i in range(len(standard_deviation)):
		for j in range(3):
			if rseed == 1:
				a = np.random.randint(1000000, size=1)
				b = int(a)
				np.random.seed(seed=b)
				if i == 0 and j == 0 and isprint:
					Prt.print_randomseed(b)
			else:
				np.random.seed(seed=rseed)
				if i == 0 and j == 0 and isprint:
					Prt.print_randomseed(rseed)
			input_velocities[i][j] = np.random.normal(loc=0, scale=standard_deviation[i], size=[1])
	return input_velocities

def maxwell_boltzmann_distribution_bis(atoms, temperature):
	input_velocities = np.zeros(shape=[len(atoms), 3])
	T = temperature
	mass = [1e-3 * M / prm.Na for M in ase.Atoms(atoms).get_masses()]
	standard_deviation = [np.sqrt((prm.kB * T) / m) for m in mass]
	for i in range(len(standard_deviation)):
		for j in range(3):
			input_velocities[i][j] = maxwell.rvs(loc=0.0, scale=standard_deviation[i], size=[1])
	return input_velocities
