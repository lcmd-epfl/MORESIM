#!/usr/bin/env python3.9
#
"""
Thermostat to apply velocity correction
	=> Bussi thtermostat: G. Bussi and M. Parrinello, "Stochastic Thermostats: Comparison of Local and Global Schemes", Computer Physics Communications, 179, 26-29 (2008)
	=> Nose-Hoover thermostat: Nose, S (1984). "A unified formulation of the constant temperature molecular-dynamics methods". Journal of Chemical Physics. 81 (1): 51-519
	=> function mdrest : remove any translational or rotational kinetic energy of the overall system center of mass
	=> function invert : invert a matrix following the Gauss Jordan Algorithm
"""

import numpy as np
from lib import Params as prm

class Bussi:
	""" A Bussi thermostat class to apply Bussi-Parrinello velocity scaling """

	def __init__(self, timestep=0.5, tautemp=200, temperature=300, verbose=False):
		self.timestep = timestep
		self.tautemp = tautemp
		self.temperature = temperature
		self.verbose = verbose
		if verbose:
			print('VERBOSE: Bussi-Parinello velocity rescaling thermostat has been activated !')

	def run(self, v, nfree, masses_kg, T_inst):
		if self.verbose:
			print('VERBOSE: Rescaling velocities using the Bussi-Parrinello thermostat')
		# Compute internal parameters from Bussi-Parrinello
		c = np.exp(-(self.timestep)/self.tautemp)
		d = (1.0 - c) * (self.temperature/T_inst) / float(nfree)
		mu, sigma = 0, 1
		#r = np.random.normal(mu, sigma, 1)
		r = normal()
		s = 0.0
		for i in range(nfree - 1):
			#si = np.random.normal(mu, sigma, 1)
			si = normal()
			s = s + si*si
		scale = c + (s+r*r)*d + 2.0*r*np.sqrt(c*d)
		scale = np.sqrt(scale)
		print(scale)
		if (r+np.sqrt(c/d)) < 0.0:
			scale = -scale
		# Rescaling velocities here
		for i in range(len(v)):
			for j in range(3):
				v[i][j] = v[i][j] * scale
		# Recompute kinetic energy and instantaneous temperature
		kinetic_energy = 0.5 * np.vdot(v * np.array(masses_kg)[:, np.newaxis], v)
		t_inst = (2 * kinetic_energy * prm.Na * 0.000239) / (float(nfree) * 1.987204259*1e-3)
		# Return all the data
		return v, kinetic_energy, t_inst

class NoseHoover:
	""" A Nose-Hoover thermostat class to apply Nose-Hoover chains """

	def __init__(self, timestep=0.5, temperature=300, verbose=False):
		self.timestep = timestep
		self.temperature = temperature
		self.verbose = verbose
		if verbose:
			print('VERBOSE: Nose-Hoover chains thermostat has been activated !')

	def run(self, v, nfree, cpt, kinetic_energy, vnh, gnh, qnh):
		if self.verbose and cpt == 1:
			print('VERBOSE: Rescaling HALF velocities using the Nose-Hoover thermostat')
		if self.verbose and cpt == 2:
			print('VERBOSE: Rescaling FULL velocities using the Nose-Hoover thermostat')
		ekt = 1.987204259*1e-3 * self.temperature
		kinetic_energy = kinetic_energy * prm.Na * prm.joule2kcal
		nc = 5
		ns = 3
		dtc = (self.timestep*1e-3) / float(nc)
		w = []
		w.append(1.0 / (2.0-2.0**(1.0/3.0)))
		w.append(1.0 - 2.0*w[0])
		w.append(w[0])
		scale = 1.0
		for i in range(nc):
			for j in range(ns):
				dts = w[j] * dtc
				dt2 = 0.5 * dts
				dt4 = 0.25 * dts
				dt8 = 0.125 * dts
				gnh[3] = (qnh[2]*vnh[2]*vnh[2]-ekt) / qnh[3]
				vnh[3] = vnh[3] + gnh[3]*dt4
				gnh[2] = (qnh[1]*vnh[1]*vnh[1]-ekt) / qnh[2]
				expterm = np.exp(-vnh[3]*dt8)
				vnh[2] = expterm * (vnh[2]*expterm+gnh[2]*dt4)
				gnh[1] = (qnh[0]*vnh[0]*vnh[0]-ekt) / qnh[1]
				expterm = np.exp(-vnh[2]*dt8)
				vnh[1] = expterm * (vnh[1]*expterm+gnh[1]*dt4)
				gnh[0] = (2.0*kinetic_energy-float(nfree)*ekt) / qnh[0]
				expterm = np.exp(-vnh[1]*dt8)
				vnh[0] = expterm * (vnh[0]*expterm+gnh[0]*dt4)
				expterm = np.exp(-vnh[0]*dt2)
				scale = scale * expterm
				kinetic_energy = kinetic_energy * expterm * expterm
				gnh[0] = (2.0*kinetic_energy-float(nfree)*ekt) / qnh[0]
				expterm = np.exp(-vnh[1]*dt8)
				vnh[0] = expterm * (vnh[0]*expterm+gnh[0]*dt4)
				gnh[1] = (qnh[0]*vnh[0]*vnh[0]-ekt) / qnh[1]
				expterm = np.exp(-vnh[2]*dt8)
				vnh[1] = expterm * (vnh[1]*expterm+gnh[1]*dt4)
				gnh[2] = (qnh[1]*vnh[1]*vnh[1]-ekt) / qnh[2]
				expterm = np.exp(-vnh[3]*dt8)
				vnh[2] = expterm * (vnh[2]*expterm+gnh[2]*dt4)
				gnh[3] = (qnh[2]*vnh[2]*vnh[2]-ekt) / qnh[3]
				vnh[3] = vnh[3] + gnh[3]*dt4
		for i in range(len(v)):
			for j in range(3):
				v[i][j] = v[i][j] * scale

def normal():
	rsq = 10.0
	while rsq >= 1.0:
		v1 = 2.0*np.random.rand() - 1.0
		v2 = 2.0*np.random.rand() - 1.0
		rsq = v1**2 + v2**2
	factor = np.sqrt(-2.0*np.log(rsq)/rsq)
	store = v1 * factor
	normal = v2 * factor
	return normal

def mdrest(masses_kg, v, position):
	# 0/ Convert the positions from A to m
	p = position.copy()
	p *= 1e-10
	totmass = 0.0
	vtot = [0.0, 0.0, 0.0]
	# 1/ compute linear velocity of the system center of mass
	for i in range(len(masses_kg)):
		totmass += masses_kg[i]
		for j in range(3):
			vtot[j] += (v[i][j] * masses_kg[i])
	# 2/ compute translational kinetic energy of overall system
	etrans = 0.0
	for i in range(3):
		vtot[i] = vtot[i] / totmass
		etrans += vtot[i]**2
	etrans = 0.5 * etrans * totmass
	# 3/ Find COM of the whole system
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for i in range(len(masses_kg)):
		xtot += p[i][0] * masses_kg[i]
		ytot += p[i][1] * masses_kg[i]
		ztot += p[i][2] * masses_kg[i]
	xtot = xtot/totmass
	ytot = ytot/totmass
	ztot = ztot/totmass
	# 4/ compute the angular momentum of the overall system
	mang = [0.0, 0.0, 0.0]
	for i in range(len(masses_kg)):
		mang[0] += (p[i][1]*v[i][2] - p[i][2]*v[i][1]) * masses_kg[i]
		mang[1] += (p[i][2]*v[i][0] - p[i][0]*v[i][2]) * masses_kg[i]
		mang[2] += (p[i][0]*v[i][1] - p[i][1]*v[i][0]) * masses_kg[i]
	mang[0] -= (ytot*vtot[2]-ztot*vtot[1]) * totmass
	mang[1] -= (ztot*vtot[0]-xtot*vtot[2]) * totmass
	mang[2] -= (xtot*vtot[1]-ytot*vtot[0]) * totmass
	# 5/ calculate the moment of inertia tensor
	xx = 0.0
	xy = 0.0
	xz = 0.0
	yy = 0.0
	yz = 0.0
	zz = 0.0
	for i in range(len(masses_kg)):
		xdel = p[i][0] - xtot
		ydel = p[i][1] - ytot
		zdel = p[i][2] - ztot
		xx += (xdel*xdel * masses_kg[i])
		xy += (xdel*ydel * masses_kg[i])
		xz += (xdel*zdel * masses_kg[i]) 
		yy += (ydel*ydel * masses_kg[i])
		yz += (ydel*zdel * masses_kg[i])
		zz += (zdel*zdel * masses_kg[i])
	tensor = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
	tensor[0][0] = yy + zz
	tensor[1][0] = -xy
	tensor[2][0] = -xz
	tensor[0][1] = -xy
	tensor[1][1] = xx + zz
	tensor[2][1] = -yz
	tensor[0][2] = -xz
	tensor[1][2] = -yz
	tensor[2][2] = xx + yy
	# 6/ diagonalize the moment of inertia tensor
	invert(3,tensor)
	#np.linalg.inv(tensor)
	# 7/ compute angular velocity and rotational kinetic energy
	erot = 0.0
	vang = [0.0, 0.0, 0.0]
	for i in range(3):
		for j in range(3):
			vang[i] += tensor[i][j]*mang[j]
		erot += vang[i]*mang[i]
	erot = 0.5 * erot
	# 8/ eliminate any translation of the overall system
	for i in range(len(masses_kg)):
		for j in range(3):
			v[i][j] -= vtot[j]
	# 9/ eliminate any rotation about the system center of mass
	for i in range(len(masses_kg)):
		xdel = p[i][0] - xtot
		ydel = p[i][1] - ytot
		zdel = p[i][2] - ztot
		v[i][0] = v[i][0] - vang[1]*zdel + vang[2]*ydel
		v[i][1] = v[i][1] - vang[2]*xdel + vang[0]*zdel
		v[i][2] = v[i][2] - vang[0]*ydel + vang[1]*xdel

def invert(size, matrix):
	""" Perform matrix inversion via the Gauss Jordan algorithm """
	ipivot = np.zeros(size, dtype=float)
	indxc = np.zeros(size, dtype=float)
	indxr = np.zeros(size, dtype=float)
	for i in range(size):
		big = 0.0
		for j in range(size):
			if ipivot[j] != 1:
				for k in range(3):
					if ipivot[k] == 0:
						if abs(matrix[j][k]) >= big:
							big = abs(matrix[j][k])
							irow = j
							icol = k
		ipivot[icol] += 1
		if irow != icol:
			for j in range(size):
				temp = matrix[irow][j]
				matrix[irow][j] = matrix[icol][j]
				matrix[icol][j] = temps
		indxr[i] = irow
		indxc[i] = icol
		pivot = matrix[icol][icol]
		matrix[icol][icol] = 1.0
		for j in range(size):
			matrix[icol][j] = matrix[icol][j]/pivot
		for j in range(size):
			if j != icol:
				temp = matrix[j][icol]
				matrix[j][icol] = 0.0
				for k in range(size):
					matrix[j][k] -= (matrix[icol][k]*temp)
		for j in range(size-1, 0, -1):
			if indxr[j] != indxc[j]:
				for k in range(size):
					temp = matrix[k][indxr[i]]
					matrix[k][indxr[i]] = matrix[k][indxc[i]]
					matrix[k][indxc[i]] = temp
	
