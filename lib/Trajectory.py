#!/usr/bin/env python3.9
#
"""
Trajectory class to store the traj
"""

import ase.io as aio
import numpy as np
import ase
from lib import Printing as Prt

class Trajectory:
	"""A class to store trajectories"""

	def __init__(self, state_traj=None, energy_traj=None, generation_details=None, flush_prefix=None, verbose=False):
		if state_traj is None:
			state_traj = []
		self.state_traj = state_traj
		if energy_traj is None:
			energy_traj = [state.energy for state in state_traj]
		self.energy_traj = energy_traj
		if generation_details is None:
			generation_details = {}
		self.generation_details = generation_details
		if verbose:
			print('VERBOSE: Well initialized Trajectory class for printing !')

	def extend(self, traj):
		if type(traj) is not type(self):
			raise ValueError('The input is not a trajectory')
		self.state_traj.extend(traj.state_traj)
		self.energy_traj.extend(traj.energy_traj)
		for key in self.generation_details.keys():
			if isinstance(self.generation_details[key], list):
				self.generation_details[key].extend(traj.generation_details[key])

	def flush(self, rep, flush_prefix, step, integration, v, verbose=False):
		if verbose:
			print('VERBOSE: Flusing trajectory ...')
		for struct in self.state_traj:
			aio.write('{}_structures.xyz'.format(flush_prefix),ase.Atoms(struct.get_chemical_symbols(),positions=struct.positions), append=True)
		self.__init__(generation_details=self.generation_details,flush_prefix=flush_prefix)
		if verbose:
			print('VERBOSE: Flushing trajectory successfully !')
		# Updating restart file here for the corresponding replica !
		Prt.print_restart_hres(rep, struct, step, integration, v, verbose=verbose)

def flush(state, flush_prefix, verbose=False):
	aio.write('{}_structures.xyz'.format(flush_prefix),ase.Atoms(state.get_chemical_symbols(),positions=state.positions), append=True)
