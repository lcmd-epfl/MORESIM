#!/usr/bin/env python3.9
#
"""
This module includes different energy functions
to compute energy and forces:
	=> class DftbEnergy
	=> class XTB
	=> class LKR
	=> class N2P2
	=> class DeepMD
	=> class MixedPotential
"""

import ase
import ase.calculators.dftb as adftb
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
from ase.constraints import FixAtoms
from ase.units import Bohr, Hartree
import sys
import time
import numpy as np
import os
import cpuinfo
import psutil
from lib import Params as prm
from lib import Printing as Prt

class DftbEnergy:
	"""docstring for dftb"""

	def __init__(self, atoms, directory, **kwargs):
		self.dftb_kwargs = kwargs
		self.atoms = ase.Atoms(atoms)
		self.directory = directory
		if 'factor' in kwargs:
			self.factor = kwargs['factor']
		else:
			self.factor = 1
		self.calc = adftb.Dftb(**kwargs)
		self.calc.directory = directory	

	def energy(self, state):
		self.calc.calculate(state)
		energy = self.calc.results['energy']
		return energy * prm.ev2kcalmol * self.factor

	def force(self, state):
		self.calc.calculate(state)
		path = os.path.abspath(os.getcwd())
		force = self.calc.read_forces()
		return force * prm.hartree2kcal * (1.0/(prm.angst2m)) * self.factor

	def energy_forces(self, state):
		self.calc.calculate(state)
		energy = self.calc.results['energy']
		force = self.calc.read_forces()
		return energy * prm.ev2kcalmol * self.factor, force * prm.ev2joule * (1.0/(prm.angst2m)) * self.factor

class XTB(FileIOCalculator):
	"""docstring for xtb"""

	implemented_properties = ["energy", "forces"]
	discard_results_on_any_change = True

	def __init__(self, gfn, id_rep, atoms, restart=None, ignore_bad_restart_file=FileIOCalculator._deprecated, verb=0):
		xtbpath = os.environ.get("XTBPATH")
		if xtbpath is None:
			mess = "The XTBPATH variable is not defined."
			raise ValueError(mess)
		self.gfn = gfn
		self.verb = verb
		label = "xtbcalc"+str(id_rep)+"/"
		self.label = label
		command="/home/afabrizi/MyQuantumChemSoftwares/xtb/bin/xtb --gfn " + str(gfn) + " --grad xtb.tmol > xtb.log"
		FileIOCalculator.__init__(self, restart=restart, ignore_bad_restart_file=ignore_bad_restart_file, label=label, atoms=atoms, command=command)
		#self.clean()

	def write_charge_mult(self):
		charge = int(self.atoms.get_initial_charges().sum().round())
		uhf = int(self.atoms.get_initial_magnetic_moments().sum().round())
		with open(self.label + ".CHRG", "w") as f1:
			f1.write(str(charge))
		with open(self.label + ".UHF", "w") as f2:
			f2.write(str(uhf))

	def write_input(self, atoms, properties=None, system_changes=None):
		if not (np.all(atoms.pbc) or not np.any(atoms.pbc)):
			raise RuntimeError("PBC must be all true or all false")
		FileIOCalculator.write_input(self, atoms, properties, system_changes)
		filename = self.label + "xtb.tmol"
		coord = atoms.get_positions()
		symbols = atoms.get_chemical_symbols()
		self.write_charge_mult()
		with open(filename, "w") as fd:
			fix_indices = set()
			if self.atoms.constraints:
				for constr in atoms.constraints:
					if isinstance(constr, FixAtoms):
						fix_indices.update(constr.get_indices())
			fix_str = []
			for i in range(len(atoms)):
				if i in fix_indices:
					fix_str.append("f")
				else:
					fix_str.append("")
			if np.all(atoms.pbc):
				fd.write("$periodic 3\n")
				fd.write("$cell angs\n")
				fd.write(" {0} {1} {2}  ".format(atoms.get_cell_lengths_and_angles()[0],atoms.get_cell_lengths_and_angles()[1],atoms.get_cell_lengths_and_angles()[2]))
				fd.write(" {0} {1} {2}\n".format(atoms.get_cell_lengths_and_angles()[3],atoms.get_cell_lengths_and_angles()[4],atoms.get_cell_lengths_and_angles()[5]))
			fd.write("$coord angs\n")
			for (x, y, z), s, fix in zip(coord, symbols, fix_str):
				fd.write("%20.14f  %20.14f  %20.14f      %2s  %2s \n"% (x, y, z, s.lower(), fix))
			fd.write("$end\n")

	def read_results(self, verb=None):
		if verb is None:
			verb = self.verb
		self.read_energy(self.directory, self.label, verb)
		self.read_forces(self.atoms, self.directory, self.label, verb)

	def read_forces(self, atoms, directory, label, verb=0):
		gradient_file = label + "gradient"
		natoms = len(atoms)
		if not os.path.exists(gradient_file):
			raise ReadError("The gradient file does not exist.")
		if os.path.isfile(gradient_file):
			with open(gradient_file, "r") as fd:
				lines = fd.readlines()
		lines = lines[-1 - natoms : -1]
		assert len(lines) == natoms
		self.results["forces"] = np.zeros((natoms, 3), float)
		for i in range(natoms):
			line = [s for s in lines[i].strip().split(" ") if len(s) > 0]
			f = -np.array([float(x.replace("D", "E")) for x in line[0:3]])
			self.results["forces"][i, :] = f * (Hartree / Bohr)
		os.remove(gradient_file)
		if verb > 1:
			print("Forces of the last atom:\n {0}".format(f * (Hartree / Bohr)))
	
	def read_energy(self, directory, label, verb=0):
		energy_file = label + "energy"
		if not os.path.exists(energy_file):
			raise ReadError("The energy file does not exist.")
		if os.path.isfile(energy_file):
			with open(energy_file, "r") as fd:
				lines = fd.readlines()
		for i in range(len(lines)):
			if lines[i].startswith("$end"):
				escf = (float(lines[i - 1].strip().split(" ")[3].replace("D", "E"))* Hartree)
				self.results["energy"] = escf
				os.remove(energy_file)
				break
		if verb > 0:
			print("Energy:\n {0}".format(escf))

	def finished_successfully(self, verb=0):
		finished = False
		mess = ""
		for line in open(self.label + "xtb.log", "r"):
			if line.rfind("finished run") > -1:
				finished = True
			if line.contains("error"):
				mess += line
		if verb > 0:
			if finished:
				print("Normal termination.")
			else:
				print("Error termination.")
		return finished, mess

	def clean(self):
		files_to_clean = [self.label + "xtb.log",self.label + "energy",self.label + "gradient",self.label + "xtbrestart",self.label + "charges"]
		for f in files_to_clean:
			try:
				os.remove(f)
			except OSError:
				pass

	def energy(self, atoms):
		e = atoms.get_potential_energy()
		return e * prm.hartree2kcal

	def energy_forces(self, atoms):
		e = atoms.get_potential_energy()
		f = atoms.get_forces()[:][:]
		return e * prm.ev2kcalmol, f * prm.ev2joule * (1.0/(prm.angst2m))

class LKR:
	"""docstring for LKR"""

	def __init__(self, representation_generator,training_representations, alpha_values, sigma, baseline=None, delta_scale=1):
		self.baseline = baseline
		self.representation_generator = representation_generator
		self.alpha_values = alpha_values
		self.sigma = sigma
		self.training_representations = training_representations
		self.delta_scale = delta_scale

	def energy(self, structure):
		delta_e = [0]
		if self.baseline is not None:
			ener = self.baseline(structure)
		else:
			ener = 0
		x = self.representation_generator.generate(structure.positions)
		k_vec = kernel(np.reshape(x, [1, len(x)]), self.training_representations, self.sigma)
		delta_e = self.delta_scale * np.dot(k_vec, self.alpha_values)
		return ener + delta_e[0]

def kernel(x, Xt, sigma):
	import qml
	return qml.kernels.gaussian_kernel(x, Xt, sigma)

class SLATMGenerator:
	"""docstring for SLATMGenerator"""

	def __init__(self, atoms):
		self.atoms = atoms
		self.atomic_numbers = ase.Atoms(symbols=atoms).get_atomic_numbers()
		self.mbtypes = qml.representations.get_slatm_mbtypes([self.atomic_numbers])

	def generate(self, structure):
		return qml.representations.generate_slatm(coordinates=structure, nuclear_charges=self.atomic_numbers,mbtypes=self.mbtypes)

class N2P2:
	"""docstring for N2P2 => BPNN model with ACSFs !"""

	def __init__(self, periodic, factor=1, verbose=False):
		self.factor = factor
		self.periodic = periodic
		self.verbose=verbose

	def energy(self, structure, id_rep, p):
		Prt.print_n2p2data(structure,id_rep,self.periodic,verbose=self.verbose)
		file_name = 'input.data.'+str(id_rep)
		p.readStructureFromFile(file_name)
		p.setupSymmetryFunctionStatistics(True, True, False, False)
		p.predict()
		s = p.structure
		e = s.energy
		return e * prm.hartree2kcal * self.factor

	def energy_forces(self, structure, id_rep, p):
		Prt.print_n2p2data(structure,id_rep,self.periodic,verbose=self.verbose)
		file_name = 'input.data.'+str(id_rep)
		p.readStructureFromFile(file_name)
		p.setupSymmetryFunctionStatistics(True, True, False, False)
		p.predict()
		s = p.structure 
		e = s.energy
		f = []
		for atom in s.atoms:
			f.append(atom.f.r)
		f = np.array(f)
		return e * prm.hartree2kcal * self.factor, f * prm.hartree2joule * (1.0/(prm.bohr2angst * prm.angst2m)) * self.factor

class DeepMD:
	"""docstring for DeepMD => BPNN model !"""

	def __init__(self, periodic, atype, factor=1):
		self.factor = factor
		self.periodic = periodic
		self.cell = None
		self.atype = atype

	def deepmd_call(self, mol, dp):
		coordinates = mol.get_positions()
		num_at = len(coordinates)
		coord = np.zeros((num_at,3),dtype=float)
		for j in range(len(coordinates)):
			for k in range(3):
				coord[j][k] = coordinates[j][k]
		coordbis = coord.reshape([1, -1])
		e, f, v = dp.eval(coordbis, self.cell, self.atype)
		return e

	def deepmd_call_f(self, mol, dp):
		coordinates = mol.get_positions()
		num_at = len(coordinates)
		coord = np.zeros((num_at,3),dtype=float)
		for j in range(len(coordinates)):
			for k in range(3):
				coord[j][k] = coordinates[j][k]
		coordbis = coord.reshape([1, -1])
		e, f, v = dp.eval(coordbis, self.cell, self.atype)
		return f[0]

	def deepmd_call_ef(self, mol, dp):
		coordinates = mol.get_positions()
		num_at = len(coordinates)
		coord = np.zeros((num_at,3),dtype=float)
		for j in range(len(coordinates)):
			for k in range(3):
				coord[j][k] = coordinates[j][k]
		coordbis = coord.reshape([1, -1])
		e, f, v = dp.eval(coordbis, self.cell, self.atype)
		return e, f[0]
	
	def energy(self, state, dp):
		energy = self.deepmd_call(state, dp)
		return energy * prm.ev2kcalmol * self.factor

	def force(self, state, dp):
		forces = self.deepmd_call_f(state, dp)
		return forces * prm.ev2kcalmol * (1.0/prm.angst2m) * self.factor

	def energy_forces(self, state, dp):
		energy, forces = self.deepmd_call_ef(state, dp)
		return energy * prm.ev2kcalmol * self.factor, forces * prm.ev2joule * (1.0/prm.angst2m) * self.factor

class MixedPotential:
	"""docstring for MixedPotential_with_forces"""

	def __init__(self, energy_functions_list, forces=False, weights=None, verbose=False):
		self.ener_list = energy_functions_list
		self.verbose = verbose
		if forces:
			self.for_list = self.ener_list
		if weights is None:
			self.weights = np.ones(len(energy_functions_list))
		else:
			self.weights = weights

	def energy_reservoir(self, mol, dp, weights=0.0, verbose=False):
		ener_baselined = self.ener_list[0](mol)
		# In case where forces are computed, we extract only energies
		try:
			ener_baselined = ener_baselined[0]
		except:
			pass
		if verbose:
			print('VERBOSE: Reservoir energy computation: is only baselined computation')
			print('VERBOSE: Structure energy: {}'.format(ener_baselined))
		return ener_baselined

	def energy_lkr(self, mol, weights=1.0, verbose=False):
		energies = []
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		energies.append(self.ener_list[0](mol))
		t2 = time.time()
		energies.append(self.ener_list[1](mol))
		t3 = time.time()
		if verbose:
			print('VERBOSE: Total time for energy = {} sec '.format(t2-t1))
		e = energies[0] + (energies[1]*weights)
		return e, t2-t1, t3-t2

	def energy_n2p2(self, mol, weights=1.0, id_rep=0, p=1, verbose=False):
		energies = []
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		energies.append(self.ener_list[0](mol))
		t2 = time.time()
		energies.append(self.ener_list[1](mol, id_rep, p))
		t3 = time.time()
		if verbose:
			print('VERBOSE: Total time for energy = {} sec '.format(t2-t1))
		e = energies[0] + (energies[1]*weights)
		return e, t2-t1, t3-t2

	def energy_hres_n2p2(self, mol, weights=1.0, id_rep=0, p=1, verbose=False):
		energies = []
		# hRES case: pynnp has to be reload as it is not pickabled
		if p == 1:
			import pynnp
			p = pynnp.Prediction()
			p.log.writeToStdout = False 
			p.setup()
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		energies.append(self.ener_list[0](mol))
		t2 = time.time()
		energies.append(self.ener_list[1](mol, id_rep, p)) 
		t3 = time.time()
		if verbose:
			print('VERBOSE: Total time for energy = {} sec '.format(t2-t1))
		e = energies[0] + (energies[1]*weights)
		return e, t2-t1, t3-t2

	def energy_dp(self, mol, dp, weights=1.0, verbose=False):
		energies = []
		times = []
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		energies.append(self.ener_list[0](mol))
		t2 = time.time()
		energies.append(self.ener_list[1](mol,dp))
		t3 = time.time()
		if verbose:
			print('VERBOSE: Total time for energy = {} sec '.format(t2-t1))
		e = sum(np.array(energies) * self.weights)
		energy = e[0][0]
		return energy, t2-t1, t3-t2

	def energy_hres_dp(self, mol, dp, weights, verbose=False):
		energies = []
		# hRES case: if dp == 1 => it comes from the REXC _smap function !
		# Need to reload here manually dp as it is not pickabled ...
		if dp == 1:
			from deepmd.infer import DeepPot
			dp = DeepPot('graph.pb')
		if weights == 1.0:
			if self.verbose:
				print('VERBOSE: dp has to reloaded as it is not pickable toward different processors !')
				print('VERBOSE: dp is now loaded and is {}'.format(dp))
		# Baselined is still first, ML correction then
		if self.verbose and weights == 1.0:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		energies.append(self.ener_list[0](mol))
		t2 = time.time()
		energies.append(self.ener_list[1](mol,dp))
		t3 = time.time()
		if self.verbose and weights == 1.0:
			print('VERBOSE: Total time for energy = {} sec '.format(t2-t1))
		scaled_energies = np.array([energies[0],energies[1]*weights], dtype=float)
		if self.verbose:
			print('For weight {} here is the old energies : {}'.format(weights, np.array(energies)))
			print('For weight {} here is the new energies : {}'.format(weights, scaled_energies))	
		return sum(scaled_energies), t2-t1, t3-t2

	def energy_forces_hres_dp(self, mol, dp, weights, verbose=False):
		energies = []
		# hRES case: if dp == 1 => it comes from the REXC _smap function !
		# Need to reload here manually dp as it is not pickabled ...
		if dp == 1:
			from deepmd.infer import DeepPot
			dp = DeepPot('graph.pb')
		if weights == 1.0:
			if self.verbose:
				print('VERBOSE: dp has to reloaded as it is not pickable toward different processors !')
				print('VERBOSE: dp is now loaded and is {}'.format(dp))
		# Baselined is still first, ML correction then
		if self.verbose and weights == 1.0:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
		t1 = time.time()
		e1, f1 = self.ener_list[0](mol)
		t2 = time.time()
		e2, f2 = self.ener_list[1](mol,dp)
		t3 = time.time()
		scaled_energies = np.array([e1,e2*weights], dtype=float)
		np.multiply(f2, weights)
		if self.verbose:
			print('For weight {} here is the old energies : {}'.format(weights, np.array(energies)))
			print('For weight {} here is the new energies : {}'.format(weights, scaled_energies))
		return sum(scaled_energies), np.add(f1, f2), t2-t1, t3-t2

	def energy_forces_dp(self, mol, dp, weights=1.0, verbose=False):
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
			print('VERBOSE: Forces checking => {} and {}'.format(self.for_list[0],self.for_list[1]))
		t1 = time.time()
		e1, f1 = self.ener_list[0](mol)
		t2 = time.time()
		e2, f2 = self.ener_list[1](mol,dp)
		t3 = time.time()
		if verbose:
			print('VERBOSE: Baselined computed in {} sec'.format(t2-t1))
			print('VERBOSE: ML correction computed in {} sec'.format(t3-t2))
			print('VERBOSE: Total time for energy and forces = {} sec'.format(t3-t1))
		np.multiply(f2, self.weights[1])
		e = e1+(e2*weights)
		energy = e[0][0]
		return energy, np.add(f1, f2), t2-t1, t3-t2

	def energy_forces_n2p2(self, mol, weights=1.0, id_rep=0, p=1, verbose=False):
		# Baselined is still first, ML correction then
		if verbose:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
			print('VERBOSE: Forces checking => {} and {}'.format(self.for_list[0],self.for_list[1])) 
		t1 = time.time()
		e1, f1 = self.ener_list[0](mol)
		t2 = time.time()
		e2, f2 = self.ener_list[1](mol,id_rep, p)
		t3 = time.time()
		if verbose:
			print('VERBOSE: Baselined computed in {} sec'.format(t2-t1))
			print('VERBOSE: ML correction computed in {} sec'.format(t3-t2))
			print('VERBOSE: Total time for energy and forces = {} sec'.format(t3-t1)) 
		np.multiply(f2, weights)
		e = e1+(e2*weights)
		return e, np.add(f1, f2), t2-t1, t3-t2

	def energy_forces_hres_n2p2(self, mol, weights=1.0, id_rep=0, p=1, verbose=False):
		# hRES case: we reload pynnp as it is not pickabled
		if p == 1:
			import pynnp
			p = pynnp.Prediction()
			p.log.writeToStdout = False
			p.setup()
		# Baselined is still first, ML correction then
		if verbose and weights == 1:
			print('VERBOSE: Should be baselined and then ML correction')
			print('VERBOSE: Energy checking => {} and {}'.format(self.ener_list[0],self.ener_list[1]))
			print('VERBOSE: Forces checking => {} and {}'.format(self.for_list[0],self.for_list[1]))
		t1 = time.time()
		e1, f1 = self.ener_list[0](mol)
		t2 = time.time()
		e2, f2 = self.ener_list[1](mol,id_rep, p)
		t3 = time.time()
		if verbose and weights == 1:
			print('VERBOSE: Baselined computed in {} sec'.format(t2-t1)) 
			print('VERBOSE: ML correction computed in {} sec'.format(t3-t2))
			print('VERBOSE: Total time for energy and forces = {} sec'.format(t3-t1))
		np.multiply(f2, weights)
		e = e1+(e2*weights)
		return e, np.add(f1, f2), t2-t1, t3-t2

	def energy_forces_mts_n2p2(self, mol, weights=1.0, id_rep=0, p=1, direct=False, baselined=True, verbose=False):
		# MTS Case: decide if it is direct or baselined !
		if direct:
			t1 = time.time()
			t2 = time.time()
			os.chdir('DIRECT/')
			e2, f2 = self.ener_list[1](mol,id_rep, p)
			os.chdir('../')
			t3 = time.time()
			e = e2
			f = np.multiply(f2, weights)
		else:
			t1 = time.time()
			e1, f1 = self.ener_list[0](mol)
			t2 = time.time()
			os.chdir('BASELINED/')
			e2, f2 = self.ener_list[1](mol,id_rep, p)
			os.chdir('../')
			t3 = time.time()
			e = e1 + e2
			np.multiply(f2, weights)
			f = np.add(f1, f2)
		if verbose and weights == 1:
			print('VERBOSE: Baselined computed in {} sec'.format(t2-t1))
			print('VERBOSE: ML correction computed in {} sec'.format(t3-t2))
			print('VERBOSE: Total time for energy and forces = {} sec'.format(t3-t1))
		return e, f, t2-t1, t3-t2

