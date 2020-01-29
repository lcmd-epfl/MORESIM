"""
This module includes different energy functions.
"""
import ase
import ase.calculators.dftb as adftb
import scipy as sp
import scipy.spatial
import numpy as np


class Energy_function:
    """Template for energy_functions to be used within this module"""

    def __init__(self, arg):
        self.arg = arg

        def energy(self, state):
            return 0

        def force(self, state):
            return 0


class Lennard_Jones(Energy_function):
    """Basic Lennard_Jones potential"""

    def __init__(self, sigma, epsilon=1):
        self.sigma = sigma
        self.epsilon = epsilon

    def energy(self, state):
        dist_mat = sp.spatial.distance_matrix(state.positions,
                                              state.positions)

        vals = dist_mat[np.triu_indices_from(dist_mat, 1)]

        repulsive = 4 * self.epsilon * (self.sigma / vals) ** 12
        attractive = - 4 * self.epsilon * (self.sigma / vals) ** 6

        return np.sum(repulsive + attractive)


class DftbEnergy(Energy_function):
    """docstring for dftb"""

    def __init__(self, atoms, directory, **kwargs):
        self.dftb_kwargs = kwargs
        self.atoms = ase.Atoms(atoms)
        self.directory = directory
        self.calc = adftb.Dftb(**kwargs)
        self.calc.directory = directory
        if 'factor' in kwargs:
            self.factor = kwargs['factor']
        else:
            self.factor = 1
        self.ev2kcal = 23.21

    def energy(self, state):

        self.calc.calculate(state)
        energy = self.calc.results['energy']
        # energy = 0
        return energy * self.ev2kcal * self.factor

    def force(self, state):

        self.calc.calculate(state)
        energy = self.calc.results['energy']
        force = self.calc.results['force']
        # energy = 0
        return energy * self.ev2kcal * self.factor, force * self.factor
