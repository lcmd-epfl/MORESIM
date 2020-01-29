"""
This module includes different useful functions.
"""

import numpy as np
import ase


def mc_prob(weight_new, weight_old):

    prob = min([1, np.exp(- weight_new + weight_old)])

    return prob


def maxwell_boltzmann_distribution(atoms, temperature):
    input_velocities = np.zeros(shape=[len(atoms), 3])
    T = temperature
    kb = 1.38064852 * 1e-23
    Na = 6.02214086 * 1e23

    mass = [1e-3 * M / Na for M in ase.Atoms(atoms).get_masses()]
    standard_deviation = [np.sqrt((kb * T) / m) for m in mass]

    for i in range(len(standard_deviation)):
        for j in range(3):
            input_velocities[i][j] = 1e-5 * np.random.normal(
                loc=0, scale=standard_deviation[i], size=[1])
    return input_velocities
