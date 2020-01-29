"""
This module contains random number generators with different distributions.
"""

import numpy as np


class RandGen:
    """Required form of the generator for this library."""

    def __init__(self):
        pass

    def generate(self, size):
        return None


class GaussianVar(RandGen):
    """Normally distributed random numbers generator."""

    def __init__(self, loc, var):
        self.loc = loc
        self.var = var

    def generate(self, size):
        return np.random.normal(self.loc, self.var, size)
