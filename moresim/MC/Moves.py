"""
Module with Monte Carlo moves.
"""


class RandomParticleMove:
    """Constructor of random particle MC moves"""

    def __init__(self, rand_num_gen, part_indexs='all'):
        self.rand_num_gen = rand_num_gen
        self.part_indexs = part_indexs

    def move(self, old_position):
        if self.part_indexs == 'all':
            new_pos = old_position + self.rand_num_gen.generate(
                old_position.shape)

        return new_pos
