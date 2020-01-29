"""
This module includes different kind of simulations.
"""

import ase.io as aio
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time
from moresim import utils
import ase


class Trajectory:
    """ A class to store trajectoreis
    """

    def __init__(self, state_traj=None, energy_traj=None,
                 generation_details=None,
                 flush_prefix=None):
        """Generate a trajectory

        Parameters
        ----------

        state_traj:  List of ase Atoms objects

        state_traj: 1D-array with potential energies

        generation_details: Dictionary storing arbitrary details of how the
        trajectory was simulated e.g. temperature

        """
        if state_traj is None:
            state_traj = []
        self.state_traj = state_traj

        if energy_traj is None:
            energy_traj = [state.energy for state in state_traj]
        self.energy_traj = energy_traj

        if generation_details is None:
            generation_details = {}
        self.generation_details = generation_details

        # if moves_used is None:
        #     moves_used = []
        # self.moves_used = moves_used

        # if moves_accepted is None:
        #     moves_accepted = []
        # self.moves_accepted = moves_accepted

    def extend(self, traj):
        if type(traj) is not type(self):
            raise ValueError('The input is not a trajectory')

        # if traj.generation_details != self.generation_details:
        #     raise ValueError(
        #         'The trajectories to merge come from different simulations.')

        self.state_traj.extend(traj.state_traj)
        self.energy_traj.extend(traj.energy_traj)
        # self.moves_used.extend(traj.moves_used)
        # self.moves_accepted.extend(traj.moves_accepted)

    # def mc_probabilities(self):
    #     probabilities = []
    #     for i in range(len(self.generation_details['move_list'])):
    #         idxs = [t for t, x in enumerate(self.moves_used) if x == i]
    #         idxs_bool = [self.moves_accepted[t] for t in idxs]

    #         probabilities.append(sum(idxs_bool) / len(idxs_bool))
    #     return probabilities

    def flush(self, flush_prefix):
        # if len(self.moves_used) > 0:
            # f = open('{}_mc_moves.dat'.format(flush_prefix), 'ab')
            # np.savetxt(f, np.array(
            #     list(zip(self.moves_used, self.moves_accepted))), fmt='%i')
            # f.close()

        f = open('{}_energies.dat'.format(flush_prefix), 'ab')
        np.savetxt(f, np.array(self.energy_traj), fmt='%.6f')
        f.close()

        for struct in self.state_traj:
            aio.write('{}_structures.xyz'.format(flush_prefix),
                      ase.Atoms(struct.get_chemical_symbols(),
                                positions=struct.positions), append=True)

        self.__init__(generation_details=self.generation_details,
                      flush_prefix=flush_prefix)


class MCSimulation:
    """ A Monte Carlo simulation constructor
    """

    def __init__(self, energy_func, temperature,
                 move_list, move_weight_list=None, kb=0.0019872041):
        """Define parameters of MC simulation

        Parameters
        ----------

        energy_func: Hamiltonian of the simulations.
        It should act as: energy_func(ase_state) -> float (energy)

        temperature: Temperature of the simulation

        move_list: List with MCmove objects. An MCmove should act as
        MCmove(atom_positions) -> new_atom_positions

        move_weight_list: probability weight for the different MCmoves

        kb: The boltzman constant. Defined to work in kcal/mols
        """
        self.temperature = temperature
        self.beta = (kb * self.temperature) ** - 1
        # self.atoms = atoms
        self.energy_func = energy_func
        self.move_list = move_list
        self.move_weight_list = move_weight_list

    def simulation_details(self):
        details = {'temperature': self.temperature}
        return details

    def _advance(self, old_state):

        move_idx = np.random.choice(
            list(range(len(self.move_list))), p=self.move_weight_list)
        move = self.move_list[move_idx]

        old_pos = old_state.positions
        old_ener = old_state.energy

        t1 = time.time()

        new_pos = move.move(old_pos)

        new_ener = self.energy_func(ase.Atoms(self.atoms, new_pos))
        t2 = time.time()

        new_weight = self.beta * new_ener
        old_weight = self.beta * old_ener

        prob = utils.mc_prob(weight_new=new_weight, weight_old=old_weight)
        accepted = np.random.rand() < prob

        # print((old_ener, new_ener))
        # print((prob, accepted))
        if not accepted:
            new_pos = old_pos
            new_ener = old_ener

        new_state = ase.Atoms(self.atoms, positions=new_pos)
        new_state.energy = new_ener

        return new_state, move_idx, accepted, t2 - t1

    def run(self, init_state, steps, stride=10,
            return_last=False, *args):
        self.atoms = init_state.get_chemical_symbols()
        np.random.seed()

        if 'energy' in dir(init_state):
            ener = init_state.energy
        else:
            ener = self.energy_func(init_state)
            init_state.energy = ener
        traj = []
        moves_used = []
        moves_accepted = []
        times = []

        state = init_state
        for i in range(steps):
            state, move_idx, accepted, tt = self._advance(state)

            moves_used.append(move_idx)
            moves_accepted.append(accepted)
            times.append(tt)

            if i % stride == 0:
                traj.append(state)

        traj = Trajectory(traj, generation_details=self.simulation_details())

        # print(time.time() - t1)
        if return_last is True:
            return [traj, state, np.mean(times)]

        else:
            return traj


class Reservoir:
    """ A constructor for Reservoirs
    """

    def __init__(self, structures_folder, energies,
                 temperature, energy_func, kb=0.0019872041):
        """Define parameters of a Reservoir

        Parameters
        ----------

        structures_folder: Directory with xyz structures of the
        reservoir stored. Each structures should be stored in a file
        starting by coord_ and followed by 5 leading zeros.

        energies: energies of each structure in the reservoir

        temperature: Temperature of the simulation

        energy_func: Hamiltonian of the simulations.
        It should act as: energy_func(ase_state) -> float (energy)

        kb: The boltzman constant. Defined to work in kcal/mols
        """
        self.structures_folder = structures_folder
        self.energies = energies
        self.size = len(energies)
        self.temperature = temperature
        self.beta = (kb * self.temperature) ** - 1
        self.energy_func = energy_func

    def simulation_details(self):
        details = {'temperature': self.temperature}
        return details

    def flush(self):
        pass

    def run(self, *args):
        np.random.seed()
        idx = np.random.choice(np.arange(self.size))

        mol = aio.read(self.structures_folder + '/coord_{:05d}'.format(idx))
        traj = Trajectory(state_traj=[mol],
                          generation_details=self.simulation_details())

        ener = self.energies[idx]
        mol.energy = ener
        time = 0
        return [traj, mol, time]


class ReplicaExchangeSimulation:
    """ A general constructor for Replica Exchange Simulation
    """

    def __init__(self, num_reps, simulations, init_states, stride, rep_steps,
                 num_processors=1):
        """Define parameters of Replica Exchange Simulation

        Parameters
        ----------

        num_reps: number of replicas in the simulation.
        It must be a pair number.

        simulations: List with simulations objects e.g. MCSimulation
        Each one should have a different temperature for T-RE, or different
        energy functions for H-RE, or either.

        init_states: List of ase.Atoms states.

        rep_steps: Number of steps before attemting to swap states

        num_processors: Number of processors used in the simulation.
        """
        directory = '.'
        self.num_reps = num_reps
        self.num_processors = num_processors
        if num_reps % 2 != 0:
            raise('Number of replicas must be pair')

        if len(simulations) != self.num_reps:
            raise('Wrong number of temperatures')

        self.temperatures = [sim.temperature for sim in simulations]

        self.energy_funcs = [sim.energy_func for sim in simulations]

        self.simulations = simulations

        self.in_rep_states = init_states

        self.par_exec = ProcessPoolExecutor(max_workers=num_processors)

        if sum(['energy' in dir(
                state) for state in init_states]) == self.num_reps:
            pass
        else:
            self.in_rep_eners = list(self.par_exec.map(
                _smap, self.energy_funcs, self.in_rep_states))

        self.rep_index = np.arange(self.num_reps)

        self.even_sims = self.rep_index[::2]

        self.odd_sims = self.rep_index[::2]

        self.accepted_exchanges = {(i, (i + 1) % self.num_reps):
                                   [] for i in range(self.num_reps)}

        self.strides = [stride for i in range(num_reps)]

        self.rep_steps = rep_steps

        for stride in self.strides:
            if self.rep_steps % stride != 0:
                raise ValueError('Rep_steps must be multiple of stride')

        self.rep_stepss = [rep_steps for i in range(self.num_reps)]

        self.directory = directory

    def run(self, num_exchanges):

        trajectories = [Trajectory(generation_details=sim.simulation_details())
                        for sim in self.simulations]

        for i in range(num_exchanges):

            t0 = time.time()

            # geenerate dynamics
            # run individual simulation in parallel
            return_last = [True for l in range(self.num_reps)]
            # print('num in -rep states', len(self.in_rep_states))
            simulation_results = list(
                self.par_exec.map(_run_simulation, self.simulations,
                                  self.in_rep_states, self.rep_stepss,
                                  self.strides, return_last))

            rep_trajs = [res[0] for res in simulation_results]
            exchange_states = [res[1] for res in simulation_results]
            sim_times = [res[2] for res in simulation_results]

            for k in range(self.num_reps):
                trajectories[k].extend(rep_trajs[k])

            self.in_rep_states = self._replica_exchange(
                exchange_states)

            self.exchange_probabilities = {key: (0.001 + sum(val)) / (len(
                val) + 0.001) for key, val in self.accepted_exchanges.items()}

            t1 = time.time()

            with open("exchange.txt", "a") as myfile:
                myfile.write(
                    'Exchange {0}, step {1}, time interval {2:.3} \n'.format(
                        i + 1, (i + 1) * self.rep_steps, t1 - t0))

                [myfile.write('{}s, '.format(
                    x)) for x in sim_times]

                myfile.write('\n')

                [myfile.write('{0}: {1:.3}\n'.format(
                    x, y)) for x, y in self.exchange_probabilities.items()]

            # print(len(trajectories[0].energy_traj))
            if i % 10 == 9:
                # t2 = time.time()
                for rep, traj in enumerate(trajectories):
                    # print(traj)
                    traj.flush(flush_prefix=(
                        self.directory + '/ms.rep{}_'.format(rep)))
                # t3 = time.time()
                # print(t3 - t2)

    def _replica_exchange(self, exchange_states):
        shift = np.random.choice([1, -1])

        rep_index = np.arange(self.num_reps)

        group1 = rep_index[::2]
        group2 = rep_index[1::2]
        exchange_structs = [xx.positions for xx in exchange_states]
        exchange_eners = [xx.energy for xx in exchange_states]

        if shift == 1:
            ex_index = np.vstack((group2, group1)).flatten(order='F')

        else:

            ex_index = np.roll(
                np.vstack((group1, np.roll(group2, 1))).flatten(
                    order='F'), -1)

        pairs = list(zip(group1, ex_index[::2]))

        old_structs = exchange_structs
        old_energies = exchange_eners
        # old_velocities = [xxx.velocities if 'velocities' in dir(
        #     xxx) else None for xxx in exchange_states]

        new_structs = [old_structs[i] for i in ex_index]
        # new_velocities = [old_velocities[i] * np.sqrt(
        #         self.temperatures[j] / self.temperatures[i])
        #                 for j, i in enumerate(
        #     ex_index)]

        # launch energy on struct for each one
        new_states = [ase.Atoms(
            symbols=exchange_states[0].get_chemical_symbols(
            ), positions=new_structs[i]) for i in range(self.num_reps)]

        new_energies = list(self.par_exec.map(
            _smap, self.energy_funcs, new_states))
        # #
        # print(['{0:.5f}'.format(x) for x in old_energies])
        # print(['{0:.5f}'.format(x) for x in new_energies])

        for pair in pairs:

            rep0 = self.simulations[pair[0]]
            rep1 = self.simulations[pair[1]]

            old_e0 = old_energies[pair[0]]
            old_e1 = old_energies[pair[1]]

            new_e0 = new_energies[pair[0]]
            new_e1 = new_energies[pair[1]]

            old_weight = rep0.beta * old_e0 + rep1.beta * old_e1
            new_weight = rep0.beta * new_e0 + rep1.beta * new_e1

            prob = utils.mc_prob(weight_new=new_weight, weight_old=old_weight)
            accepted = np.random.rand() < prob

            # print(pair[0], pair[1], '{0:.2f}'.format(prob), accepted)

            if shift == 1:
                self.accepted_exchanges[(pair[0], pair[1])].append(accepted)
            else:
                self.accepted_exchanges[(pair[1], pair[0])].append(accepted)

            if accepted:
                pass
            else:
                new_structs[pair[0]] = old_structs[pair[0]]
                new_structs[pair[1]] = old_structs[pair[1]]

                new_energies[pair[0]] = old_energies[pair[0]]
                new_energies[pair[1]] = old_energies[pair[1]]

                # new_velocities[pair[0]] = old_velocities[pair[0]]
                # new_velocities[pair[1]] = old_velocities[pair[1]]

        new_states = [ase.Atoms(
            symbols=exchange_states[0].get_chemical_symbols(
            )) for i in range(self.num_reps)]

        for i, state in enumerate(new_states):
            state.positions = new_structs[i]
            state.energy = new_energies[i]
            # if new_velocities[i] is None or old_velocities[i] is None:
            #     pass
            # else:
            #     state.velocities = new_velocities[i]

        return new_states


def _smap(f, *args):
    return f(*args)


def _run_simulation(simulation, *args):
    return simulation.run(*args)
