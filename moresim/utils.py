"""
This module includes useful functions for the analysis of replica exchange and other sampling simulations.
"""
import ase
import numpy as np
import math
from numpy.linalg import norm
from functools import reduce
import os
import sp

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


def factors(n):
    """Return factors of number n."""
    return set(
        reduce(list.__add__,
               ([i, n // i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))


def normal_to_svd_plane(points):
    """Return the normal vector to the rank-2 SVD plane."""
    points = np.asarray(points)
    centroid = np.mean(points, axis=0)
    rpoints = np.subtract(points, centroid)
    U, s, V = np.linalg.svd(rpoints)
    normal_vector = V[len(V) - 1, :]
    return normal_vector


def svd_vectors(points):
    """Return the normal vector to the rank-2 SVD plane."""
    points = np.asarray(points)
    centroid = np.mean(points, axis=0)
    rpoints = np.subtract(points, centroid)
    U, s, V = np.linalg.svd(rpoints)
    normal_vector = V[- 1, :]
    v1 = V[0, :]
    v2 = V[1, :]

    return v1, v2, normal_vector


def point_projection_to_plane(point, v1, v2, planep):
    """Projection of point to plane."""
    proj1 = np.dot(point, v1)
    proj2 = np.dot(point, v1)

    proj_point = v1 * proj1 + v2 * proj2 + planep

    return proj_point


def make_right_hand(v1, v2, v3):
    """Return a right hand basis."""
    xx = np.cross(v1, v2)
    return np.sign(np.dot(xx, v3)) * v3


def plane_angle(normal_vecor1, normal_vecor2):
    """Docstring for plane_angle."""
    angle = math.degrees(np.arccos(np.dot(normal_vecor1, normal_vecor2)) /
                         (norm(normal_vecor1) * norm(normal_vecor2)))
    return angle


def unit_vector(vector):
    """Return the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Return the angle in radians between vectors 'v1' and 'v2'."""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def compute_weights(exp_energies, cheap_energies,
                    temperature, method='metropolis', units='kcal'):
    """Weights for free energy perturbation."""
    cheap_energies = np.asarray(cheap_energies)
    exp_energies = np.asarray(exp_energies)

    if units == 'atomic':
        kb = 1
    elif units == 'kcal':
        kb = 0.0019872041

    beta = 1 / (kb * temperature)
    if method == 'metropolis':
        weights = np.exp(- beta * (exp_energies - cheap_energies))

    elif method == 'bennet':
        weights = 1 / (1 + np.exp(- beta * (cheap_energies - exp_energies)))

    else:
        Exception('Select valid method')

    return weights / sum(weights)


def rel_free_energy(labels, label, label0=None, weights=None,
                    temp=300, kb=0.0019872041):
    """Relative free energy of cluster with specific label."""
    if weights is None:
        weights = np.ones(len(labels))

    counts = dict()

    for i in np.unique(labels):
        counts[i] = sum(weights[labels == i])

    fes = {k: - kb * temp * np.log(v) for k, v in counts.items()}

    rel_fe = fes[label]
    if label0 is not None:
        rel_fe -= fes[label0]
    return rel_fe


def bootstraping(coord1, coord2, bins, num_iter=100, weights=None, len_bt=1):
    """Compute histograms with bootstraping"""
    if len(coord1) != len(coord2):
        raise ValueError('Length of coord1 must be equal to length of coord2')

    num_items = len(coord1)
    len_bt = int(len_bt * num_items)
    # bt_indices_list = []

    if weights is None:
        weights = np.ones(num_items)

    histogram_list = []

    for i in range(num_iter):
        bt_indices = np.random.choice(num_items, len_bt)
        # bt_indices_list.append(bt_indices)

        bt_coord1 = [coord1[i] for i in bt_indices]
        bt_coord2 = [coord2[i] for i in bt_indices]
        bt_weigths = [weights[i] for i in bt_indices]

        min1 = min(coord1)
        max1 = max(coord1)

        min2 = min(coord2)
        max2 = max(coord2)

        histogram, x, y = np.histogram2d(bt_coord1, bt_coord2, bins=bins,
                                         range=[[min1, max1], [min2, max2]],
                                         weights=bt_weigths)

        histogram_list.append(histogram)

    return np.asarray(histogram_list)


def bootstraping_indices(num_items, len_bt=1, num_iter=100):
    """Obtain indices of partitions from bootstraping."""
    len_bt = int(len_bt * num_items)

    bt_indices_list = []
    for i in range(num_iter):
        bt_indices = np.random.choice(num_items, len_bt)
        bt_indices_list.append(bt_indices)

    return bt_indices_list


def analyze_bt_results(indices_list, coords1, bins, coords2, weights=None):
    """Return 2d histograms of the number of counts using the bt results."""
    num_iter = indices_list.shape[0]
    len_bt = indices_list.shape[1]

    if len(coords1) != len(coords2):
        raise ValueError(
            'Length of coords1 must be equal to length of coords2')
    if weights is None:
        weights = np.ones(len_bt)

    histogram_list = []
    for i in range(num_iter):
        bt_indices = indices_list[i, :]

        bt_coord1 = [coords1[i] for i in bt_indices]
        bt_coord2 = [coords2[i] for i in bt_indices]
        bt_weigths = [weights[i] for i in bt_indices]

        min1 = min(coords1)
        max1 = max(coords1)

        min2 = min(coords2)
        max2 = max(coords2)

        histogram, x, y = np.histogram2d(bt_coord1, bt_coord2, bins=bins,
                                         range=[[min1, max1], [min2, max2]],
                                         weights=bt_weigths)

        histogram_list.append(histogram)

    return np.asarray(histogram_list)


def counts_bt(labels, weights=None, iters=100, len_bt=1):
    """Compute diccionary of labels and weighted counts of labels."""
    if weights is None:
        weights = np.ones(len(labels))

    num_items = len(labels)
    len_bt = int(len_bt * num_items)

    clusters = np.unique(labels)
    num_clusters = len(clusters)

    cluster_counts_dict_list = []

    for bt_iter in range(iters):

        idx = np.random.choice(num_items, len_bt)
        bt_weights = weights[idx]
        bt_labels = labels[idx]

        cluster_counts_dict = {}
        for i in range(num_clusters):
            counts = sum(bt_weights[bt_labels == clusters[i]])
            cluster_counts_dict[clusters[i]] = counts
        cluster_counts_dict_list.append(cluster_counts_dict)

    return cluster_counts_dict_list


def remove_file_if_exists(file):
    try:
        os.remove(file)
    except OSError:
        pass


def blocking_method(signal):
    """Blocking method to compute convergence and error estimate of free energies."""
    block_stds = []
    block_stds_stds = []

    while(len(signal) > 5):

        n = len(signal)

        c0 = np.var(signal)

        block_std = np.sqrt(c0 / (n - 1))
        block_std_std = block_std / np.sqrt(2 * (n - 1))

        block_stds.append(block_std)
        block_stds_stds.append(block_std_std)
        signal = block_renormalization(signal)

    return block_stds, block_stds_stds


def block_renormalization(signal):
    """Renormalization step of blocking method."""
    signal = [0.5 * (signal[2 * i] + signal[2 * i + 1])
              for i in range(len(signal) // 2)]

    return signal


def block_averages(signal, num_dec, min_num_blocks=5):
    """Compute bloack averages for blocking method"""
    num_samples = len(signal)
    ptbs = []
    block_sizes = []
    for num_blocks in range(min_num_blocks, num_dec, 1):

        block_size = num_samples // num_blocks

        blocks = [signal[i:i + block_size] for i in range(
            0, len(signal), block_size)]

        block_means = [np.mean(block) for block in blocks]

        var_blocks = np.var(block_means)
        ptbs.append(block_size * var_blocks)
        block_sizes.append(block_size)
    ptbs = np.asarray(ptbs) / np.var(signal)

    return ptbs, np.asarray(block_sizes)


def p_convergence(labels, labelx, init=0, step=10):
    """Probability convergence of each labelled cluster."""
    num_items = len(labels)
    converging_p = []

    for i in range(init, num_items, step):
        temp = labels[0:i]
        numx = temp.count(labelx)
        converging_p.append(numx / i)

    return converging_p


def p_convergence_weights(labels, labelx, weights, init=100, step=10):
    """Probability convergence of each labelled cluster with weights."""
    num_items = len(labels)
    converging_p = []

    for i in range(init, num_items, step):
        temp_lab = labels[0:i]
        temp_wei = weights[0:i]

        count = np.sum(temp_wei[temp_lab == labelx])

        converging_p.append(count / np.sum(temp_wei))

    return converging_p


def moving_p(labels, labelx, window, init=0, step=10):
    """Moving average of probabily for each labelled cluster."""
    num_items = len(labels)
    moving_p = []

    for i in range(init, num_items - window, step):
        temp = labels[i:(i + window)]
        numx = temp.count(labelx)
        moving_p.append(numx / window)

    return moving_p


def mol_cm(mol):
    """Return center of mass of ase Atoms system"""
    cm = np.mean(mol.positions, axis=0)
    return cm


def mol_cm_distance(mol1, mol2):
    """
    Return the distance betwen the center of masses of two ase Atoms systems
    """

    cm1 = np.mean(mol1.positions, axis=0)
    cm2 = np.mean(mol2.positions, axis=0)

    return np.linalg.norm(cm1 - cm2)


def mol_min_distance(mol1, mol2):
    """Compute the minimum distance between two atoms in two ase Atoms systems"""
    dist_mat = sp.spatial.distance_matrix(mol1.positions, mol2.positions)
    return np.min(dist_mat)
