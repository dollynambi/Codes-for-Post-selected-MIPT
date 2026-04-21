import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import qr
from scipy.linalg import norm
from numpy import random
import cmath
from itertools import groupby, count
import time
from numpy import loadtxt
import sys
import os
import math
import argparse


def create_bell_tensor_sqrt(D):
    d_frac = ((D-1)/2)**1
    bell_final = np.zeros((3, 3))
    bell_final[0, 0] = (1-d_frac) + d_frac/np.sqrt(3)
    bell_final[1, 1] = d_frac/np.sqrt(3)
    bell_final[2, 2] = d_frac/np.sqrt(3)
    return np.sqrt(bell_final/np.linalg.norm(bell_final))


def create_random_tensor(n_legs):
    return np.random.randn(*[3 for _ in range(n_legs)])


def gate_action(control, target, gate, reg):
    """Apply a two-site gate to the control and target indices of ``reg``."""
    result = np.tensordot(gate, reg, ((2, 3), (control, target)))
    result = np.moveaxis(result, (0, 1), (control, target))
    return result


def create_tensor_grid(N, t):
    """Return an ``N x t`` grid filled with one fixed rank-4 qutrit tensor."""
    tensor_grid = np.empty((N, t), dtype=object)
    tensor_fixed = create_random_tensor(4)
    for j in range(t):
        for i in range(N):
            tensor_grid[i, j] = tensor_fixed
    return tensor_grid


def bell_contract_tensors_sqrt(N, t, tensors, D):
    """Dress each leg of every local tensor with the square-root Bell factor."""
    bell = create_bell_tensor_sqrt(D)
    for j in range(t):
        for i in range(N):
            tensors[i, j] = np.tensordot(tensors[i, j], bell, ([0], [0]))
            tensors[i, j] = np.moveaxis(tensors[i, j], 3, 0)
            tensors[i, j] = np.tensordot(tensors[i, j], bell, ([1], [0]))
            tensors[i, j] = np.moveaxis(tensors[i, j], 3, 1)
            tensors[i, j] = np.tensordot(tensors[i, j], bell, ([2], [0]))
            tensors[i, j] = np.moveaxis(tensors[i, j], 3, 2)
            tensors[i, j] = np.tensordot(tensors[i, j], bell, ([3], [0]))
    return tensors


def ee_bipartite(wave, la, l):
    """Compute the bipartite entropy across a cut of size ``la``."""
    lb = l-la
    temp = np.reshape(wave, (3**la, 3**lb))
    sp = np.linalg.svd(temp, compute_uv=False)
    sp = sp[np.nonzero(sp)]
    el = sp**2
    von = -np.dot(el, (np.log2(el)))
    return von


def ee_general(psi, sites, N):
    psi = np.moveaxis(psi, sites, np.arange(len(sites)))
    ee = ee_bipartite(psi, len(sites), N)
    return ee


def tmi(psi, qlen):
    """Compute TMI for four equal contiguous blocks."""
    A = np.arange(qlen//4)
    S_A = ee_general(psi, A, qlen)
    B = np.arange(qlen//4, 2*qlen//4)
    S_B = ee_general(psi, B, qlen)
    C = np.arange(2*qlen//4, 3*qlen//4)
    S_C = ee_general(psi, C, qlen)
    D = np.arange(3*qlen//4, qlen)
    AB = np.concatenate([A, B])
    S_AB = ee_general(psi, AB, qlen)
    BC = np.concatenate([B, C])
    S_BC = ee_general(psi, BC, qlen)
    AC = np.concatenate([A, C])
    S_AC = ee_general(psi, AC, qlen)
    ABC = np.concatenate([A, B, C])
    S_ABC = ee_general(psi, ABC, qlen)
    tripmi = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    return tripmi


def tmi_routine(N, t, tensors, D):
    """Evolve to the final depth and average the last-step TMI."""
    obj = np.zeros((3,) * N)
    obj[(0,)*N] = 1
    nq1 = range(N)
    oq1 = nq1[::2]
    eq1 = nq1[1::2]
    for t1 in range(0, t, 2):
        for l1 in oq1:
            obj = gate_action(l1, l1+1, tensors[l1//2, t1], obj)
            obj /= np.linalg.norm(obj)

        if t1 == (t-2):
            tmivst1 = tmi(obj, N)

        for l2 in eq1:
            obj = gate_action((l2-2) %
                              N, l2-1, tensors[(l2-1)//2, t1+1], obj)
            obj /= np.linalg.norm(obj)

        if t1 == (t-2):
            tmivst2 = tmi(obj, N)
            tmivst = (tmivst1+tmivst2)/2

    return tmivst

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run TMI Calculation")
    parser.add_argument('--N', type=int, required=True, help="System size N")
    parser.add_argument('--D', type=float, required=True,
                        help="Bond dimension D")
    parser.add_argument('--output', type=str, required=True,
                        help="Output file path")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    N = args.N
    D = args.D
    output_file = args.output

    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Simulation parameters
    t = 8*N  # circuit depth
    n_realizations = 5000  # number of realizations
    int_D = int(D) if D == math.floor(D) else int(math.floor(D) + 1)
    frac_D = D - math.floor(D)

    tmi_record = np.zeros(n_realizations)

    start_time = time.perf_counter()
    with open(output_file, 'w') as f:
        for iter in range(n_realizations):
            print(iter)
            tensors = create_tensor_grid(N//2, t)
            tensors_partial = bell_contract_tensors_sqrt(N//2, t, tensors, D)
            tmi_record[iter] = tmi_routine(N, t, tensors_partial, D)
            f.write(f"{tmi_record[iter]}\n")
            f.flush()

    end_time = time.perf_counter()
    print("total time", end_time - start_time)
