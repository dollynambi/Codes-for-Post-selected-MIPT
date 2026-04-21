import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import qr
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
    """Return an ``N x t`` grid of random rank-4 qutrit tensors."""
    tensor_grid = np.empty((N, t), dtype=object)
    for j in range(t):
        for i in range(N):
            tensor_grid[i, j] = create_random_tensor(4)
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


def bip_mi(psi, N, B, C):
    """Compute the mutual information between intervals ``B`` and ``C``."""
    S_B = ee_general(psi, B, N)
    S_C = ee_general(psi, C, N)
    BC = np.concatenate([B, C])
    S_BC = ee_general(psi, BC, N)
    bip_mi = S_B + S_C - S_BC
    return bip_mi

def mi_routine(N, t, tensors):
    obj = np.zeros((3,) * N)
    obj[(0,)*N] = 1
    nq1 = range(N)
    oq1 = nq1[::2]
    eq1 = nq1[1::2]
    for t1 in range(0, t, 2):
        for l1 in oq1:
            obj = gate_action(l1, l1+1, tensors[l1//2, t1], obj)
            obj /= np.linalg.norm(obj)

        for l2 in eq1:
            obj = gate_action((l2-2) %
                              N, l2-1, tensors[(l2-1)//2, t1+1], obj)
            obj /= np.linalg.norm(obj)

    return obj


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run mutual information Calculation")
    parser.add_argument('--N', type=int, required=True, help="System size N")
    parser.add_argument('--D', type=float, required=True,
                        help="Bond dimension D")
    parser.add_argument('--output', type=str, required=True,
                        help="Output folder")
    parser.add_argument('--a', type=int, required=True, help="ath eta")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    N = args.N
    D = args.D
    output_dir = args.output
    a = args.a

    selected_data = []
    x1 = 0
    x3 = N//2

    for x2 in np.arange(x1+2, x3, 1):
        for x4 in np.arange(x3+2, N-1, 1):
            x12 = np.abs(x2 - x1)
            x34 = np.abs(x4 - x3)
            x13 = np.abs(x3 - x1)
            x24 = np.abs(x4 - x2)
            l_1 = list(range(x1, x2+1))
            l_2 = list(range(x3, x4+1))

            eta = (
                np.sin(np.pi * x12 / N) * np.sin(np.pi * x34 / N)
                / (np.sin(np.pi * x13 / N) * np.sin(np.pi * x24 / N))
            )
            counter = 0
            if all(np.abs(eta - stored_eta[-1]) >= 0.0001 for stored_eta in selected_data) and 0 < eta < 1:
                selected_data.append((x1, x2, x3, x4, eta))

    x1, x2, x3, x4, eta = selected_data[a]

    random_number = random.randint(1000, 9999)
    output_file = os.path.join(output_dir, f"mi_N_{N}_D_{
                               D}_eta_{eta}_{random_number}.txt")

    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)

    # Simulation parameters
    t = 4*N  # circuit depth
    n_realizations = 1  # number of realizations
    N_list = list(range(N))
    B = list(range(x1, x2+1))
    C = list(range(x3, x4+1))
    print(print(len(selected_data), x1, x2, x3, x4, eta, B, C))

    int_D = int(D) if D == math.floor(D) else int(math.floor(D) + 1)
    frac_D = D - math.floor(D)

    mi_record = np.zeros(n_realizations)

    start_time = time.perf_counter()
    with open(output_file, 'w') as f:
        for iter in range(n_realizations):
            print(iter)
            tensors = create_tensor_grid(N//2, t)
            tensors_partial = bell_contract_tensors_sqrt(N//2, t, tensors, D)
            result = mi_routine(N, t, tensors_partial)
            mi_record[iter] = bip_mi(result, N, B, C)
            f.write(f"{mi_record[iter]}\n")
            f.flush()

    end_time = time.perf_counter()
    print("total time", end_time - start_time)
