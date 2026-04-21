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


def ee_bipartite(wave, la, l, n=1):
    """Compute the von Neumann or Renyi entropy across a cut of size ``la``."""
    lb = l-la
    temp = np.reshape(wave, (3**la, 3**lb))
    sp = np.linalg.svd(temp, compute_uv=False)
    tol = 1e-20
    sp[abs(sp) < tol] = 0.0
    sp = sp[np.nonzero(sp)]
    el = sp**2
    if n != 1:
        ren = (1 / (1 - n)) * np.log2(np.sum(el**(n)))
        return ren
    else:
        von = -np.dot(el, np.log2(el))
        return von


def ee_routine(N, t, tensors, D, n):
    """Evolve to the final depth and average the entropy across the last step."""
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
            eevst1 = ee_bipartite(obj, N//2, N, n)

        for l2 in eq1:
            obj = gate_action((l2-2) %
                              N, l2-1, tensors[(l2-1)//2, t1+1], obj)
            obj /= np.linalg.norm(obj)

        if t1 == (t-2):
            eevst2 = ee_bipartite(obj, N//2, N, n)
            eevst = (eevst1+eevst2)/2

    return eevst

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run EE Calculation")
    parser.add_argument('--N', type=int, required=True, help="System size N")
    parser.add_argument('--D', type=float, required=True,
                        help="Bond dimension D")
    parser.add_argument('--n', type=float, required=True,
                        help="Renyi index n")
    parser.add_argument('--output', type=str, required=True,
                        help="Output file path")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    N = args.N
    D = args.D
    n = args.n
    output_file = args.output

    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Simulation parameters
    t = 4*N  # circuit depth
    n_realizations = 5000  # number of realizations
    int_D = int(D) if D == math.floor(D) else int(math.floor(D) + 1)
    frac_D = D - math.floor(D)

    ee_record = np.zeros(n_realizations)

    start_time = time.perf_counter()
    with open(output_file, 'w') as f:
        for iter in range(n_realizations):
            print(iter)
            tensors = create_tensor_grid(N//2, t)
            tensors_partial = bell_contract_tensors_sqrt(N//2, t, tensors, D)
            ee_record[iter] = ee_routine(N, t, tensors_partial, D, n)
            f.write(f"{ee_record[iter]}\n")
            f.flush()

    end_time = time.perf_counter()
    print("total time", end_time - start_time)
