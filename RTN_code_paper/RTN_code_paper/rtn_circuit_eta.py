import numpy as np
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


c = np.exp(2j * np.pi / 3)
bell_3 = (1/np.sqrt(3))*np.transpose(np.array([[1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 1, 0, 0, 0, 1, 1, 0, 0], [0, 0, 1, 1, 0, 0, 0, 1, 0], [1, 0, 0, 0, c, 0, 0, 0, c**2], [
    0, 1, 0, 0, 0, c, c**2, 0, 0], [0, 0, 1, c, 0, 0, 0, c**2, 0], [1, 0, 0, 0, c**2, 0, 0, 0, c], [0, 1, 0, 0, 0, c**2, c, 0, 0], [0, 0, 1, c**2, 0, 0, 0, c, 0]]))
bell_qubit_qutrit = (1/np.sqrt(2))*np.array([[1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [
    0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
bell_gate_3 = np.reshape(bell_3, (3, 3, 3, 3))


def create_bell_tensor(D):
    d_frac = ((D-1)/2)**1
    bell_final = np.zeros((3, 3))
    bell_final[0, 0] = (1-d_frac) + d_frac/np.sqrt(3)
    bell_final[1, 1] = d_frac/np.sqrt(3)
    bell_final[2, 2] = d_frac/np.sqrt(3)
    return bell_final/np.linalg.norm(bell_final)


def create_bell_tensor_sqrt(D):
    d_frac = ((D-1)/2)**1
    bell_final = np.zeros((3, 3))
    bell_final[0, 0] = (1-d_frac) + d_frac/np.sqrt(3)
    bell_final[1, 1] = d_frac/np.sqrt(3)
    bell_final[2, 2] = d_frac/np.sqrt(3)
    return np.sqrt(bell_final/np.linalg.norm(bell_final))


def gate_action(control, target, gate, reg):
    """Apply a two-site gate to the control and target indices of ``reg``."""
    result = np.tensordot(gate, reg, ((2, 3), (control, target)))
    result = np.moveaxis(result, (0, 1), (control, target))
    return result


def create_random_tensor(n_legs):
    return np.random.randn(*[3 for _ in range(n_legs)])


def create_tensor_grid(N, t):
    """Return an ``N x t`` grid of random rank-4 qutrit tensors."""
    tensor_grid = np.empty((N, t), dtype=object)
    for j in range(t):
        for i in range(N):
            tensor_grid[i, j] = create_random_tensor(4)
    return tensor_grid


def ee_bipartite(wave, la, l):
    """Compute the base-3 bipartite entropy across a cut of size ``la``."""
    lb = l-la
    temp = np.reshape(wave, (3**la, 3**lb))
    sp = np.linalg.svd(temp, compute_uv=False)
    sp = sp[np.nonzero(sp)]
    el = sp**2
    von = -np.dot(el, (np.log2(el)))
    von *= 1/np.log2(3)
    return von


def ee_general(psi, sites, N):
    psi = np.moveaxis(psi, sites, np.arange(len(sites)))
    ee = ee_bipartite(psi, len(sites), N)
    return ee


def bip_mi(psi, sites, N):
    """Compute the mutual information between the two sites."""
    B = [sites[0]]
    C = [sites[1]]
    S_B = ee_general(psi, B, N)
    S_C = ee_general(psi, C, N)
    BC = np.concatenate([B, C])
    S_BC = ee_general(psi, BC, N)
    bip_mi = S_B + S_C - S_BC
    return bip_mi


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


projectorsbell3 = [
    np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
    np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]]),
    np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]]),
]


def mutual_info_routine(N, t, steady_t, tensors, D, x1, x2):
    """Insert two ancillas after reaching steady state and track their mutual information."""
    choice1 = random.choice([0, 1, 2])
    choice2 = random.choice([0, 1, 2])
    projector1 = projectorsbell3[choice1]
    projector2 = projectorsbell3[choice2]
    ancillas = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]])

    obj = np.zeros((3,) * N)
    obj[(0,)*N] = 1
    nq1 = range(N)
    oq1 = nq1[::2]
    eq1 = nq1[1::2]

    for t1 in range(0, steady_t, 2):
        for l1 in oq1:
            obj = gate_action(l1, l1+1, tensors[l1//2, t1], obj)
            obj /= np.linalg.norm(obj)

        for l2 in eq1:
            obj = gate_action((l2-2) % N, l2-1, tensors[(l2-1)//2, t1+1], obj)
            obj /= np.linalg.norm(obj)

    obj = np.moveaxis(np.tensordot(projector1, obj, (1, x1)), 0, x1)
    obj = np.moveaxis(np.tensordot(projector2, obj, (1, x2)), 0, x2)
    obj = np.tensordot(obj, ancillas, axes=0)
    obj = gate_action(x1, N, bell_gate_3, obj)
    obj = gate_action(x2, N+1, bell_gate_3, obj)
    obj /= np.linalg.norm(obj)
    mutual_info = np.zeros((t-steady_t)+1)
    mutual_info[0] = bip_mi(obj, [N, N+1], N+2)
    print("mutual info 0 is ", mutual_info[0])
    for t1 in range(steady_t, t, 2):
        for l1 in oq1:
            obj = gate_action(l1, l1+1, tensors[l1//2, t1], obj)
            obj /= np.linalg.norm(obj)

        mutual_info[t1-steady_t+1] = bip_mi(obj, [N, N+1], N+2)
        for l2 in eq1:
            obj = gate_action((l2-2) % N, l2-1, tensors[(l2-1)//2, t1+1], obj)
            obj /= np.linalg.norm(obj)
        mutual_info[t1-steady_t+2] = bip_mi(obj, [N, N+1], N+2)
    return mutual_info


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run Dynamical Exp Calculation")
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
        os.makedirs(output_dir, exist_ok=True)

    # Simulation parameters
    t = 10*N  # circuit depth
    steady_t = 4*N  # circuit depth to steady state
    n_realizations = 2000  # number of realizations

    x1 = 1
    x2 = x1 + N//2

    start_time = time.perf_counter()

    for iter in range(n_realizations):
        print(iter)
        tensors = create_tensor_grid(N//2, t)
        partial_tensors = bell_contract_tensors_sqrt(N//2, t, tensors, D)
        mutual_inf_data = mutual_info_routine(
            N, t, steady_t, partial_tensors, D, x1, x2)
        base_name = os.path.splitext(os.path.basename(output_file))[0]
        ext = os.path.splitext(output_file)[1]
        iteration_output_file = os.path.join(
            output_dir,
            f"{base_name}_{iter}{ext}"
        )

        with open(iteration_output_file, 'w') as f:
            for val in mutual_inf_data:
                f.write(f"{val}\n")

    end_time = time.perf_counter()
    print("total time", end_time - start_time)
