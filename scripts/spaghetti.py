import sys
from io import StringIO
import re

import numpy as np
import scipy.constants as C
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['figure.titlesize'] = 'xx-large'
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'


def parse(fmt, s):
    types = dict(d=int, f=float, E=float, s=str)

    pat = fmt \
        .replace('{:d}', r'\s*(\d+)') \
        .replace('{:f}', r'\s*([+-]?[\d.]+)') \
        .replace('{:E}', r'\s*([+-]?[\d.]+(?:E[+-]\d{2})?)') \
        .replace('{:s}', r'([A-Z]*)\s*')

    return (f(ss) for f, ss in zip(
        [types[t] for t in re.findall(r'\{:(.)\}', fmt)],
        re.match(pat, s).groups()
    ))


def parse_struct(case):
    FMT_A = '{:f} {:f} {:f} {:f} {:f} {:f}'

    with open(f"{case}.struct") as f:
        for _ in range(3):
            f.readline()  # skip
        # FIXME 立方晶系と正方晶系以外にも対応する
        # case.outputirに変換行列が書いてある
        a, b, c, _, _, _ = parse(FMT_A, f.readline())

    return np.array([a, b, c])


def parse_klist_band(case):
    n_k = -1
    k_ticks = []
    k_labels = []

    with open(f"{case}.klist_band") as f:
        for line in f:
            n_k += 1
            s = line[:3]
            if s.strip() and s != 'END':
                k_ticks += [n_k]
                k_labels += s.strip()

    return n_k, k_ticks, k_labels


def parse_energy(case, n_k, a):
    MAX_N_E = 200
    FMT_K = '{:E}{:E}{:E}{:s}{:d}{:d}{:f}'
    FMT_E = '{:d}{:E}'

    k = np.zeros((n_k, 3))
    ene = np.zeros((n_k, MAX_N_E))
    n_ene = np.zeros(n_k, int)

    with open(f"{case}.energy_band") as f:
        for _ in range(8):
            f.readline()  # skip

        for i_k in range(n_k):
            x, y, z, _, _, n_band, _ = parse(FMT_K, f.readline())
            k[i_k, :] = np.array([x, y, z]) / a

            for i in range(n_band):
                _, e = parse(FMT_E, f.readline())
                ene[i_k, i] = e

            n_ene[i_k] = n_band

    return k, ene, n_ene


def parse_scf(case):
    RE_FERMI = re.compile(r'^:FER')

    with open(f'{case}.scf') as f:
        for line in f:
            if RE_FERMI.match(line):
                ef = float(line.split('=')[-1])
    return ef  # get last one value


def parse_spaghetti_ene(case, n_k):
    MAX_I = 200

    # n_k vector
    k = None
    # n_k x bands matrix
    ene = np.zeros((n_k, MAX_I))

    n_band = 0
    with open(f'{case}.spaghetti_ene') as f:
        while True:
            if not f.readline():  # skip:
                break
            with StringIO() as buf:
                for i in range(n_k):
                    buf.write(f.readline() + '\n')
                lines = buf.getvalue()
            with StringIO(lines) as buf:
                a = np.loadtxt(buf)
            k = a[:, 3]
            ene[:, n_band] = a[:, 4]
            n_band += 1

    return k, ene[:, :n_band]


def main():
    case = sys.argv[1]

    # a = parse_struct(case)  # lattice const.
    n_k, i_ticks, k_labels = parse_klist_band(case)

    print(f"reading {n_k} kpoints")

    # k: n_k x 3 matrix
    # ene: n_k x MAX_N_E matrix
    # n_ene: n_k vector
    # k, ene, n_ene = parse_energy(case, n_k, a)

    # get Fermi energy in Ry
    # ef = parse_scf(case)
    # ene -= ef

    # Ry -> eV
    # ry = C.physical_constants['Rydberg constant times hc in eV'][0]
    # ene *= ry

    # k: n_k vector
    # ene: n_k x n_band matrix
    k, ene = parse_spaghetti_ene(case, n_k)

    # calc k-path
    # k_path = np.zeros(n_k)
    # k_from = k[0, :]
    # for i in range(1, n_k):
    #     k_to = k[i, :]
    #     k_path[i] = k_path[i - 1] + np.linalg.norm(k_to - k_from)
    #     k_from = k_to
    # k_ticks = k_path[i_ticks]
    k_ticks = k[i_ticks]

    # flatten
    # k = np.concatenate([k_path[i] * np.ones(n_ene[i]) for i in range(n_k)])
    # ene = np.concatenate([ene[i, :n_ene[i]] for i in range(n_k)])

    np.savez_compressed('spaghetti.npz',
                        x=k, y=ene, x_ticks=k_ticks, x_labels=k_labels)
    print(f"written data file: spaghetti.npz")


if __name__ == "__main__":
    main()
