import sys
import re
import pickle
from functools import reduce

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


def plot(k, ene, k_ticks, k_labels):
    YLIM = (-10, 10)  # eV

    fig, ax = plt.subplots()
    ax.axhline(0, c='k', ls='--', lw=1)

    ax.scatter(k, ene, s=1, c='C0', marker='.')

    ax.set_xlim((0, k[-1]))
    ax.set_ylim(YLIM)

    ax.set_xticks(k_ticks)
    ax.set_xticklabels(k_labels)

    ax.grid(True, axis='x')

    ax.set_ylabel('Energy [eV]')

    fig.tight_layout()

    return fig


def dump(obj, fname):
    with open(fname, 'wb') as f:
        pickle.dump(obj, f)


def main():
    case = sys.argv[1]

    a = parse_struct(case)  # lattice const.
    n_k, i_ticks, k_labels = parse_klist_band(case)

    print(f"reading {n_k} kpoints")

    # k: n_k x 3 matrix
    # ene: n_k x MAX_N_E matrix
    # n_ene: n_k vector
    k, ene, n_ene = parse_energy(case, n_k, a)

    # get Fermi energy in Ry
    ef = parse_scf(case)
    ene -= ef

    # Ry -> eV
    ry = C.physical_constants['Rydberg constant times hc in eV'][0]
    ene *= ry

    # calc k-path
    k_path = np.zeros(n_k)
    k_from = k[0, :]
    for i in range(1, n_k):
        k_to = k[i, :]
        k_path[i] = k_path[i - 1] + np.linalg.norm(k_to - k_from)
        k_from = k_to
    k_ticks = k_path[i_ticks]

    # flatten
    k = np.concatenate([k_path[i] * np.ones(n_ene[i]) for i in range(n_k)])
    ene = np.concatenate([ene[i, :n_ene[i]] for i in range(n_k)])

    fig = plot(k, ene, k_ticks, k_labels)
    fig.savefig('spaghetti.pdf')
    fig.savefig('spaghetti.png')
    dump(fig, 'spaghetti.pkl')
    print(f"written spaghetti.pdf, spaghetti.png, spaghetti.pkl")

    np.savez_compressed('spaghetti.npz',
                        x=k, y=ene, x_ticks=k_ticks, x_labels=k_labels)
    print(f"written data file: spaghetti.npz")


if __name__ == "__main__":
    main()
