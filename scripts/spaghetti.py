import sys
from io import StringIO
import re

import numpy as np
import matplotlib as mpl

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


FMT_A = '{:f} {:f} {:f} {:f} {:f} {:f}'
RE_ATOM = re.compile(r'^ATOM\s+')
RE_ATOM_X = re.compile(r'^([A-Z][a-z]?)\s+')
def parse_struct(case):
    with open(f"{case}.struct") as f:
        for _ in range(3):
            f.readline()  # skip
        # FIXME 立方晶系と正方晶系以外にも対応する
        # case.outputirに変換行列が書いてある
        a, b, c, _, _, _ = parse(FMT_A, f.readline())

    atoms = []
    with open(f'{case}.struct') as f:
        b = False
        for line in f:
            if RE_ATOM.match(line):
                b = True
            m = RE_ATOM_X.match(line)
            if b and m:
                atoms += [m[1]]
                b = False

    return np.array([a, b, c]), atoms


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


MAX_N_E = 200
FMT_K = '{:E}{:E}{:E}{:s}{:d}{:d}{:f}'
FMT_E = '{:d}{:E}'
def parse_energy(case, n_k, a):
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


RE_FERMI = re.compile(r'^:FER')
def parse_scf(case):
    with open(f'{case}.scf') as f:
        for line in f:
            if RE_FERMI.match(line):
                ef = float(line.split('=')[-1])
    return ef  # get last one value


MAX_I = 200
def parse_spaghetti_ene(case, n_k):
    # n_k vector
    k = None
    # n_bands x n_k matrix
    ene = np.zeros((MAX_I, n_k))

    n_band = 0
    with open(f'{case}.spaghetti_ene') as f:
        while True:
            if not f.readline():  # skip:
                break
            with StringIO() as buf:
                for i in range(n_k):
                    buf.write(f.readline() + '\n')
                buf.seek(0)
                a = np.loadtxt(buf)
            k, ene[n_band, :] = a[:, 3], a[:, 4]
            n_band += 1

    return k, ene[:n_band, :]


N_ORB = 11
RE_NAT = re.compile(r'NAT=\s*(\d+)')
LABELS_ORB = ['total', 's', 'p', 'pz', 'px+py',
              'd', 'dz2', 'dxy', 'dx2y2', 'dxz+dyz', 'f']
def parse_qtl(case, atoms, n_bands, n_k):
    with open(f'{case}.qtl') as f:
        for _ in range(3):
            f.readline()  # skip

        n_atoms = int(re.search(RE_NAT, f.readline())[1])
        if n_atoms != len(atoms):
            print(f"ERROR invalid number of atoms: {n_atoms} != len({atoms})")
            os.exit(1)

        dtype = [
            (atom, [
                (orb, float)
                for orb in LABELS_ORB
            ])
            for atom in atoms
        ]
        chg = np.zeros((n_atoms, N_ORB, n_bands, n_k), dtype=dtype)

        for _ in range(n_atoms):
            f.readline()  # skip
        for i_band in range(n_bands):
            f.readline()  # skip
            with StringIO() as buf:
                for _ in range(n_k):
                    for j in range(n_atoms):
                        buf.write(f.readline() + '\n')
                    f.readline()  # skip
                buf.seek(0)
                a = np.loadtxt(buf)
            a = a[:, 2:].reshape((n_k, n_atoms, N_ORB)).transpose((1, 2, 0))
            chg[:, :, i_band, :] = a

    return chg


def main(case):
    _, atoms = parse_struct(case)  # lattice const.
    print(f"reading {len(atoms)} atoms")

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
    # ene: n_band x n_k matrix
    k, ene = parse_spaghetti_ene(case, n_k)
    n_bands = ene.shape[0]

    # calc k-path
    # k_path = np.zeros(n_k)
    # k_from = k[0, :]
    # for i in range(1, n_k):
    #     k_to = k[i, :]
    #     k_path[i] = k_path[i - 1] + np.linalg.norm(k_to - k_from)
    #     k_from = k_to
    # k_ticks = k_path[i_ticks]
    k_ticks = k[i_ticks]

    chg = parse_qtl(case, atoms, n_bands, n_k)

    # flatten
    # k = np.concatenate([k_path[i] * np.ones(n_ene[i]) for i in range(n_k)])
    # ene = np.concatenate([ene[i, :n_ene[i]] for i in range(n_k)])

    np.savez_compressed('spaghetti.npz',
                        k=k, ene=ene, chg=chg,
                        k_ticks=k_ticks, k_labels=k_labels)
    print(f"written data file: spaghetti.npz")


if __name__ == "__main__":
    main(sys.argv[1])
