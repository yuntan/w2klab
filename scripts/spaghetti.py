# input: case.struct, case.klist_band, case.spaghetti_ene, case.qtl
# output: spaghetti.npz
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


# const for parse_struct
FMT_A = '{:f} {:f} {:f} {:f} {:f} {:f}'
RE_ATOM = re.compile(r'^ATOM\s+')
RE_ATOM_X = re.compile(r'^([A-Z][a-z]?)\s+')

# const for parse_energy
MAX_N_E = 200
FMT_K = '{:E}{:E}{:E}{:s}{:d}{:d}{:f}'
FMT_E = '{:d}{:E}'

# const for parse_scf
RE_FERMI = re.compile(r'^:FER')

# const for parse_spaghetti_ene
RE_BANDINDEX = re.compile(r'^\s*bandindex:')

# const for parse_qtl
N_ORBS = 11
RE_NAT = re.compile(r'NAT=\s*(\d+)')
RE_BAND = re.compile(r'^\s*BAND')
ORBS = ['total', 's', 'p', 'pz', 'px+py',
        'd', 'dz2', 'dxy', 'dx2y2', 'dxz+dyz', 'f']


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


def parse_scf(case):
    with open(f'{case}.scf') as f:
        for line in f:
            if RE_FERMI.match(line):
                ef = float(line.split('=')[-1])
    return ef  # get last one value


def parse_spaghetti_ene(case, n_k):
    with open(f'{case}.spaghetti_ene') as f:
        n_bands = sum(1 for line in f if RE_BANDINDEX.match(line))

    # n_k vector
    k = None
    # n_bands x n_k matrix
    ene = np.zeros((n_bands, n_k))

    with open(f'{case}.spaghetti_ene') as f:
        for i_band in range(n_bands):
            f.readline()  # skip
            with StringIO() as buf:
                for i in range(n_k):
                    buf.write(f.readline() + '\n')
                buf.seek(0)
                a = np.loadtxt(buf)
            k, ene[i_band, :] = a[:, 3], a[:, 4]

    return k, ene


def parse_qtl(case, atoms, n_k):
    n_atoms = len(atoms)

    with open(f'{case}.qtl') as f:
        n_bands = sum(1 for line in f if RE_BAND.match(line))

    ene = np.zeros((n_bands, n_k))
    dt = np.dtype([(atom, [(orb, float) for orb in ORBS]) for atom in atoms])
    chg = np.zeros((n_bands, n_k), dtype=dt)

    with open(f'{case}.qtl') as f:
        [f.readline() for _ in range(4 + n_atoms)]  # skip

        for i_band in range(n_bands):
            f.readline()  # skip

            with StringIO() as buf:
                for _ in range(n_k):
                    for _ in range(n_atoms):
                        buf.write(f.readline() + '\n')
                    f.readline()  # skip
                buf.seek(0)
                a = np.loadtxt(buf)

            ene[i_band, :] = a[::n_atoms, 0]
            for i_atom in range(n_atoms):
                for i_orb in range(N_ORBS):
                    atom = atoms[i_atom]
                    orb = ORBS[i_orb]
                    chg[atom][orb][i_band, :] = a[i_atom::n_atoms, i_orb + 2]

    return ene, chg


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
    k_ticks = k[i_ticks]
    print(f"reading {ene.shape[0]} bands")

    # calc k-path
    # k_path = np.zeros(n_k)
    # k_from = k[0, :]
    # for i in range(1, n_k):
    #     k_to = k[i, :]
    #     k_path[i] = k_path[i - 1] + np.linalg.norm(k_to - k_from)
    #     k_from = k_to
    # k_ticks = k_path[i_ticks]

    _, chg = parse_qtl(case, atoms, n_k)

    # flatten
    # k = np.concatenate([k_path[i] * np.ones(n_ene[i]) for i in range(n_k)])
    # ene = np.concatenate([ene[i, :n_ene[i]] for i in range(n_k)])

    np.savez_compressed('spaghetti.npz',
                        k=k, ene=ene, chg=chg,
                        k_ticks=k_ticks, k_labels=k_labels)
    print(f"written data file: spaghetti.npz")


if __name__ == "__main__":
    main(sys.argv[1])
