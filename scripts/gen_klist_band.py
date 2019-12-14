import sys

import numpy as np

# Pearson symbol
# ==============
# see https://en.wikipedia.org/wiki/Pearson_symbol
#
# first letter
# a: triclinic = anorthic
# m: monoclinic
# o: orthorhombic
# t: tetragonal
# h: hexagonal
# c: cubic
#
# second letter
# P: Primitive
# I: Body centred
# F: All faces centred

BZ = dict(
    cP = dict(  # cubic P
        G = [  0,   0,   0],
        X = [1/2,   0,   0],
        M = [1/2, 1/2,   0],
        R = [1/2, 1/2, 1/2],
    ),
    tP = dict(  # tetra P
        G = [  0,   0,   0],
        X = [  0, 1/2,   0],
        M = [1/2, 1/2,   0],
        Z = [  0,   0, 1/2],
        R = [  0, 1/2, 1/2],
        A = [1/2, 1/2, 1/2],
    ),
    # TODO 他の格子のも追加する
)

# convert List to numpy array
for d in BZ.values():
    for k, v in d.items():
        d[k] = np.array(v)

FMT = "{:10}{:5d}{:5d}{:5d}{:5d}{:5.1f}\n"
WEIGHT = 2


def main():
    case, pearson, syms, num = sys.argv[1:]
    num = int(num)

    n = len(syms)
    bz = BZ[pearson]
    n_k = 0

    with open(f"{case}.klist_band", 'w') as f:
        for i in range(n):
            sym = syms[i]
            k = bz[sym]
            i_x = np.round(k * num).astype(int)
            f.write(FMT.format(sym, *i_x, num, WEIGHT))
            n_k += 1

            if i == n - 1:
                break

            k_to = bz[syms[i + 1]]
            v = k_to - k
            m = int(np.round(np.abs(v * num).max()))

            for j in range(1, m):
                i_x = np.round((k + v * j / m) * num).astype(int)
                f.write(FMT.format('', *i_x, num, WEIGHT))
                n_k += 1

        f.write("END\n")

    print(f"generated {n_k} kpoints")
    print(f"file written: {case}.klist_band")


if __name__ == "__main__":
    main()
