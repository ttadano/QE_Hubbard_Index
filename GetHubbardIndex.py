from ase.io import read
import numpy as np
import argparse


# Set consistent with the QE value (Modules/parameters.f90)
SC_SIZE = 1


def get_supercell_index(structure_base):
    """
    Get the supercell index of the structure
    """
    natom_base = structure_base.get_global_number_of_atoms()
    natom_super = (2 * SC_SIZE + 1) ** 3 * natom_base
    map_p2s = np.zeros((natom_base, (2 * SC_SIZE + 1)**3), dtype=int)
    map_s2p = np.zeros(natom_super, dtype=int)

    xf_super = np.zeros((natom_super, 3), dtype=float)

    xf_super[:natom_base] = structure_base.get_scaled_positions()
    map_p2s[:natom_base, 0] = np.arange(natom_base)
    map_s2p[:natom_base] = np.arange(natom_base)

    iatom = natom_base

    icell = 1
    for nx in range(-SC_SIZE, SC_SIZE + 1):
        for ny in range(-SC_SIZE, SC_SIZE + 1):
            for nz in range(-SC_SIZE, SC_SIZE + 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                for i in range(natom_base):
                    xf_super[iatom] = xf_super[i] + np.array([nx, ny, nz], dtype=float)
                    map_p2s[i, icell] = iatom
                    map_s2p[iatom] = i
                    iatom += 1
                icell += 1

    xc_super = np.dot(xf_super, structure_base.cell)

    return xc_super, map_p2s, map_s2p


def compute_all_distances(xc_super):
    """
    Compute all distances
    """
    natom_super = len(xc_super)
    distances = np.zeros((natom_super, natom_super), dtype=float)

    for i in range(natom_super):
        for j in range(i, natom_super):
            distances[i, j] = np.linalg.norm(xc_super[i] - xc_super[j])
            distances[j, i] = distances[i, j]

    return distances

def print_intersiteV(structure, distances, map_s2p, cutoff, elements):

    natom_base = structure.get_global_number_of_atoms()
    chemical_symbols = structure.get_chemical_symbols()

    for i in range(natom_base):
        if elements is not None:
            if structure.get_chemical_symbols()[i] not in elements:
                continue

        dist_tmp = distances[i, :]
        index_sort = np.argsort(dist_tmp)
        for j in index_sort:
            if dist_tmp[j] > cutoff:
                continue
            symbol_i = chemical_symbols[i]
            symbol_j = chemical_symbols[map_s2p[j]]
            print("V {:4s} {:4s} {:5d} {:5d} {:.3f}".format(symbol_i,
                                                            symbol_j,
                                                            i + 1,
                                                            j + 1,
                                                            dist_tmp[j]))


def get_argoptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('--elements', type=str, nargs='+',
                        help='List of elements', default=None)
    parser.add_argument('--input', type=str, help="structure file",
                        default=None)
    parser.add_argument('--cutoff', type=float, help="cutoff radius",
                        default=5.0)
    return parser.parse_args()


def main():
    argopt = get_argoptions()
    structure = read(argopt.input, format='espresso-in')
    xc_super, _, map_s2p = get_supercell_index(structure)
    distances = compute_all_distances(xc_super)

    print_intersiteV(structure, distances, map_s2p,
                     argopt.cutoff, argopt.elements)


if __name__ == '__main__':
    main()
