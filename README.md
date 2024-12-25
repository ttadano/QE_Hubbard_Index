# What is this?

This repository contains a simple Python script for generating Hubbard indices for the new DFT+Hubbard input of Quantum ESPRESSO (ver. 7.1 and later).

# Requirements

- ase
- numpy

Please install ase using pip or conda.

# Example

You can get the indices of Hubbard $V_{ij}$ (intersite Coulomb) between atoms $i$ and $j$, for example, by the following command:
```bash
python GetHubbardIndex.py --elements Ni Fe Co  --input scf.in  --cutoff 2.5
```

- `--input`: set an Quantum ESPRESSO input file
- `--elements`: set a list of elements for atomic site $i$
- `--cutoff`: cutoff radius to list the $(i,j)$ pair (in unit of Å)

