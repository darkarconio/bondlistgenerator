#!/Users/michalplucinski/Virtualenvs/ase-3.6.0.2515/bin/python

from ase.lattice.spacegroup import crystal

a = 5.431
diamond = crystal('C', [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])

print diamond.positions
