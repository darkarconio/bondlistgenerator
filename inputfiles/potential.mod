# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential
bond_style harmonic
bond_coeff 1 5.15667 2.35169198
angle_style harmonic 
angle_coeff 1 1.171 109.5

# Setup neighbor style
atom_modify sort 0 2
communicate single cutoff 2
#neighbor 1.0 nsq
#neigh_modify once no every 1 delay 0 check yes 

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no

