##########
### Compute several different versions of the hyperoctahedral group
###

from cgt import *
from sage.all_cmdline import *  # import sage library

### Some constants
n 		= 3 						# number of (signed) letters
N 		= list(range(1, n+1))		# The set {1, ..., n}
NN		= list(range(1, 2*n+1))		# The set {1, ..., 2n}
pmN 	= list(range(-n, 0)) + N 	# The set {+-1, ..., +-n}

### Symmetric groups
## on {1, ..., n}
S_n = SymmetricGroup(N)
## on {1, ..., 2n}
S_2n = SymmetricGroup(NN)
## on {+-1, ..., +-n}
S_pmn = SymmetricGroup(pmN)

### Hyperoctahedral group as a group of signed permutations
## One-row notation
HO = SignedPermutations(n)
## Cycle notation
HO_cycles = S_pmn.subgroup([
	conversions.signed_permutation_to_cycles(n, sigma) for sigma in HO.gens()
])

print(S_2n.subgroup(["(3,4)", "(1,2)", "(5,6)(1,3)(2,4)"]).is_isomorphic(DihedralGroup(4)))

### Hyperoctahedral group as a wreath product
C2_wr_Sn = gap.WreathProduct(gap.CyclicGroup(2), gap.SymmetricGroup(n)) # Group isomorphic to HO and HO_cycles but in S_{2n}

### Hyperoctahedral group in S_2n (As converted from wreath product by GAP)
HO_in_S2n = gap.WreathProduct(SymmetricGroup(2), SymmetricGroup(n)) # Group isomorphic to HO and HO_cycles but in S_{2n}
HO_in_S2n = PermutationGroup(gap_group = HO_in_S2n)
print(HO_cycles.is_isomorphic(HO_in_S2n))

### Hyperoctahedral group as a matrix group (canonical representation)
OnZ = HO.matrix_group()

### We can check that each of these is in fact the hyperoctahedral group
print(HO_cycles.is_isomorphic(HO_in_S2n) and OnZ.is_isomorphic(HO_cycles))

##########
### Dihedral group in the hyperoctahedral group
###

### Dihedral group (cycle notation)
r, s = rotation_permutation(n), reflection_permutation(n) # Generators for D_n
Dn_cycles = HO_cycles.subgroup([s,r])
###

### Dihedral group (as a matrix group in O_n(Z))
s = conversions.cycles_to_signed_permutation(n, s).to_matrix()
r = conversions.cycles_to_signed_permutation(n, r).to_matrix()
Dn_in_OnZ = OnZ.subgroup([s, r])

### These are of course isomorphic as well
print(Dn_in_OnZ.is_isomorphic(Dn_cycles))