from cgt import *
from sage.all_cmdline import *  # import sage library

n = 3
S = list(range(-n, 0)) + list(range(1, n+1)) # The set {+-1, ..., +-n}

SS = SymmetricGroup(S) 	# The symmetric group on S isomorphic to S_n
HO 	= SignedPermutations(n) # Signed permutations of {1,...,n}, a subgroup of S_S above
HO_cycles = SS.subgroup([conversions.signed_permutation_to_cycles(n, sigma) for sigma in HO.gens()]) # Above group HO as a subgroup of S_S above (cycle notation)
S2_W_Sn = PermutationGroup(gap_group = gap.WreathProduct(SymmetricGroup(2), SymmetricGroup(n))) # Group isomorphic to HO and HO_cycles but in S_{2n}

gens = []
for sigma in HO_cycles.list():
	gens.append(conversions.signed_cycle_to_unsigned(n, str(sigma)))
print(HO_cycles.list())
print(gens)
print(S2_W_Sn.list())

# Constructing D_n in HO_cycles
r, s = rotation_permutation(n), reflection_permutation(n) # Generators for D_n
Dn_cycles = HO_cycles.subgroup([s,r])
###

# Constructing D_n as a subset of HO (signed permutations)
s = conversions.cycles_to_signed_permutation(n, s).to_matrix()
r = conversions.cycles_to_signed_permutation(n, r).to_matrix()
Dn = HO.matrix_group().subgroup([s, r])
Dn = {conversions.matrix_to_signed_permutation(n, sigma.list()) for sigma in Dn.list() }
###