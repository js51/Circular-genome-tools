##########
### Computing irreps for S_n and H_n in Sage/GAP
###

from cgt import *
from sage.all_cmdline import *  # import sage library

### We can compute irreducible representations of the symmetric group using only sage.
##### we don't need to define the group first.
sigma = Permutation("(1,2)(3)(4)") # Just using (1,2) doesn't work.
print(sigma)
irrep = SymmetricGroupRepresentations(4, 'specht')([3,1]) # [3,1] represents the partition and sets n=4
print(irrep.representation_matrix(sigma))

### We can do the same thing with GAP and the package repsn
gap.LoadPackage('"repsn"')
irrep_characters = gap.Irr(SymmetricGroup(4))
irrep = gap.IrreducibleAffordingRepresentation(irrep_characters[2])
m = gap.Image(irrep, sigma)
print(matrix(QQ, m))

### These two matrices may or may not be the same! 
##### But each method gives a system of irreducible representations.

### The second method can be used for general groups (Like the hyperoctahedral group!)

n = 4
# Hyperoctahedral group as signed cycles
S_pmn = SymmetricGroup(list(range(-n, 0)) + list(range(0,n+1)))
HO = SignedPermutations(n)
HO_cycles = S_pmn.subgroup([
	conversions.signed_permutation_to_cycles(n, sigma) for sigma in HO.gens()
])
### Hyperoctahedral group as a wreath product
C2_wr_Sn = gap.WreathProduct(gap.CyclicGroup(2), gap.SymmetricGroup(n)) # Group isomorphic to HO and HO_cycles but in S_{2n}
S2_wr_Sn = gap.WreathProduct(gap.SymmetricGroup(2), gap.SymmetricGroup(n)) # Group isomorphic to HO and HO_cycles but in S_{2n}

# HO as cycles
irrep_characters = gap.Irr(HO_cycles)
irrep = gap.IrreducibleAffordingRepresentation(irrep_characters[13])
image = gap.Image(irrep, HO_cycles("(1,-1)"))
print(image)
#show(matrix(CC, image))

# HO as C2_wr_Sn
irrep_characters = gap.Irr(C2_wr_Sn)
irrep = gap.IrreducibleAffordingRepresentation(irrep_characters[13])
image = gap.Image(irrep, gap.GeneratorsOfGroup(C2_wr_Sn)[1])
print(image)
#show(matrix(CC, image))

# HO as S2_wr_Sn
irrep_characters = gap.Irr(S2_wr_Sn)
irrep = gap.IrreducibleAffordingRepresentation(irrep_characters[13])
image = gap.Image(irrep, gap.GeneratorsOfGroup(S2_wr_Sn)[1])
print(image)
#show(matrix(CC, image))