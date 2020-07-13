##########
### Drawing oriented genomes
###

from cgt import *
from sage.all_cmdline import *  # import sage library

### Hyperoctahedral group as signed permutations (cycles)
HO = SignedPermutations(4)
permutation = "(1,-4,-3)(-1,4,3)"

print("Permutation: ", permutation)
permutation = conversions.cycles_to_signed_permutation(4, permutation)
print("...in one-row (region->position): ", list(permutation))
drawing.draw_genome(permutation = list(permutation))


# Now reverse regions in the segment (1,2)
new =  HO(conversions.cycles_to_signed_permutation(4, "(1,-2)(-1,2)")) * HO(list(permutation))
print(new)
drawing.draw_genome(permutation = list(new))

HO = SignedPermutations(5)

# Extra stuff
new = HO([-1,4,3,-5,-2]) * HO([-3,-2,-1,4,5])
print(new)
drawing.draw_genome(permutation = list(new))

# Extra stuff
new = HO([-1,4,3,-5,-2]) * HO(conversions.cycles_to_signed_permutation(5, "(2,-5)(-2,5)(3,-4)(-3,4)"))
print(new)
drawing.draw_genome(permutation = list(new))
