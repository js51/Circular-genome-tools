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

