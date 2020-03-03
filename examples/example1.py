import cgt

permutation = "(1,-3)(-1,3)"
print(permutation)
permutation = list(cgt.conversions.cycles_to_signed_permutation(4, permutation))
print(permutation)
cgt.drawing.draw_genome(permutation = permutation)