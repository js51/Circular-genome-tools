from cgt import *
from sage.all_cmdline import *  # import sage library

n = 5
s = reflection_permutation(n)
r = rotation_permutation(n)
H = SignedPermutations(n)
S = list(range(-1*n, 0)) + list(range(1, n+1)) # The set {+-1, ..., +-n}
S_S = SymmetricGroup(S) 	# The symmetric group on S isomorphic to S_n
H 	= SignedPermutations(n) # Signed permutations of {1,...,n}, a subgroup of S_S above
H_CYCLES = S_S.subgroup([
	conversions.signed_permutation_to_cycles(n, sigma) for sigma in H.gens()
]) # Above group HO as a subgroup of S_S above (cycle notation)
D	= H_CYCLES.subgroup([r,s])


for n in range(2,10):
	s = reflection_permutation(n)
	r = rotation_permutation(n)
	H = SignedPermutations(n)
	S = list(range(-1*n, 0)) + list(range(1, n+1)) # The set {+-1, ..., +-n}
	S_S = SymmetricGroup(S) 	# The symmetric group on S isomorphic to S_n
	H 	= SignedPermutations(n) # Signed permutations of {1,...,n}, a subgroup of S_S above
	H_CYCLES = S_S.subgroup([
		conversions.signed_permutation_to_cycles(n, sigma) for sigma in H.gens()
	]) # Above group HO as a subgroup of S_S above (cycle notation)
	D	= H_CYCLES.subgroup([r,s])
	SETS_OF_CLASSES_2 = {
		frozenset.union(
			frozenset({ d_1*sigma*d_2 for d_1 in D for d_2 in D}),
			frozenset({ d_1*sigma.inverse()*d_2 for d_1 in D for d_2 in D})
		) 
		for sigma in H_CYCLES
	}
	
	
	# print(H_CYCLES.character_table())
	
	print(f'N={n}')
	print("group:\torder:")
	print(f'H\t{H.order()}')
	print(f'H_CY\t{H_CYCLES.order()}')
	print(f'D\t{D.order()}')
	print(f'Equiv. Classes:\t{len(SETS_OF_CLASSES)}')
	print(f'Equiv. Classes:\t{len(SETS_OF_CLASSES_2)}')


