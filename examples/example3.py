from cgt import *
from sage.all_cmdline import *  # import sage library

n = 4
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

NEW_REARRANGEMENTS = H_CYCLES.subgroup([
	H_CYCLES("(1,-2)(-1,2)"), 
	H_CYCLES("(2,-3)(-2,3)"),
	H_CYCLES("(3,-4)(-3,4)"),
	H_CYCLES("(1,-1)"),
	H_CYCLES("(2,-2)"),
	H_CYCLES("(3,-3)"),
	H_CYCLES("(4,-4)"),
])

#print(NEW_REARRANGEMENTS.order())

#print(D.is_normal(H_CYCLES))
show(H_CYCLES.cosets(D))


def print_stuff():
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
		SETS_OF_CLASSES = {
			#frozenset({ d_1*sigma*d_2 for d_1 in D for d_2 in D}) for sigma in H_CYCLES
		}
		SETS_OF_CLASSES_2 = set()
		#for sigma in H_CYCLES:
		#	A = set()
		#	B = set()
		#	for d_1 in D:
		#		for d_2 in D:
		#			A.add(d_1*sigma*d_2)
		#			B.add(d_1*sigma.inverse()*d_2)
		#	SETS_OF_CLASSES_2.add(frozenset.union(frozenset(A), frozenset(B)))
		# print(H_CYCLES.character_table())
		print(f'N={n}')
		print("group:\torder:")
		print(f'H\t{H.order()}')
		print(f'H_CY\t{H_CYCLES.order()}')
		print(f'D\t{D.order()}')
		print(f'Equiv. Classes:\t{len(SETS_OF_CLASSES)}')
		print(f'Equiv. Classes:\t{len(SETS_OF_CLASSES_2)}')
	
print_stuff()
