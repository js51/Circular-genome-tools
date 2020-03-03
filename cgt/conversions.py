from sage.all_cmdline import *  # import sage library

def cycles_to_signed_permutation(n, disjoint_cycle_string):
	"""Converts a cycle notation string into two-row notation"""
	g = disjoint_cycle_string
	S = list(range(-n, 0)) + list(range(1, n+1))
	S_N = SymmetricGroup(S)
	g = S_N(g)
	g = [g(i) for i in range(1, n+1)]
	H = SignedPermutations(n)
	return H(g)
	
def signed_permutation_to_cycles(n, signed_permutation):
	g = list(signed_permutation)
	S = list(range(-n, 0)) + list(range(1, n+1))
	S_N = SymmetricGroup(S)
	g = [g[i-1] for i in range(1,n+1)]
	g = list(reversed(list(map(lambda x: -1*x, g)))) + g
	return S_N(g)