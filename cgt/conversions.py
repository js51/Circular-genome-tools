from sage.all_cmdline import *

def cycles_to_signed_permutation(n, disjoint_cycle_string):
	"""Converts a cycle notation string into two-row notation"""
	g = disjoint_cycle_string
	S = list(range(-n, 0)) + list(range(1, n+1))
	S_N = SymmetricGroup(S)
	g = S_N(g)
	g = [g(i) for i in range(1, n+1)]
	H = SignedPermutations(n)
	return H(g)
	
def matrix_to_signed_permutation(n, M):
	sigma = [None for n in range(n)]
	for i in range(n):
		for j in range(n):
			if M[i][j] != 0:
				sigma[j] = (i+1)*M[i][j]
	return SignedPermutations(n)(sigma)

def signed_permutation_to_cycles(n, signed_permutation, signed=True):
	g = list(signed_permutation)
	if signed:
		S = list(range(-n, 0)) + list(range(1, n+1))
	else:
		S = list(range(1, n+1))
	S_N = SymmetricGroup(S)
	g = [g[i-1] for i in range(1,n+1)]
	if signed:
		g = list(reversed(list(map(lambda x: -1*x, g)))) + g
	return S_N(g)	
	
def signed_cycle_to_unsigned(n, disjoint_cycle_string): # Need to test this
	ans = []
	for a in disjoint_cycle_string.split(')('):
		ans.append((a.strip('(').strip(')')).split(','))
	solution = []
	return ans
