from sage.all_cmdline import *
from cgt import hyperoctahedral_groups, conversions
import numpy as np
import itertools


def c_permutation(n, G):
	c_string = f'({",".join(str(i) for i in list(range(1, n+1)))})({",".join(str(i) for i in list(range(-n, 0)))})'
	c = G(c_string)
	return c

def c_range(x,y,n,G):
	c = c_permutation(n,G)
	i=x
	while i != y:
		yield i
		i = c(i)

def rearrangements_with_k_cuts(k,n):
	N = (i for i in range(1,n+1))
	for pattern in itertools.combinations(N,k):
		yield from rearrangements_with_cuts_set(pattern, n)
		
def to_pos_paradigm(block_graph, blocks, cuts_set, G):
	H_c = hyperoctahedral_groups.HyperoctahedralGroup(n)
	H = hyperoctahedral_groups.HyperoctahedralGroup(n, hyperoctahedral_groups.SET.ONE_ROW)
	D = hyperoctahedral_groups.DihedralSubgroup(H_c, n)
	perm = []
	c = c_permutation(max(max(blocks)), G)
	for b in block_graph[1::2]: #every second element
		if b>0:
			perm += list(blocks[b-1])
		elif b<0:
			perm += list(reversed(list(-x for x in blocks[-b-1])))
	return perm

def rearrangements_with_cuts_set(cuts_set, n):
	pmN = list(range(-n, 0)) + list(range(1, n+1))	# The set {+-1, ..., +-n}
	S_pmn = SymmetricGroup(pmN)
	blocks = []
	c = c_permutation(n, S_pmn)
	cs = cuts_set
	for i in range(len(cs)):
		blocks.append(tuple(j for j in c_range(c(cs[i]), c(cs[(i+1)%len(cs)]), n, S_pmn)))
	m = len(blocks)
	for b in blocks:
		if 2 in b:
			ind = blocks.index(b)
	blocks = blocks[ind:] + blocks[:ind]
	print(f"cuts set {cs} leads to blocks {blocks}")
	pmM = list(range(-m, 0)) + list(range(1, m+1))	# The set {+-1, ..., +-n}
	S_pmm = SymmetricGroup(pmM)
	nat_ord = S_pmm(sum(((-i,i) for i in range(1,m+1)), tuple()))
	perm = [-1,1]
	def _next_block(perm):
		for option in range(-m,m+1):
			if option not in perm+[0] and option != nat_ord(perm[-1]) and nat_ord(option) != perm[-1]:
				if len(perm)+1 < 2*m-2 or nat_ord(-option) != perm[0]:
					perm_copy = perm.copy()
					perm_copy += [option, -option]
					if len(perm_copy) == 2*m:
						print(perm_copy)
						yield to_pos_paradigm(perm_copy, blocks, cuts_set, S_pmn)
					else:
						yield from _next_block(perm_copy)
	yield from _next_block(perm)
			
	
#	
	# Choose orientations
		#	orientations = keywords = [i for i in itertools.product([1,-1], repeat = len(blocks))]
	
	
		#	for orientation in orientation:
#		#for b,block in enumerate(blocks):
		#	o = orientation[b]
		#	# Decide the orientation, add to the perm
		#	_block = block[::-1] if o=-1 else block
		#	p = tuple(o*r for r in _block))
			
			
	# 	Choose a following block (if current one is not inverted, choose a non-adjacent one)
	# Convert to position paradigm.
		
	
n=6
pmN = list(range(-n, 0)) + list(range(1, n+1))	# The set {+-1, ..., +-n}
S_pmn = SymmetricGroup(pmN)

def cuts(sigma):
	c = c_permutation(n, S_pmn)
	sigma = list(sigma)
	return set([ 
		tuple(([i+1,c(i+1)])) for i in range(len(sigma)) 
		if (sigma[c(i + 1)-1] != c(sigma[i]))
	])


cs = [1,3,4]
rearrangements = list(rearrangements_with_cuts_set(cs,n))
print([cuts(a) for a in rearrangements])

rearrangements = list(rearrangements_with_k_cuts(3, n))
H = hyperoctahedral_groups.HyperoctahedralGroup(n)
D = hyperoctahedral_groups.DihedralSubgroup(H, n)
rearrangements_cycles = [H(conversions.signed_permutation_to_cycles(n, a)) for a in rearrangements]


rearrangements_cycles = {a for a in list(H) if len(cuts(conversions.cycles_to_signed_permutation(n, a))) == 3}
cosets = {
	frozenset([d_1*a*d_2 for d_1 in D for d_2 in D]) for a in rearrangements_cycles
}
for coset in cosets:
	print({frozenset(cuts(conversions.cycles_to_signed_permutation(n, x))) for x in coset})
