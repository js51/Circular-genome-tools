from sage.all_cmdline import *
from cgt import hyperoctahedral_groups
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

def to_pos_paradigm(block_graph, blocks, G):
	# Position paradigm (for num regions n):	
	perm = []
	for b in block_graph[::2]: #every second element
		if b>0:
			perm += list(blocks[b-1])
		elif b<0:
			perm += list(reversed(list(-x for x in blocks[-b-1])))
	return G(str(tuple(perm)) + str(tuple(reversed(list(-p for p in perm)))))

def rearrangements_with_cuts_set(cuts_set, n):
	pmN = list(range(-n, 0)) + list(range(1, n+1))	# The set {+-1, ..., +-n}
	S_pmn = SymmetricGroup(pmN)
	blocks = []
	c = c_permutation(n, S_pmn)
	cs = cuts_set
	for i in range(len(cs)):
		blocks.append(tuple(j for j in c_range(c(cs[i]), cs[(i+1)%len(cs)]+1, n, S_pmn)))
	m = len(blocks)
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
						yield to_pos_paradigm(perm_copy, blocks, S_pmn)
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
		
	

cs = [1,2,3,4,5]
print(len(list(rearrangements_with_cuts_set(cs, 1000))))