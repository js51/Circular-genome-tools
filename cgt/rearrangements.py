"""
"""

from enum import Enum, auto
from copy import deepcopy
from sage.all_cmdline import *
from . import structures
from .enums import *

def css(cuts_set):
	string = ''.join(sorted([str(cut[0]) for cut in cuts_set]))
	if len(string) != len(cuts_set): print("something went very wrong here...")
	return string

def cuts(framework, sigma):
	sigma = deepcopy(framework.one_row(sigma))
	n = framework.n
	pmN = structures.set_plus_minus(n)
	G = SymmetricGroup(pmN)
	c_string = f'({",".join(str(i) for i in list(range(1, n+1)))})({",".join(str(i) for i in list(range(-n, 0)))})'
	c = G(c_string)
	sigma = list(sigma)
	return set([ 
		tuple(([i+1,c(i+1)])) for i in range(len(sigma)) 
		if (sigma[c(i + 1)-1] != c(sigma[i]))
	])

def __all_canonical_inversions(framework, num_regions=None):
	if not framework.oriented:
		raise NotImplementedError(f"not yet implemented for {str(framework)}")
	if framework.symmetry not in {SYMMETRY.circular, SYMMETRY.linear}:
		raise NotImplementedError(f"not implemented for framework with {str(framework.symmetry)} symmetry")
	n = framework.n
	G = framework.genome_group()
	if num_regions == 1:
		return {G(f'({i},-{i})') for i in range(1, n+1)} if framework.oriented else {}
	elif num_regions == 2:
		if framework.oriented:
			perms = {G(f'({i},{-1*(i+1)})({-1*i},{i+1})') for i in range(1,n)}
			if framework.symmetry is not SYMMETRY.linear:
				perms.union({G(f'({n},-1)(-{n},1)')})
		else:
			perms = {G(f'({i},{i+1})') for i in range(1,n)} 
			if framework.symmetry is not SYMMETRY.linear:
				perms.union({G(f'(1,{n})')})
		return perms
	elif num_regions == None: # Return all inversions up to length floor(n/2)
		up_to_length = floor(framework.n/2)
		if framework.oriented: 
			perms = set()
			for permutation in G:
				cycle_type = list(permutation.cycle_type())
				if set(cycle_type)=={1,2} and list(cycle_type).count(2)<=up_to_length and len(cuts(framework, permutation))==2:
					perms.add(permutation)
			if framework.symmetry is SYMMETRY.linear:
				perms = {perm for perm in perms if not ('1' in str(perm) and str(n) in str(perm))}
			return perms
		else:
			raise NotImplementedError(f"not yet implemented for {framework}")
	else:
		raise NotImplementedError(f"inversions of length {num_regions} not yet implemented")
	
def __one_region_adjacent_transposition_reps(framework):
	warnings.warn("This function is currently untested! Generators might be incorrect")
	if not framework.oriented or framework.symmetry != SYMMETRY.circular:
		raise NotImplementedError(f"not yet implemented for {str(framework)}")
	G = framework.genome_group()
	return { G('(-2,-1)(1,2)'), G('(-2,1,2,-1)'), G('(-2,-1,2,1)'), G('(-2,2)(-1,1)') }

def __representatives(framework, set_of_permutations, classes=CLASSES.double_cosets):
	if not framework.oriented:
		raise NotImplementedError(f"not yet implemented for {str(framework)}")
	Z = framework.symmetry_group()
	perms = {x for x in set_of_permutations}
	cosets = []
	while len(perms):
		perm = perms.pop()
		if classes is CLASSES.double_cosets:
			coset = { d1 * perm * d2 for d1 in Z for d2 in Z }
		elif classes is CLASSES.conjugacy_classes:
			coset = { d.inverse() * perm * d for d in Z }
		elif classes is CLASSES.cosets:
			coset = { perm * d for d in Z }
		for element in coset:
			perms.discard(element)
		cosets.append(coset)
	reps = set()
	for coset in cosets:
		reps.add(sorted(list(coset), key=lambda x: (sum(set(x.cycle_type())), list(x.cycle_type()).count(2), str(x)) )[0])
	return reps

def all_inversions_representatives(framework, num_regions=None):
	return __representatives(framework, __all_canonical_inversions(framework, num_regions=num_regions), classes=CLASSES.double_cosets)