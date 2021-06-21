from enum import Enum, auto
from sage.all_cmdline import *
from . import conversions
from . import helper_functions
from .enums import *

def HyperoctahedralGroup(n, as_set_of=SET.SIGNED_CYCLES):
	### Some constants
	N 		= list(range(1, n+1))		# The set {1, ..., n}
	NN		= list(range(1, 2*n+1))		# The set {1, ..., 2n}
	pmN 	= list(range(-n, 0)) + N 	# The set {+-1, ..., +-n}
	### Symmetric groups
	S_n = SymmetricGroup(N)
	S_2n = SymmetricGroup(NN)
	S_pmn = SymmetricGroup(pmN)

	if as_set_of in [SET.SIGNED_CYCLES, SET.ONE_ROW]:
		HO = SignedPermutations(n)
		if as_set_of == SET.ONE_ROW:
			return HO
		else:
			HO_cycles = S_pmn.subgroup([
				conversions.signed_permutation_to_cycles(n, sigma) for sigma in HO.gens()
			])
			HO_cycles._as_set_of = SET.SIGNED_CYCLES
			HO_cycles._n = n
			return HO_cycles
	else:
		return "Not yet implemented"

def DihedralSubgroup(G, n=None):
	try:
		if G._as_set_of == SET.SIGNED_CYCLES:
			n = G._n
			r, s = helper_functions.rotation_permutation(n), helper_functions.reflection_permutation(n) # Generators for D_n
			D_n = G.subgroup([s,r])
			D_n._as_set_of = SET.SIGNED_CYCLES
			D_n._n = G._n
			return D_n
	except:
		print('Assuming group is S_n')
		return DihedralGroup(n)
		

def EquivalenceClasses(G, n=None, symmetry_group='dihedral', classes='double-cosets-and-inverses', classes_as="counts", sort_classes=False):
	'''Example use: EquivalenceClasses(HyperoctahedralGroup(5, as_set_of=SET.SIGNED_CYCLES))'''
	H_n = G
	try:
		n = H_n._n
	except:
		print("Assuming group is S_n")
	if symmetry_group == 'dihedral':
		D_n = DihedralSubgroup(H_n, n)
	else:
		D_n = symmetry_group
	cards = {} if classes_as == 'dict' else []
	H_n_elements = Set(H_n)
	while H_n_elements.cardinality()>0:
		if classes == 'genomes':
			T = Set([H_n_elements[0]*d for d in D_n])
		if classes == 'double-cosets':
			T = Set([d_1*H_n_elements[0]*d_2 for d_1 in D_n for d_2 in D_n])
		elif classes == 'double-cosets-and-inverses':
			T = Set([d_1*H_n_elements[0]*d_2 for d_1 in D_n for d_2 in D_n])
			T = T.union(Set([g.inverse() for g in T]))
		H_n_elements = H_n_elements.difference(T)
		perms=[conversions.cycles_to_signed_permutation(n, str(g)) for g in T]
		perms = sorted(perms, key = lambda perm : str(perm).replace('-', 'Z')) # always sort these. Why not!?
		if classes_as == 'counts':
			cards.append((perms[0],len(perms)))
		elif classes_as == 'lists':
			cards.append(perms)
		elif classes_as == 'dict':
			cards[perms[0]] = perms
	if sort_classes and classes_as != 'dict':
		cards = sorted(cards, key = lambda card : str(card[0]).replace('-', 'Z'))
	return cards