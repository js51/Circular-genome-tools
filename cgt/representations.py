from sage.all_cmdline import *
from cgt import hyperoctahedral_groups
import numpy as np

def regular_representation(G, g): 
	"""regular representation """
	C = QQbar # Complex field
	CG =  GroupAlgebra(G,C)
	return matrix(QQbar, [(CG(g)*G(h)).to_vector(QQbar) for h in G]).transpose()
	
def irreducible_representations(n, signed=False):
	"""Return a complete system of pairwise irreducible representations of either S_n or H_n"""
	if not signed:
		representatives = SymmetricGroupRepresentations(n).list()
		return representatives
	elif signed:
		G = hyperoctahedral_groups.HyperoctahedralGroup(n)
		characters = gap.Irr(G)
		representations = []
		for character in characters:
			irrep = gap.IrreducibleAffordingRepresentation(character)
			def representation(sigma, as_gap_matrix=False, _irrep=irrep):	
				image = gap.Image(_irrep, sigma)
				return image if as_gap_matrix else matrix(UniversalCyclotomicField(), image)
			representations.append(representation)
		return representations