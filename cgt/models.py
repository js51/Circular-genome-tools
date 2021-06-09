## Models 

from enum import Enum, auto
from sage.all_cmdline import *
from cgt import conversions, helper_functions

class MODEL(Enum):
	ONE_AND_TWO_REGION_INVERSIONS 	= auto()
	TWO_REGION_INVERSIONS 			= auto()
	ALL_INVERSIONS					= auto()

# TODO: add the cuts function here, enable creation of more models
# a model class which consists of a dictionary (?) of model elements mapping to their probability
# write two checks for model properties, and the ability to return a model but closed under certain
# algebraic properties

def css(cuts_set):
		string = ''.join(sorted([str(cut[0]) for cut in cuts_set]))
		if len(string) != len(cuts_set): print("something went very wrong here...")
		return string

def cuts(sigma, n):
	pmN = list(range(-n, 0)) + list(range(1, n+1))	# The set {+-1, ..., +-n}
	G = SymmetricGroup(pmN)
	c_string = f'({",".join(str(i) for i in list(range(1, n+1)))})({",".join(str(i) for i in list(range(-n, 0)))})'
	c = G(c_string)
	sigma = list(sigma)
	return set([ 
		tuple(([i+1,c(i+1)])) for i in range(len(sigma)) 
		if (sigma[c(i + 1)-1] != c(sigma[i]))
	])

def model(G, n, signed=True, model_type = MODEL.ONE_AND_TWO_REGION_INVERSIONS):
	# print("WARNING: only working for subgroups of sage 'SymmetricGroup' G (signed or unsigned)")
	if model_type == MODEL.TWO_REGION_INVERSIONS:
		return set(sorted(__all_inversions_model(G, n, signed, num_regions = 2)))
	elif model_type == MODEL.ONE_AND_TWO_REGION_INVERSIONS:
		return set(sorted(__all_inversions_model(G, n, signed, num_regions = 1) + __all_inversions_model(G, n, signed, num_regions = 2)))
	elif model_type == MODEL.ALL_INVERSIONS:
		return set(sorted(__all_inversions_model(G, n, signed)))

def __all_inversions_model(G, n, signed, num_regions='any'):
	if num_regions == 1:
		if signed:
			return [G(f'({i},-{i})') for i in range(1, n+1)]
		else:
			return []
	elif num_regions == 2:
		if signed:
			return [G(f'({i},{-1*(i+1)})({-1*i},{i+1})') for i in range(1,n)] + [G(f'({n},-1)(-{n},1)')]
		else:
			return [G(f'({i},{i+1})') for i in range(1,n)] + [G(f'(1,{n})')]
	elif num_regions == 'any':
		if signed: 
			M = []
			for permutation in G:
				if set(permutation.cycle_type()) == {1,2} and cuts(conversions.cycles_to_signed_permutation(n, permutation), n) == 2:
					M.append(permutation)
			return M
		else:
			return 'not yet implemented'
		
	