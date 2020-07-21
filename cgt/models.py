## Models 

from enum import Enum, auto
from sage.all_cmdline import *
from cgt import conversions, helper_functions

class MODEL(Enum):
	TWO_REGION_INVERSIONS = auto()

def model(G, n, signed=True, model_type = MODEL.TWO_REGION_INVERSIONS):
	print("WARNING: only working for sage 'SymmetricGroup' G (signed or unsigned)")
	if model_type == MODEL.TWO_REGION_INVERSIONS:
		return __all_inversions_model(G, n, signed, num_regions = 2)
	
def __all_inversions_model(G, n, signed, num_regions):
	if num_regions == 2:
		if signed:
			return [G(f'({i},{-1*(i+1)})({-1*i},{i+1})') for i in range(1,n)] + [G(f'({n},-1)(-{n},1)')]
		else:
			return [G(f'({i},{i+1})') for i in range(1,n)] + [G(f'(1,{n})')]