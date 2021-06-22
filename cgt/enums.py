from enum import Enum, auto

class MODEL(Enum):
    one_region_inversions  = auto()
    two_region_inversions  = auto()
    all_inversions		   = auto()
    all_inversions_		   = auto()
    one_region_swaps       = auto()

class TYPE(Enum):
    reg_to_signed_pos = auto()
    pos_to_signed_reg = auto()
    signed_reg_to_pos = auto()
    signed_pos_to_reg = auto()

class FORMAT(Enum):
    formal_sum     = auto()
    equiv_classes  = auto()
    dictionary     = auto()
    only_reps      = auto()

class SET(Enum):
	signed_cycles 	= auto()
	unsigned_cycles = auto()
	one_row			= auto()
	wreath			= auto()
	wreath_s2		= auto()
	
class SYMMETRY(Enum):
    circular = dihedral = D_n              = auto()
    linear = S_2 = C_2 = flip = reflection = auto()