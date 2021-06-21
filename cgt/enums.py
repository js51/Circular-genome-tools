from enum import Enum, auto

class SYMMETRY(Enum):
    circular = dihedral = D_n = auto()
    linear = S_2 = C_2 = flip = reflection = auto()

class MODEL(Enum):
	ONE_AND_TWO_REGION_INVERSIONS 			= auto()
	TWO_REGION_INVERSIONS 					= auto()
	ALL_INVERSIONS							= auto()
	ALL_INVERSIONS_HALF_CIRCLE				= auto()
	SHORT_INVERSIONS_AND_ONE_REGION_SWAPS   = auto()

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
	SIGNED_CYCLES 	= auto()
	UNSIGNED_CYCLES = auto()
	ONE_ROW			= auto()
	WREATH			= auto()
	WREATH_S2		= auto()
	