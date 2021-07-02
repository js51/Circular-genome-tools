"""
"""

from enum import Enum, auto

class MODEL(Enum):
    one_region_inversions  = "inversions of a single region"
    two_region_inversions  = "inversions of two adjacent regions"
    all_inversions		   = "all inversions"
    one_region_swaps       = "adjacent transpositions"

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

class CLASSES(Enum):
    conjugacy_classes = auto()
    cosets            = auto()
    double_cosets     = auto()

class DISPLAY(Enum):
    arrows  = auto()
    one_row = auto()
    cycles  = auto()

class ALGEBRA(Enum):
    group        = auto()
    genome       = auto()
    genome_class = auto()