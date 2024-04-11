"""
"""

from enum import Enum, auto

class MODEL(Enum):
    # Inversions
    all_inversions         = "all inversions"
    one_region_inversions  = "inversions of a single region"
    two_region_inversions  = "inversions of two adjacent regions"
    # Transpositions
    all_transpositions  = "all transpositions"
    ## "Swap" interpretation
    two_region_transpositions    = "Two region adjacent transpositions, with or without inversion"
    two_region_transpositions_without_inversions = "Two region adjacent transpositions, without inversion"
    two_region_revrevs = "Take two adjacent regions and invert both of them"
    ## "Move" interpretation
    one_region_moves  = "Move a single region, with or without inversion"
    one_region_moves_without_inversions  = "Move a single region, without inversion"
    # Other transpositions
    three_region_transpositions = "transpositions involving a segment of three regions, with or without inversion"
    # Deprecated options
    two_region_adjacent_transpositions = "transpositions of two adjacent regions"

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

class IRREP_TYPE(Enum):
    specht = speck = auto()
    orthogonal     = auto()
    seminormal     = auto()

class DISTANCE(Enum):
    min = minimum = minimum_distance = auto()
    min_weighted = minimum_weighted = auto()
    MFPT = mean_first_passage_time = auto()
    MLE = maximum_likelihood_distance = maximum_likelihood_estimate = auto()
    discrete_MFPT = DMFPT = discrete_mean_first_passage_time = auto()

class DATA(Enum):
    eig_data = 'eigen_data'
    eigval_lists = 'eigenvalue_lists'
    eigval_sets = 'eigenvalue_sets'
    eigvec_lists = 'eigenvector_lists'
    eigvec_mat_inv = "eigenvector_matrix_inverses"
    projections = "projections"
    eigval_sets_old = 'eig_lists'
