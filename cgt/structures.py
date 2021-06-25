"""
A place for functions which return important but general algebraic structures requiring only a few parameters
"""

from sage.all_cmdline import *

def set_plus_minus(n):
    return tuple(range(1,n+1)) + tuple(range(-n, 0))

def HyperoctahedralGroup(n):
    """Return the hyperoctahedral group of order (2^n)n! as a subgroup of SymmetricGroup(2n) with two generators."""
    s_pmn = SymmetricGroup(set_plus_minus(n))
    if n==1: return s_pmn
    h_gens = [
        s_pmn([(1,2),(-1,-2)]),
        s_pmn(tuple(range(1,n+1)) + tuple(range(-1, -(n+1), -1)))
    ]
    return s_pmn.subgroup(h_gens)