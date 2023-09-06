"""
"""

from copy import deepcopy
from sage.all_cmdline import SymmetricGroup, floor, ceil
import warnings
from . import structures
from .enums import *


def inversion(framework, about_position, length):
    """
    Return an instance of the inversion that inverts a segment in a genome.

    Args:
        framework (PositionParadigmFramework): The framework to create the inversion in
        about_position (int): The position to invert about
        length (int): The length of the segment to invert

    Returns:
        Permutation: An instance representing the inversion that inverts the segment
    """
    if not framework.oriented or framework.symmetry != SYMMETRY.circular:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")
    n = framework.n
    r = about_position
    string_rep = ""
    if length % 2 == 0:
        k = int(length / 2)
        for i in range(1, k + 1):
            string_rep += (
                f"({(r-(i-1)) % n}, -{(r+i) % n})(-{(r-(i-1)) % n}, {(r+i) % n})"
            )
    else:
        k = int((length - 1) / 2)
        string_rep += f"({about_position},{-about_position})"
        for i in range(1, k + 1):
            string_rep += f"({(r-i) % n},-{(r+i) % n})(-{(r-i) % n},{(r+i) % n})"
    string_rep = string_rep.replace("0", str(n))
    return framework.cycles(string_rep)


signed_inversion = inversion  # old name


def all_inversion_instances(framework, length=None):
    """
    Return all instances that represent an inversion of the given length

    Args:
        framework (PositionParadigmFramework): The framework to create the inversion in
        length (int, optional): The length of the inversion. Defaults to None (all inversions).

    Returns:
        set: The set of all instances that represent an inversion of the given length
    """
    if not framework.oriented or framework.symmetry != SYMMETRY.circular:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")

    n = framework.n
    rearrangements = []

    for r in range(1, n + 1):
        for l in (length,) if length is not None else range(1, n + 1):
            rearrangements.append(inversion(framework, r, l))

    return rearrangements


def transposition(framework, sec_1, sec_2, inv_1=False, inv_2=False, revrev=False):
    """
    Return an instance of the transposition that swaps two segments in a genome.

    Args:
        framework (PositionParadigmFramework): The framework to create the transposition in
        sec_1 (tuple): The first section to transpose (tuple of start (included) and end (excluded) positions)
        sec_2 (tuple): The second section to transpose (tuple of start (included) and end (excluded) positions)
        inv_1 (bool, optional): Whether to also invert the first section. Defaults to False.
        inv_2 (bool, optional): Whether to also invert the second section. Defaults to False.

    Returns:
        Permutation: An instance representing the transposition that swaps the two segments
    """
    if not framework.oriented or framework.symmetry != SYMMETRY.circular:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")
    if all((inv_1, inv_2)) and not revrev:
        raise ValueError(
            "for circular genomes, inverting both segments is not allowed. Explicitly allow revrevs by setting allow_revrev=True"
        )
    if revrev:
        if not all((inv_1, inv_2)):
            raise ValueError("revrevs are only allowed if both segments are inverted")
        else:
            inv_1 = inv_2 = False
    n = framework.n

    if sec_1[1] != sec_2[0]:
        raise ValueError("sections must be adjacent")
    full_inversion = (
        signed_inversion(framework, segment_midpoint(n, sec_1[0], sec_2[1]), segment_length(n, sec_1[0], sec_2[1]))
        if not revrev
        else framework.cycles("()")
    )
    left_inversion = (
        signed_inversion(framework, segment_midpoint(n, *sec_1), segment_length(n, *sec_1))
        if not inv_1
        else framework.cycles("()")
    )
    right_inversion = (
        signed_inversion(framework, segment_midpoint(n, *sec_2), segment_length(n, *sec_2))
        if not inv_2
        else framework.cycles("()")
    )
    return left_inversion * right_inversion * full_inversion


def all_transposition_instances(
        framework, scope_limit = None, single_segment_limit = None, with_inversion = True, include_revrevs = False, only_revrevs=False, canonical_reps_only = False):
    """
    Return all instances that represent a transposition of the given length

    Args:
        framework (PositionParadigmFramework): The framework to create the transposition in
    Returns:
        set: The set of all instances that represent a transposition
    """
    if not framework.oriented or framework.symmetry != SYMMETRY.circular:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")

    n = framework.n
    rearrangements = set()

    if scope_limit is None:
        scope_limit = (n-1)

    for full_length in range(2, scope_limit + 1): # Scope of the rearrangement
        for start in range(1, n+1): # Start of first segment

            if single_segment_limit is not None:
                first_segment_lengths = list(range(1, single_segment_limit + 1)) + list(range(full_length - single_segment_limit, full_length))
            else:
                first_segment_lengths = range(1, full_length)

            for l1 in first_segment_lengths: # Length of first segment
                l2 = full_length - l1 # length of second segment
                middle = (start + l1 - 1) % n + 1
                end = (middle + l2 - 1) % n + 1
                s, m, e = start, middle, end
                possible_transpositions = { transposition(framework, (s, m), (m, e)) } if not only_revrevs else set()
                if with_inversion and not only_revrevs:
                    possible_transpositions |= {
                        transposition(framework, (s, m), (m, e), inv_1=True),
                        transposition(framework, (s, m), (m, e), inv_2=True)
                    }
                if include_revrevs or only_revrevs:
                    possible_transpositions |= {
                        transposition(framework, (s, m), (m, e), inv_1=True, inv_2=True, revrev=True)
                    }
                rearrangements = rearrangements.union(possible_transpositions)

    if canonical_reps_only:
        rearrangements = __representatives(framework, rearrangements, prioritise_string_length=True)

    return rearrangements

def c_perm(n):
    """Return the permutation (1,...,n)(-n,...,1))"""
    pmN = structures.set_plus_minus(n)
    G = SymmetricGroup(pmN)
    c_string = f'({",".join(str(i) for i in list(range(1, n+1)))})({",".join(str(i) for i in list(range(-n, 0)))})'
    c = G(c_string)
    return c


def cut_positions(framework, sigma):
    """
    The set of cuts required to perform a rearrangement represented by sigma

    Args:
        framework (PositionParadigmFramework): The framework to create the transposition in
        sigma (Permutation): The permutation representing the rearrangement

    Returns:
        set: The set of cuts required to perform the rearrangement represented by sigma
    """
    sigma = deepcopy(framework.one_row(sigma))
    n = framework.n
    c = c_perm(n)
    sigma = list(sigma)
    return set(
        [
            tuple(([i + 1, c(i + 1)]))
            for i in range(len(sigma))
            if (sigma[c(i + 1) - 1] != c(sigma[i]))
        ]
    )


def __all_canonical_inversions(framework, num_regions=None):
    if not framework.oriented:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")
    if framework.symmetry not in {SYMMETRY.circular, SYMMETRY.linear}:
        raise NotImplementedError(
            f"not implemented for framework with {str(framework.symmetry)} symmetry"
        )
    n = framework.n
    G = framework.genome_group()
    if num_regions == 1:
        return {G(f"({i},-{i})") for i in range(1, n + 1)} if framework.oriented else {}
    elif num_regions == 2:
        if framework.oriented:
            perms = {G(f"({i},{-1*(i+1)})({-1*i},{i+1})") for i in range(1, n)}
            if framework.symmetry is not SYMMETRY.linear:
                perms.union({G(f"({n},-1)(-{n},1)")})
        else:
            perms = {G(f"({i},{i+1})") for i in range(1, n)}
            if framework.symmetry is not SYMMETRY.linear:
                perms.union({G(f"(1,{n})")})
        return perms
    elif num_regions == None:  # Return all inversions up to length floor(n/2)
        up_to_length = floor(framework.n / 2)
        if framework.oriented:
            perms = set()
            for permutation in G:
                cycle_type = list(permutation.cycle_type())
                if (
                    set(cycle_type) == {1, 2}
                    and list(cycle_type).count(2) <= up_to_length
                    and len(cut_positions(framework, permutation)) == 2
                ):
                    perms.add(permutation)
            if framework.symmetry is SYMMETRY.linear:
                perms = {
                    perm
                    for perm in perms
                    if not ("1" in str(perm) and str(n) in str(perm))
                }
            return perms
        else:
            raise NotImplementedError(f"not yet implemented for {framework}")
    else:
        raise NotImplementedError(
            f"inversions of length {num_regions} not yet implemented"
        )


def __two_region_adjacent_transposition_reps(framework):
    if not framework.oriented or framework.symmetry != SYMMETRY.circular:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")
    G = framework.genome_group()
    return {G("(-2,-1)(1,2)"), G("(-2,1,2,-1)"), G("(-2,-1,2,1)"), G("(-2,2)(-1,1)")}


def double_coset(framework, perm):
    """Return the double coset of perm in the symmetry group of framework"""
    Z = framework.symmetry_group()
    return {d1 * perm * d2 for d1 in Z for d2 in Z}


def conjugacy_class(framework, perm):
    """Return the conjugacy class of perm in the symmetry group of framework"""
    Z = framework.symmetry_group()
    return {d.inverse() * perm * d for d in Z}


def single_coset(framework, perm):
    """Return the single coset of perm in the symmetry group of framework"""
    Z = framework.symmetry_group()
    return {perm * d for d in Z}

def segment_length(n, start, end):
    """Return the length of the segment from start to end"""
    return (end - start - 1) % n + 1

def segment_midpoint(n, start, end):
    """Return the midpoint of the segment from start to end"""
    return ((start + floor((segment_length(n,start, end) - 1) / 2) - 1) % n) + 1

def __representatives(framework, set_of_permutations, classes=CLASSES.double_cosets, prioritise_string_length=False):
    if not framework.oriented:
        raise NotImplementedError(f"not yet implemented for {str(framework)}")
    Z = framework.symmetry_group()
    perms = {x for x in set_of_permutations}
    cosets = []
    while len(perms):
        perm = perms.pop()
        if classes is CLASSES.double_cosets:
            coset = double_coset(framework, perm)
        elif classes is CLASSES.conjugacy_classes:
            coset = conjugacy_class(framework, perm)
        elif classes is CLASSES.cosets:
            coset = single_coset(framework, perm)
        for element in coset:
            perms.discard(element)
        cosets.append(coset)
    reps = set()

    if prioritise_string_length:
        sort_key = lambda x: (len(str(x).replace('-','')), str(x))

    else:
        sort_key = lambda x: (
            sum(set(x.cycle_type()) - {1}),
            list(x.cycle_type()).count(2),
            str(x),
        )

    for coset in cosets:
        reps.add(
            sorted(
                list(coset),
                key=sort_key,
            )[0]
        )
    return reps

representatives = __representatives


def all_inversions_representatives(framework, num_regions=None):
    """
    Return the representatives of all inversions in the symmetry group of framework

    Args:
        framework (PositionParadigmFramework): The framework to create the transposition in
        num_regions (int, optional): The number of regions to invert. Defaults to None (all inversions).

    Returns:
        set: The set of representatives of all inversions in the symmetry group of framework
    """
    return __representatives(
        framework,
        __all_canonical_inversions(framework, num_regions=num_regions),
        classes=CLASSES.double_cosets,
    )


def all_adjacent_transpositions_representatives(framework, num_regions=None):
    if num_regions == 2:
        return __two_region_adjacent_transposition_reps(framework)
    else:
        raise NotImplementedError(f"model not yet implemented")


def permutation_with_cuts(framework, cuts, perm=None, start=None):
    """
    A generator for all permutations with the specified set of cuts.

    Args:
        framework (PositionParadigmFramework): The framework to create the transposition in
        cuts (set): The set of cuts to use

    Yields:
        Permutation: A permutation with the specified set of cuts
    """
    n = framework.n
    c = c_perm(n)
    if perm is None:
        if not (framework.oriented and framework.symmetry == SYMMETRY.circular):
            raise NotImplementedError(f"not yet implemented for {str(framework)}")
        perm = {}
        perm[1] = 1  # making the canonical instance
        results = permutation_with_cuts(cuts, perm, 2)
        for result in results:
            if result:
                yield list(result.values())
    else:  # recurse
        possibilities = set(range(1, n + 1)) | set(range(-n, 0))
        for val in perm.values():
            possibilities.remove(val)
            possibilities.remove(-val)
        if start - 1 in cuts:
            if c(perm[start - 1]) in possibilities:
                possibilities.remove(c(perm[start - 1]))
        else:
            possibilities = possibilities & {c(perm[start - 1])}
        if start == n:  # also check if *start* is a cut
            if start in cuts:
                # add exclusions
                if start in possibilities:
                    possibilities.remove(start)  # p[1]=1 so we can't have p[n]=n
            else:
                possibilities = possibilities & {n}  # we must choose n
            for possibility in possibilities:
                bperm = perm.copy()
                bperm[start] = possibility
                yield bperm
        else:
            for possibility in possibilities:
                bperm = perm.copy()
                bperm[start] = possibility
                for result in permutation_with_cuts(cuts, bperm, start + 1):
                    yield result
