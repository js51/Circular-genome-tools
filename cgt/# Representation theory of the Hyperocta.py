#%% Representation theory of the Hyperoctahedral Group
from cgt import *
from sage.all import (
    SymmetricGroup as SymmetricGroup,
    cartesian_product,
    GroupSemidirectProduct,
    Permutation,
    CyclicPermutationGroup,
    SignedPermutations,
)

#%% Definitions


def action(s, c):
    try:
        _ = len(s)  # works if it is a "tuple"
    except:
        pass
    else:
        Sn = SymmetricGroup(n)
        s = Sn(s[0]) * Sn(s[1])
    return C2n((c[s(k) - 1] for k in range(1, n + 1)))


class as_set(Enum):
    SemidirectProduct = auto()
    WreathProduct = auto()
    SignedPermutations = auto()
    SubgroupOfSymmetric = auto()


def semidirect_product(H, G, action):
    HG = GroupSemidirectProduct(H, G, act_to_right=False, twist=action)

    def random_element():
        C, Sn = HG.cartesian_factors()
        return HG((C.random_element(), Sn.random_element()))

    def order():
        return H.order() * G.order()

    def elements():
        return [HG((c, s)) for c in H for s in G]

    HG.random_element = random_element
    HG.action = action
    HG.order = order
    HG.elements = elements
    return HG


def signed_permutations(n):
    S = tuple(range(1, n + 1)) + tuple(range(-n, 0))
    s_pmn = SymmetricGroup(S)
    if n == 1:
        return s_pmn
    h_gens = [
        s_pmn([(1, 2), (-1, -2)]),
        s_pmn(tuple(range(1, n + 1)) + tuple(range(-1, -(n + 1), -1))),
    ]
    return s_pmn.subgroup(h_gens)


def hyperoctahedral_group(n, compute_as=as_set.SemidirectProduct):
    """
    Returns the hyperoctahedral group of order 2^n * n!
    """
    C2n = cartesian_product([SymmetricGroup(2) for _ in range(n)])
    Sn = SymmetricGroup(n)
    if compute_as is as_set.SemidirectProduct:
        return semidirect_product(C2n, Sn, action)
    elif compute_as is as_set.SignedPermutations:
        return signed_permutations(n)
    else:
        raise NotImplementedError


def cycles_to_one_row(n, cycles):
    row = list(cycles.dict().values())[0:n]
    row = SignedPermutations(n)(row)
    return row


def Phi(g, one_row=True):
    """Maps semidirect product elements to signed permutations."""
    c, s = g
    n = len(c)
    Hn = hyperoctahedral_group(n, compute_as=as_set.SignedPermutations)
    positive_row = [c[k - 1].sign() * s(k) for k in range(1, n + 1)]
    other_half = list(-int(positive_row[n - i]) for i in range(1, n + 1))
    result = Hn(positive_row + other_half)
    if one_row:
        return cycles_to_one_row(n, result)
    else:
        return result


def Phi_inv(sig):
    n = len(list(sig.tuple())) / 2
    C2 = SymmetricGroup(2)
    Sn = SymmetricGroup(n)
    c = tuple(C2(() if sig(k) > 1 else (1, 2)) for k in range(1, n + 1))
    s = Sn(Permutation([abs(sig(k)) for k in range(1, n + 1)]))
    Hn = hyperoctahedral_group(n, compute_as=as_set.SemidirectProduct)
    return Hn((c, s))


def prod_to_one_row_string(n, g):
    c, s = g
    return str(c) + " * " + str(cycles_to_one_row(n, s))


def little_subgroup(Sn, l, m, as_direct_product=False):
    n = Sn.degree()
    if l + m != n or l < 0 or m < 0 or n < 0:
        raise ValueError("l,m must be non-negative and l+m must equal n")
    if l == 0 or m == 0 and not as_direct_product:
        return Sn
    Sl = SymmetricGroup(list(range(1, l + 1)))
    if as_direct_product:
        Sm = SymmetricGroup(list(range(1, n - l + 1)))
        return cartesian_product((Sl, Sm))
    Sm = SymmetricGroup(list(range(l + 1, n + 1)))
    return Sn.subgroup(Sl.gens() + Sm.gens())


def little_subgroup_pairs(n):
    return [(l, n - l) for l in range(n + 1)]


def little_subgroups(n):
    Sn = SymmetricGroup(n)
    return [little_subgroup(Sn, l, m) for l, m in little_subgroup_pairs(n)]


def orbit_representative(l, m):
    C2 = SymmetricGroup(2)
    C2n = cartesian_product([C2 for _ in range(l + m)])
    return C2n(tuple(() for _ in range(l)) + tuple((1, 2) for _ in range(m)))


def gap_wreath_to_sage_semidirect(n, gap_elt):
    elt_string = str(gap_elt)[21:-1].split(",")
    c = tuple((1, 2) if x == "f1" else () for x in elt_string[0:n])
    if elt_string[n:][0] == "()":
        s = ()
    else:
        s = ",".join(elt_string[n:])
    return hyperoctahedral_group(n)((c, s))


def stabiliser_of(object, group, action):
    stabiliser = set()
    for elt in group:
        if action(elt, object) == object:
            stabiliser.add(elt)
    return stabiliser


def right_transversal(group, subgroup):
    transversal = [subgroup.one()]
    num_cosets = group.order() / subgroup.order()
    while len(transversal) < num_cosets:
        elt = group.random_element()
        if all(elt.inverse() * t not in subgroup for t in transversal):
            transversal.append(elt)
    return transversal


def character_for_tuple(c):
    C2 = SymmetricGroup(2)
    return tuple(
        (lambda x: 1) if entry == C2(()) else (lambda x: x.sign()) for entry in c
    )


def symmetric_group_representations(n):
    partitions_of_n = Partitions(n)
    return {p: SymmetricGroupRepresentation(p) for p in partitions_of_n}


def subgroup_elt_to_direct_product(elt, l, m):
    """
    Takes elt as an element of the little subgroup of Sn of order l!m! and returns an element of the direct product of S_l and S_m.
    """
    Sn = SymmetricGroup(l + m)
    Sl = SymmetricGroup(l)
    Sm = SymmetricGroup(m)
    SlxSm = cartesian_product((Sl, Sm))
    cycles = list(elt)
    left = Sl.one()
    right = Sm.one()
    for cycle in cycles:
        if cycle in Sl:
            left *= Sl(cycle)
        else:
            string = str(cycle)
            for new, old in enumerate(range(l + 1, n + 1)):
                string = string.replace(str(old), str(new + 1))
            right *= Sm(string)
    return SlxSm((left, right))


#%% Construct Hyperoctahedral Group as a wreath product
n = 8
C2, Sn = SymmetricGroup(2), SymmetricGroup(n)
C2n = cartesian_product([C2 for _ in range(n)])
Hn = hyperoctahedral_group(n)
action = Hn.action

# %%
Hn = hyperoctahedral_group(n, compute_as=as_set.SemidirectProduct)
pairs = little_subgroup_pairs(n)

partition_pair = (Partitions(6)[4], Partitions(2)[1])


# %%


def rep_factory(partition_pair):
    rep_left = SymmetricGroupRepresentation(partition_pair[0])
    rep_right = SymmetricGroupRepresentation(partition_pair[1])
    l, m = sum([0] + list(partition_pair[0])), sum([0] + list(partition_pair[1]))
    n = l + m
    left_group = SymmetricGroup(l)
    right_group = SymmetricGroup(m)
    c = orbit_representative(l, m)
    x = character_for_tuple(c)

    def rep(g):
        try:
            _ = len(g[1])  # works if it is a "tuple"
        except:  # g[1] is not a "tuple" of permutations
            g = (g[0], subgroup_elt_to_direct_product(g[1], l, m))
        prod = rep_left(g[1][0]).tensor_product(rep_right(g[1][1]))
        character = (x[i](ci) for i, ci in enumerate(g[0]))
        coeff = 1
        for char in character:
            coeff *= char
        return coeff * prod

    return rep


def induced_rep_factory(partition_pair):
    rep = rep_factory(partition_pair)
    l, m = sum([0] + list(partition_pair[0])), sum([0] + list(partition_pair[1]))
    subgroup = little_subgroup(Sn, l, m)
    C2n_wr_SlxSm = semidirect_product(C2n, subgroup, action)
    tv = right_transversal(Hn, C2n_wr_SlxSm)

    def irrep(g):
        id_rep = rep(C2n_wr_SlxSm.one())
        rep_dimension = (id_rep.nrows(), id_rep.ncols())
        Y = [[None for _ in range(len(tv))] for _ in range(len(tv))]
        for i in range(len(tv)):
            for j in range(len(tv)):
                element = tv[i] * g * tv[j].inverse()
                if element in C2n_wr_SlxSm:
                    Y[i][j] = rep(element)
                else: # A matrix of zeros
                    Y[i][j] = matrix(*rep_dimension)
        return block_matrix(Y, subdivide=False)

    return irrep


# %%
