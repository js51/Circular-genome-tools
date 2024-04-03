"""
This file contains the class for the hyperoctahedral group.
"""

# %%
from cgt import structures
import sage
import numpy as np
from sage.all_cmdline import (
    SymmetricGroup,
    cartesian_product,
    GroupSemidirectProduct,
    SymmetricGroupRepresentation,
    SignedPermutations,
    matrix,
    block_matrix,
    Partitions,
    Permutation,
)


def HyperoctahedralGroupRepresentations(n):
    """
    This function returns the irreducible representations of the hyperoctahedral group of order 2^n * n!, indexed by pairs of partitions of l and m, where l + m = n.

    Args:
        n (int): The order of the hyperoctahedral group.

    Returns:
        dict: A dictionary of the irreducible representations of the hyperoctahedral group of order 2^n * n!
    """
    irreps = {}
    Hn = hyperoctahedral_group(n)
    for little_subgroup_pair in Hn.little_subgroup_pairs():
        for partition_left in Partitions(little_subgroup_pair[0]):
            for partition_right in Partitions(little_subgroup_pair[1]):
                partition_pair = (partition_left, partition_right)
                irrep = Hn.irrep(partition_pair)

                def wrapped_irrep(elt, _irrep=irrep):
                    return _irrep(Hn.Phi(elt))

                irreps[partition_pair] = wrapped_irrep
    return irreps


class hyperoctahedral_group:
    """
    This class represents the hyperoctahedral group of order 2^n * n!. It stores the group as a number of different isomorphic structures, as well as a number of other things that are helpful to keep track of in the context of genome rearrangements. It implements functions to compute the irreducible representations of the hyperoctahedral group.
    """

    def __init__(self, n, default_element_type="signed_permutation"):
        """
        Args:
            n (int): The order of the hyperoctahedral group.
            default_element_type (str, optional): The type of element to use when generating random elements. Defaults to "signed_permutation".
        """
        self._n = n
        self._signed_perms = structures.HyperoctahedralGroup(self._n)
        self._sage_signed_perms = SignedPermutations(self._n)
        self._Sn = SymmetricGroup(self._n)
        self._C2 = SymmetricGroup(2)
        self._N = self._signed_perms.subgroup(
            [self._signed_perms((i, -i)) for i in range(1, self._n + 1)]
        )
        self._Sn_in_Hn = self._signed_perms.subgroup(
            [
                self._signed_perms(
                    f'({",".join(str(x) for x in range(1,n+1))})({",".join(str(-x) for x in range(1,n+1))})'
                ),
                self._signed_perms(f"(1,2)(-1,-2)"),
            ]
        )
        self._Sn_wr_N = GroupSemidirectProduct(self._Sn_in_Hn, self._N)
        self._C2n = cartesian_product(tuple(SymmetricGroup(2) for _ in range(self._n)))
        self._semidirect_product = self.semidirect_product_with_hyperoctahedral_twist(
            self._Sn, self._C2n
        )

    def hyperoctahedral_twist(self, s, c):
        """
        This function implements the twist of the hyperoctahedral group, which is the action of S_n on (C_2)^n by permuting the factors.

        Args:
            s (Permutation): A permutation in S_n.
            c (tuple): A tuple of elements of C_2.

        Returns:
            an element of (C_2)^n: The tuple of elements of C_2 obtained by applying s to the tuple c.
        """
        try:
            _ = len(s)  # works if it is a "tuple"
        except:
            pass  # We have the usual hyperoctahedral group
        else:  # We have a wreath product of a little subgroup (as a direct product)
            #s = self._Sn(s[0]) * self._Sn(s[1])
            raise ValueError("This is not implemented yet")
        return self._C2n((c[s(k) - 1] for k in range(1, self._n + 1)))

    def orbit_representative(self, l, m):
        return self._C2n(tuple(() for _ in range(l)) + tuple((1, 2) for _ in range(m)))

    def little_subgroup_pairs(self):
        return tuple((l, self._n - l) for l in range(self._n + 1))

    def little_subgroup(self, l, m, as_direct_product=False):
        if l + m != self._n or l < 0 or m < 0 or self._n < 0:
            raise ValueError("l,m must be non-negative and l+m must equal n")
        if l == 0 or m == 0 and not as_direct_product:
            return self._Sn
        Sl = SymmetricGroup(list(range(1, l + 1)))
        if as_direct_product:
            Sm = SymmetricGroup(list(range(1, self._n - l + 1)))
            return cartesian_product((Sl, Sm))
        Sm = SymmetricGroup(list(range(l + 1, self._n + 1)))
        return self._Sn.subgroup(Sl.gens() + Sm.gens())

    def little_subgroup_pairs(self):
        return [(l, self._n - l) for l in range(self._n + 1)]

    def little_subgroups(self):
        return [
            self.little_subgroup(l, m) for l, m in self.little_subgroup_pairs(self._n)
        ]

    def orbit_representative(self, l, m):
        return self._C2n(tuple(() for _ in range(l)) + tuple((1, 2) for _ in range(m)))

    def character_for_tuple(self, c):
        return tuple(
            (lambda x: 1) if entry == self._C2(()) else (lambda x: x.sign())
            for entry in c
        )

    def representation_little_subgroup(self, partition_pair):
        rep_left = SymmetricGroupRepresentation(partition_pair[0])
        rep_right = SymmetricGroupRepresentation(partition_pair[1])
        l, m = sum([0] + list(partition_pair[0])), sum([0] + list(partition_pair[1]))
        n = l + m
        left_group = SymmetricGroup(l)
        right_group = SymmetricGroup(m)
        c = self.orbit_representative(l, m)
        x = self.character_for_tuple(c)

        def rep(g):
            try:
                _ = len(g[0])  # works if it is a "tuple"
            except:  # g[0] is not a "tuple" of permutations
                g = (self.element_in_little_subgroup(g[0], l, m), g[1])
            prod = rep_left(g[0][0]).tensor_product(rep_right(g[0][1]))
            character = (x[i](ci) for i, ci in enumerate(g[1]))
            coeff = 1
            for char in character:
                coeff *= char
            return coeff * prod

        return rep


    def _elt_in_subgroup(self, elt, l):
        return all(
            (
                all(int(x) <= l for x in str(cycle).strip("()").split(","))
                or all(int(x) > l for x in str(cycle).strip("()").split(","))
            )
            for cycle in elt
        )

    def irrep(self, partition_pair):
        """
        Returns the irreducible representation of the hyperoctahedral group of order 2^n * n! indexed by the pair of partitions of l and m where l + m = n. The representation is given as a function that takes an element of the hyperoctahedral group as input and returns the corresponding matrix.

        Args:
            partition_pair (tuple): A pair of partitions of l and m where l + m = n.

        Returns:
            function (group element -> real matrix): A function that takes an element of the hyperoctahedral group as input and returns the corresponding matrix.
        """
        rep = self.representation_little_subgroup(partition_pair)
        l, m = sum([0] + list(partition_pair[0])), sum([0] + list(partition_pair[1]))
        subgroup = self.little_subgroup(l, m)
        C2n_wr_SlxSm = self.semidirect_product_with_hyperoctahedral_twist(
            subgroup, self._C2n
        )
        tv = self.right_transversal_of_little_subgroup(l)
        id_rep = rep(C2n_wr_SlxSm.one())
        rep_dimension = (id_rep.nrows(), id_rep.ncols())

        def _irrep(g):
            Y = [[None for _ in range(len(tv))] for _ in range(len(tv))]
            for row in range(len(tv)):
                for col in range(len(tv)):
                    element = tv[col] * g * tv[row].inverse()
                    if self._elt_in_subgroup(element[0], l):
                        Y[row][col] = np.array(rep(element), dtype=object)
                    else:  # A matrix of zeros
                        Y[row][col] = np.zeros(rep_dimension, dtype=sage.rings.rational.Rational)

            return matrix(np.block(Y))

        return _irrep

    def element_in_little_subgroup(self, elt, l, m):
        """
        Takes elt as an element of the little subgroup of Sn of order l!m! and returns an element of the direct product of S_l and S_m.
        """
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
                for new, old in enumerate(range(l + 1, self._n + 1)):
                    string = string.replace(str(old), str(new + 1))
                right *= Sm(string)
        return SlxSm((left, right))

    def semidirect_product_with_hyperoctahedral_twist(self, G, N):
        """
        Compute a semi-direct product of H and G with the hyperoctahedral twist (G acts on H^n by permuting factors). The group returned is a group of pairs (c, s) where c is an element of H^n and s is an element of G. The group operation is given by (c1, s1) * (c2, s2) = (c1 * s1(c2), s1 * s2).

        Args:
            H (group): the direct product of n copies of some group.
            G (group): A group.

        Returns:
            group: The wreath product of H and G with the hyperoctahedral twist.
        """
        semidirect_product = GroupSemidirectProduct(
            G, N, act_to_right=True, twist=lambda s, c: self.hyperoctahedral_twist(s, c)
        )

        def random_element():
            Sn, C = semidirect_product.cartesian_factors()
            return semidirect_product((Sn.random_element(), C.random_element()))

        def order():
            return G.order() * N.order()

        def elements():
            return [semidirect_product((s, c)) for c in N for s in G]

        semidirect_product.random_element = random_element
        semidirect_product.action = self.hyperoctahedral_twist
        semidirect_product.order = order
        semidirect_product.elements = elements
        return semidirect_product

    def right_transversal(self, group, subgroup):
        transversal = [subgroup.one()]
        num_cosets = group.order() / subgroup.order()
        while len(transversal) < num_cosets:
            elt = group.random_element()
            if all(t * elt.inverse() not in subgroup for t in transversal):
                transversal.append(elt)
        return transversal

    def right_transversal_of_little_subgroup(self, k):
        from itertools import combinations

        n = self._n
        H = self
        combs = list(combinations(tuple((i for i in range(1, n + 1))), k))
        t_new = []
        for comb in combs:
            sets = ((i for i in range(1, k + 1)), (i for i in range(k + 1, n + 1)))
            elt = Permutation(
                [next(sets[0]) if i in comb else next(sets[1]) for i in range(1, n + 1)]
            )
            t_new.append(
                H._semidirect_product(
                    (H._Sn(elt.cycle_string()), (() for _ in range(n)))
                ).inverse()
            )
        return t_new

    def Phi(self, elt):
        c = tuple(self._C2(() if elt.inverse()(k) > 0 else (1, 2)) for k in range(1, self._n + 1))
        s = self._Sn(Permutation([abs(elt(k)) for k in range(1, self._n + 1)]))
        return self._semidirect_product((s, c))


# %%
