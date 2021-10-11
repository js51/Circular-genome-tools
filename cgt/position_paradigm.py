"""
position_paradigm is the most important module in cgt. PositionParadigmFramework keeps track of important objects, can
convert group elements and genomes between different forms, and instances can be passed to functions from other modules so
that they can acces information without it needing to be recreated.
"""

from sage.all_cmdline import QQ, FormalSum, SymmetricGroup, SymmetricGroupRepresentations, gap, UniversalCyclotomicField, Permutations, matrix, QQbar
from sage.combinat.colored_permutations import SignedPermutations
from .enums import *
import numpy as np
import warnings
from copy import deepcopy
from .structures import HyperoctahedralGroup
from scipy.sparse import dok_matrix as dok
from random import choice
from functools import cache

class PositionParadigmFramework:
    """Everything you need for working with genomes under the position paradigm"""

    def __init__(self, num_regions, oriented=True, symmetry=SYMMETRY.circular, genome_type=TYPE.reg_to_signed_pos):
        """Instantiate a system of position paradigm genomes with given set of symmetries and number of regions."""
        if genome_type != TYPE.reg_to_signed_pos:
            raise NotImplementedError("Genome type not yet supported. For now use TYPE.reg_to_signed_pos")
        self.n = num_regions
        self.oriented = oriented
        self.symmetry = symmetry
        
    def __eq__(self, other): # String representation is unique and repr calls str
        return self.__repr__() == other.__repr__()

    def __str__(self):
        string = f"Framework for {str(self.symmetry).replace('SYMMETRY.', '')} genomes"
        string += f" with {self.n}{' oriented' if self.oriented else ''} regions"
        return string

    def __repr__(self):
        return self.__str__() # String representation is unique

    def __call__(self, x):
        """Return an object as a group algebra element if possible"""
        return self.genome_algebra()(self.cycles(x))

    def __hash__(self):
        return self.__str__().__hash__() # String rep is unique, so just hash that!

    def cycles(self, element):
        try: # See if it's already a group element
            return self.genome_group()(element)
        except ValueError:
            pass # Not able to be directly converted, try something else!
        if type(element) is np.ndarray: # If it's a matrix, let's convert it!
            return self._permutation_group_element_from_matrix(element)
        else: # Not a matrix, assume it's a list
            try:
                elt = list(element)
            except TypeError:
                raise ValueError("Invalid data for constructing a permutation")
        if (not self.oriented) or len(elt) == 2*self.n:
            return self.genome_group()(elt)
        elif len(elt) == self.n:
            other_half = list(-int(elt[self.n - i]) for i in range(1, self.n+1))
            return self.genome_group()(elt + other_half)
        else:
            raise ValueError("Invalid data for constructing a permutation")

    def one_row(self, element, as_list=False):
        """Return a given genome instance in one-row notation.

        Args:
            element: the genome instance to represent in one-row notation. Multiple formats accepted.
            as_list: if true, return as a list of images, otherwise return a sage permutation.

        EXAMPLES::
            sage: import cgt
            sage: ppf = cgt.PositionParadigmFramework(3)
            sage: ppf.cycles(ppf.one_row('(1,-2)(-1,2)')) == ppf.cycles('(1,-2)(-1,2)')
            True
        """
        elt = self.cycles(element)
        row = list(elt.dict().values())[0:self.n]
        if not as_list:
            if self.oriented:
                row = SignedPermutations(self.n)(row)
            else:
                row = Permutations(self.n)(row)
        return row

    @cache
    def genome_group(self):
        """Return the permutation group containing genome instances."""
        if self.oriented:
            return HyperoctahedralGroup(self.n)
        else:
            return SymmetricGroup(self.n)

    def canonical_instance(self, instance):
        """Return the 'canonical' instance of the genome represented by the permutation if there is one.

        TESTS::
            sage: import cgt; ppf = cgt.PositionParadigmFramework(3, symmetry=cgt.SYMMETRY.circular)
            sage: ppf.canonical_instance('(1,-2)(-1,2)')
            [1, 2, -3]
        """
        instance = deepcopy(self.one_row(self.cycles(instance)))
        if self.symmetry in {SYMMETRY.circular, SYMMETRY.linear}:
            f = self.one_row(self.standard_reflection())
            if list(instance)[0] < 0:
                instance = instance*f
            if self.symmetry is SYMMETRY.circular:
                r = self.one_row(self.standard_rotation())
                instance = instance*(r**((self.n-list(instance)[0]+1)%self.n))
            return instance
        else:
            raise NotImplementedError("No canonical form exists in the current framework")

    def random_instance(self, genome=None):
        """Return a random permutation corresponding to a given genome, or a random genome if none is given."""
        if genome:
            return choice(list(genome))[0]
        else:
            return self.genome_group().random_element()

    def random_genome(self, format=FORMAT.formal_sum):
        """Return a random genome"""
        return self.genome(self.random_instance(), format=format)

    @cache
    def symmetry_group(self):
        """Return the symmetry group of the genomes."""
        if self.symmetry == SYMMETRY.circular:
            gens = [self.standard_reflection(), self.standard_rotation()]
        elif self.symmetry == SYMMETRY.linear:
            gens = [self.standard_reflection()]
        else:
            gens = [self.genome_group().one()]
        return self.genome_group().subgroup(gens)

    @cache
    def group_algebra(self):
        """Return the group alegbra, where the group is the group containing genome instances."""
        return self.genome_group().algebra(QQ)

    def symmetry_element(self):
        """Return the symmetry element: a convex sum of elements from symmetry_group()."""
        try:
            return self.z
        except AttributeError:
            self.z = 1/self.symmetry_group().order() * sum(self.group_algebra()(d) for d in self.symmetry_group())
        return self.z

    def num_genomes(self):
        """Return the number of distinct genomes up to symmetries."""
        return self.genome_group().order()/self.symmetry_group().order()

    def _sort_key(self, one_row_perm):
        return str(one_row_perm).replace('-', 'Z')

    def _sort_key_cycles(self, cycles_perm):
        return str(self.one_row(cycles_perm)).replace('-', 'Z')

    def _genome_coset(self, instance):
        coset = set(instance*d for d in self.symmetry_group())
        return sorted([self.one_row(g) for g in coset], key=self._sort_key)

    def _double_coset(self, instance):
        coset = set(d2 * instance * d1 for d1 in self.symmetry_group() for d2 in self.symmetry_group())
        return sorted(coset, key=self._sort_key_cycles)

    def genome_equivalence_classes(self, combine_inverse_classes=False, sort_classes=True):
        """Return double cosets of instances---genomes in the same class will have the same likelihood functions"""
        instances = set(self.genome_group())
        classes = {}
        while len(instances) > 0:
            instance = instances.pop()
            if combine_inverse_classes:
                dcoset = sorted(set(self._double_coset(instance)+self._double_coset(instance.inverse())), key=self._sort_key_cycles)
            else:
                dcoset = self._double_coset(instance)
            instances -= set(dcoset)
            classes[dcoset[0]] = dcoset
        if sort_classes:
            classes = dict(sorted(classes.items(), key=lambda x: self._sort_key_cycles(x[0])))
        return classes

    def genome(self, instance, format=None):
        coset = self._genome_coset(instance)
        if format == FORMAT.formal_sum:
            Z = self.symmetry_group()
            A = self.group_algebra()
            return sum(1/Z.order()*A(self.cycles(dx)) for dx in coset)
        else:
            return coset

    def genomes(self, format=FORMAT.dictionary, sort_genomes=True):
        if format not in {FORMAT.dictionary, FORMAT.formal_sum}:
            raise NotImplementedError("Not yet implemented! Convert manually to other formats from FORMAT.dictionary")
        instances = set(self.genome_group())
        genomes = {}
        while len(instances) > 0:
            instance = instances.pop()
            coset = self._genome_coset(instance)
            instances -= set(coset)
            genomes[coset[0]] = coset
        if format == FORMAT.formal_sum:
            Z = self.symmetry_group()
            A = self.group_algebra()
            genomes = { rep : sum(1/Z.order()*A(self.cycles(dx)) for dx in coset) for rep, coset in genomes.items() }
        if sort_genomes:
            genomes = dict(sorted(genomes.items(), key=lambda x: self._sort_key(x[0])))
        return genomes

    def standard_reflection(self):
        """Return a permutation which reflects a genome instance about the center region (n odd), or the center two regions (n even)."""
        if self.oriented:
            string = f'({1},-{self.n})({self.n},-{1})'
            for i in range(1, int(self.n/2)):
                string += f'({i+1}, -{self.n-(i-1)-1})({self.n-(i-1)-1}, -{i+1})'
            if self.n % 2 != 0:
                mid = int(self.n/2) + 1
                string += f'({mid},-{mid})'
        else:
            string = f'(1,{self.n})'
            for i in range(1, int(self.n/2)):
                string += f'({i+1},{self.n-(i-1)-1})'
        return self.genome_group()(string)

    def standard_rotation(self):
        """Return a permutation which rotates positions clockwise by one position."""
        positive = tuple(i for i in range(1,self.n+1))
        string = str(positive)
        if self.oriented:
            negative = tuple(-i for i in positive)
            string += str(negative)
        string = string.replace(',)', ')')
        return self.genome_group()(string)

    def draw_instance(self, instance, shortened=False):
        """Return a one-line drawing of a genome using arrows indicating orientation and elipsis indicating rotational symmetry"""
        permutation = instance.inverse()
        if shortened:
            left_tail, right_head = '\u2758', '\u276D'
            left_head, right_tail = '\u276C', '\u2758'
        else:
            left_tail, right_head = '|-', '->'
            left_head, right_tail = '<-', '-|'
        if self.symmetry in {SYMMETRY.circular, SYMMETRY.linear}:
            string = "..." if self.symmetry is SYMMETRY.circular else ""
            for position in range(1, self.n + 1):
                signed_region = permutation(position)
                if signed_region < 0:
                    region = f"{left_head}{abs(signed_region)}{right_tail}"
                else:
                    region = f"{left_tail}{abs(signed_region)}{right_head}"
                if position == 1:
                    first_region = region
                string += region
            string += f"{first_region}..." if self.symmetry is SYMMETRY.circular else ""
            string = string if self.symmetry is SYMMETRY.circular else string[1:]
            return string if self.oriented else string.replace(right_head, right_tail).replace(left_head, left_tail)
        else:
            raise NotImplementedError(f"Can't draw genome instance with symmetry group {str(self.symmetry)}")

    def _permutation_group_element_from_matrix(self, matrix):
        sigma = []
        for c, col in enumerate(matrix.T):
            image = np.nonzero(col)[0][0]
            sigma.append(matrix[image,c]*(image+1))
        return self.cycles(sigma)

    def matrix(self, element):
        """Return the defining representation matrix. Note that matrix(x)matrix(y) = matrix(yx)"""
        with warnings.catch_warnings(): # Sage uses deprecated objects in to_matrix for coloured permutations
            warnings.simplefilter("ignore", category=PendingDeprecationWarning)
            return np.array(self.one_row(self.cycles(element)).to_matrix())

    def reg_rep_of_zs(self, model, to_adjacency_matrix=False, sparse=True):
        warnings.warn("Untested! Use at your own risk!")
        # TODO: #14 re-write to directly use model element from the genome algebra
        if self is not model.framework:
            if self != model.framework:
                raise ValueError(f"Current framework and model framework are not the same!")
            else:
                warnings.warn("Current framework and model framework reference different objects! This might cause problems.")
        Z = self.symmetry_group()
        genomes = self.genomes()
        genome_list = list(genomes.values())
        reps = list(genomes.keys()) # Thousands of these
        #model_classes = list({frozenset({ d.inverse() * a * d for d in Z }) for a in model.generating_dictionary.keys()})
        model_generators_cycles = list(model.generating_dictionary.keys()) #[sorted(list(model_class))[0] for model_class in model_classes]
        model_generators = [
            self.one_row(elt) for elt in model_generators_cycles
        ]
        #model_classes = { rep : list(model_classes[r]) for r, rep in enumerate(model_generators) }
        genome_lookup = { rep : r for r, rep in enumerate(reps) }
        num_genomes = len(genomes.keys())
        matrix = dok((num_genomes, num_genomes), dtype=np.float32)
        for g, genome in enumerate(genome_list): # For each genome
            print(f'\rComputing class {g+1} with rep {genome[0]}', end="")
            for permutation in genome:
                for rearrangement in model_generators:
                    coeff = 1
                    if not to_adjacency_matrix:
                        coeff = 1/Z.order() * model.generating_dictionary[model_generators_cycles[model_generators.index(rearrangement)]] # * len(model_classes[rearrangement])
                    matrix[g, genome_lookup[self.canonical_instance(permutation * rearrangement)]] += coeff
        print('...done!')
        return matrix if sparse else matrix.toarray()

    def irreps(self, element=None):
        """
        Return a complete list of pairwise irreducible representations of the genome group.

        Args:
            element: return the image of element for each irreducible representation.

        Returns:
            a complete list of pairwise irreducible representations of the genome group, or a list of matrices if element is not None. 
            Output is cached so that irreps are computed only once per session.
            Matrix entries sit inside the UniversalCyclotomicField and are exact expressions.
        """
        representations = self._cached_irreps()
        if element is not None:
            return [irrep(element) for irrep in representations]
        else:
            return representations

    @cache
    def _cached_irreps(self):
        representations = []
        def irrep_function_factory(irrep, signed):
            def representation(sigma, _irrep=irrep, _signed=signed):
                    result = 0
                    if sigma in self.group_algebra(): # sigma is an algebra element
                        for term in sigma:
                            perm, coeff = term
                            result += coeff * (gap.Image(_irrep, perm) if _signed else _irrep(perm))
                    else: # sigma is a group element
                        result = (gap.Image(_irrep, sigma) if _signed else _irrep(sigma))
                    mat = matrix(UniversalCyclotomicField(), result)
                    return mat.transpose() if _signed else mat
            return representation
        if not self.oriented:
            irreps = SymmetricGroupRepresentations(self.n)
        else:
            irreps = (gap.IrreducibleAffordingRepresentation(character) for character in gap.Irr(self.genome_group()))
        for irrep in irreps:
            representations.append(irrep_function_factory(irrep, self.oriented))
        return representations

    def regular_representation(self, g):
        """Return the regular representation of a single element"""
        warnings.warn("this function is untested! Use at your own risk.", DeprecationWarning)
        return matrix(QQbar, [(self.group_algebra()(g)*self.genome_group()(h)).to_vector(QQbar) for h in self.genome_group()]).transpose()

    def coefficient_in(self, more_terms, fewer_terms):
        """For example, the coefficient of (a/2 + b/2) in x=(a/3 + b/3 + c/3) is 2/3, since x=2/3*(a/2 + b/2) + c/3"""
        coefficients_in_larger_sum = {}
        for term in fewer_terms.terms():
            perm  = self.genome_group()(list(term)[0][0]) # a
            coeff = list(term)[0][1] # 1/2
            coefficients_in_larger_sum[perm] = more_terms.coefficient(perm) * (1/coeff) # 1/3 * 2 = 2/3
        terms = set(coefficients_in_larger_sum.values())
        if len(terms) != 1:
            warnings.warn("Note that extra terms from the smaller sum remain in the larger sum!")
        return min(terms)

    def collect_genome_terms(self, formal_sum, display=DISPLAY.one_row):
        """Return a weighted sum of genomes, represented as [canonical representation]z"""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            formal_sum = deepcopy(formal_sum)
        terms = formal_sum.terms()
        new_sum = [] # [(coeff, rep), ...]
        while len(terms) > 0:
            term = list(terms[0])[0][0]
            canonical_instance = self.canonical_instance(term)
            genome = self.genome(self.cycles(canonical_instance), format=FORMAT.formal_sum)
            coeff = self.coefficient_in(formal_sum, genome)
            formal_sum -= coeff * genome
            terms = formal_sum.terms()
            if display is DISPLAY.one_row:
                genome_string = str(canonical_instance).replace(' ', '') + 'z'
            elif display is DISPLAY.cycles:
                genome_string = str(self.cycles(canonical_instance)) + 'z'
            elif display is DISPLAY.arrows:
                genome_string = '(' + self.draw_instance(self.cycles(canonical_instance), shortened=True) + ')'
            new_sum.append((coeff, genome_string))
        return FormalSum(new_sum, parent=QQ)
