from enum import Enum, auto
import numpy as np

class SYMMETRY(Enum):
    # Circular genomes
    circular    = auto()
    dihedral    = circular
    D_n         = circular
    # Linear genomes
    linear      = auto()
    S_2         = linear
    C_2         = linear

class MODEL(Enum):
    pass

class TYPE(Enum):
    reg_to_signed_pos = auto()
    pos_to_signed_reg = auto()
    signed_reg_to_pos = auto()
    signed_pos_to_reg = auto()

class FORMAT(Enum):
    formal_sum      = auto()
    equiv_classes   = auto()
    dictionary      = auto()
    only_reps       = auto()

def HyperoctahedralGroup(n):
    pmn = tuple(range(1,n+1)) + tuple(range(-n, 0))
    s_pmn = SymmetricGroup(pmn)
    if n==1: return s_pmn
    h_gens = [
        s_pmn([(1,2),(-1,-2)]),
        s_pmn(tuple(range(1,n+1)) + tuple(range(-1, -(n+1), -1)))
    ]
    return s_pmn.subgroup(h_gens)

class PositionParadigmFramework:
    """Everything you need for working with genomes under the position paradigm"""

    def __init__(self, num_regions, oriented=True, symmetry=None, genome_type=TYPE.reg_to_signed_pos):
        """Instantiate a system of position paradigm genomes with given set of symmetries and number of regions."""
        if genome_type != TYPE.reg_to_signed_pos:
            raise NotImplementedError("Genome type not yet supported. For now use TYPE.reg_to_signed_pos")
        self.n = num_regions
        self.oriented = oriented
        self.symmetry = symmetry

    def __str__(self):
        string = f"Framework for {str(self.symmetry).replace('SYMMETRY.', '')} genomes"
        string += f" with {self.n}{' oriented' if self.oriented else ''} regions"
        return string

    def __repr__(self):
        return self.__str__()

    ### Conversion functions

    def __call__(self, x):
        """Return an object as a group or group algebra element if possible, oherwise return None"""
        if x in self.genome_group():
            return self.genome_group()(x)
        elif type(x) is list:
            return self.cycles(x)
        elif x in self.group_algebra():
            return self.group_algebra()(x)
        else:
            return None

    def cycles(self, element):
        try: # See if it's already a group element
            return self.genome_group()(element)
        except ValueError:
            pass # Not able to be directly converted, try something else!
        if type(element) is np.ndarray: # If it's a matrix, let's convert it!
            return self.__permutation_group_element_from_matrix(element)
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

        Examples:
            >>> PositionParadigmFramework(3).one_row('(1,-2)(-1,2)')
            [-2, -1, 3]
        """
        elt = self.cycles(element)
        row = list(elt.dict().values())[0:self.n]
        if not as_list:
            if self.oriented:
                row = SignedPermutations(self.n)(row)
            else:
                row = Permutations(self.n)(row)
        return row

    ### Algebraic structures

    def genome_group(self):
        """Return the permutation group containing genome instances."""
        try:
            return self.G
        except AttributeError:
            if self.oriented:
                self.G = HyperoctahedralGroup(self.n)
            else:
                self.G = SymmetricGroup(self.n)
        return self.G

    def symmetry_group(self):
        """Return the symmetry group of the genomes."""
        try:
            return self.Z
        except AttributeError:
            if self.symmetry == SYMMETRY.circular:
                gens = [self.standard_reflection(), self.standard_rotation()]
            elif self.symmetry == SYMMETRY.linear:
                gens = [self.standard_reflection()]
            else:
                gens = [self.genome_group().one()]
            self.Z = self.G.subgroup(gens)
        return self.Z
    
    def group_algebra(self):
        """Return the group alegbra, where the group is the group containing genome instances."""
        try:
            return self.CG
        except AttributeError:
            self.GA = self.genome_group().algebra(CDF)
        return self.GA
        
    def symmetry_element(self):
        """Return the symmetry element: a convex sum of elements from symmetry_group()."""
        try:
            return self.z
        except AttributeError:
            self.z = 1/self.symmetry_group().order() * sum(self.group_algebra()(d) for d in self.symmetry_group())
        return self.z

    def set_model(self):
        raise(NotImplementedError())
    
    def num_genomes(self):
        """Return the number of distinct genomes up to symmetries."""
        return self.genome_group().order()/self.symmetry_group().order()
    
    def sort_key(self, one_row_perm):
        return str(one_row_perm).replace('-', 'Z')

    def genomes(self, format=FORMAT.dictionary, sort_genomes=True):
        if format not in {FORMAT.dictionary, FORMAT.formal_sum}:
            raise NotImplementedError("Not yet implemented! Convert manually to other formats from FORMAT.dictionary")
        instances = set(self.genome_group())
        genomes = {}
        while len(instances) > 0:
            instance = instances.pop()
            coset = set(instance*d for d in self.symmetry_group())
            instances -= coset
            coset = sorted([self.one_row(g) for g in coset], key=self.sort_key)
            genomes[coset[0]] = coset
        if format == FORMAT.formal_sum:
            Z = self.symmetry_group()
            A = self.group_algebra()
            genomes = { rep : sum(1/Z.order()*A(self.cycles(dx)) for dx in coset) for rep, coset in genomes.items() }
        if sort_genomes:
            genomes = dict(sorted(genomes.items(), key=lambda x: self.sort_key(x[0])))
        return genomes

    def num_rearrangements(self):
        raise(NotImplementedError())

    def standard_reflection(self):
        """Return a permutation which reflects a genome instance about the center region (n odd), or 
        the center two regions (n even).
        """
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

    def __permutation_group_element_from_matrix(self, matrix):
        sigma = []
        for c, col in enumerate(matrix.T):
            image = np.nonzero(col)[0][0]
            sigma.append(matrix[image,c]*(image+1))
        return self.cycles(sigma)

    def matrix(self, element):
        """Return the defining representation matrix. Note that matrix(x)matrix(y) = matrix(yx)"""
        return np.array(self.one_row(self.cycles(element)).to_matrix())

    def make_inversion(self, a, b):
        raise(NotImplementedError())