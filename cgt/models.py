"""
"""

from sage.all_cmdline import QQ
from .enums import *
from . import rearrangements
from scipy.sparse import dok_matrix as dok
import numpy as np


class Model:
    """Defines a model. A model consists of some collection of permutations and a map from these permutations to probabilities [0,1]"""
    def __init__(self, framework, generating_dictionary):
        """Define a model from a dictionary of single permutations, with their probabilities as the values."""
        self.framework = framework
        self.generating_dictionary = generating_dictionary
        self.names = []
        self._reg_rep_of_zs = None # For caching the regular representation of zs
        self.data_bundle = {} # Used to cache data for several functions
        # Check for time-reversibility
        #rearrangements = set(self.generating_dictionary.keys())
        #for rearrangement in rearrangements:
        #    if rearrangement.inverse() not in self.generating_dictionary:
        #        print(f"Model does not contain the inverse of {rearrangement}. Adding it #with equal probability so that the model is time-reversible. This means #that the generating dictionary is no longer the same as the one provided to #the constructor.")
        #        prob = self.generating_dictionary[rearrangement] / 2
        #        self.generating_dictionary[rearrangement.inverse()] = prob
        #        self.generating_dictionary[rearrangement] = prob

        # TODO: Implement checks for certain model properties, for example symmetry and missing rearrangements

    def __repr__(self):
        return f"Model({str(self.framework)}, {str(self.generating_dictionary)})"

    def __str__(self):
        """Return a descriptive string for the model"""
        # TODO: #21 Make the model able to describe which permutations it allows after it is created
        string = f"Model under the {str(self.framework).lower()}"
        string += f" containing {' and '.join(str(name.value) for name in self.names) if self.names else 'generating perms: ' + str(self.generating_dictionary)}."
        return string

    @classmethod
    def named_model_with_relative_probs(cls, framework, named_model_dictionary):
        """
        Create a model from a dictionary of named classes of rearrangements and their relative probabilities. 

        Args:
            framework (PositionParadigmFramework): The framework to create the model in
            named_model_dictionary (dict): A dictionary of named classes of rearrangements and their relative probabilities.

        Returns:
            Model: A model with the specified named classes of rearrangements and their relative probabilities.
        """
        if not framework.oriented:
            raise NotImplementedError(f"not yet implemented for {str(framework)}")
        # look through all specified model names and handle them appropriately
        model = {}
        if abs(sum(named_model_dictionary.values()) - 1) > 0.00001: 
            raise ValueError("supplied probabilities do not sum to 1")
        
        special_model = False
        # Check for special named models
        if MODEL.all_inversions_larger_less_likely in named_model_dictionary:
            special_model = True
            if len(named_model_dictionary) > 1:
                raise ValueError(f"Cannot have {MODEL.all_inversions_larger_less_likely} with other models")
            elements = rearrangements.fast_all_inversion_reps(framework)
            elements = sorted(list(elements), key = lambda x: len(x.cycles()))
            denom = 2**len(elements) - 1
            for e, element in enumerate(elements):
                model[element] = (2**(len(elements) - e - 1))/denom

        if not special_model:
            for model_name, relative_prob in named_model_dictionary.items():
                match model_name:
                    # Inversions
                    case MODEL.all_inversions:
                        gens = rearrangements.all_inversions_representatives(framework)
                    case MODEL.one_region_inversions:
                        gens = rearrangements.all_inversions_representatives(framework, num_regions=1)
                    case MODEL.two_region_inversions:
                        gens = rearrangements.all_inversions_representatives(framework, num_regions=2)
                    # Transpositions
                    case MODEL.all_transpositions:
                        gens = rearrangements.all_transposition_instances(framework, canonical_reps_only = True)
                    case MODEL.two_region_transpositions:
                        gens = rearrangements.all_transposition_instances(framework, scope_limit = 2, canonical_reps_only = True)
                    case MODEL.two_region_transpositions_without_inversions:
                        gens = rearrangements.all_transposition_instances(framework, scope_limit = 2, with_inversion = False, canonical_reps_only = True)
                    case MODEL.two_region_revrevs:
                        gens = rearrangements.all_transposition_instances(framework, scope_limit = 2, only_revrevs = True, canonical_reps_only = True)
                    case MODEL.one_region_moves:
                        gens = rearrangements.all_transposition_instances(framework, single_segment_limit = 1, canonical_reps_only = True)
                    case MODEL.one_region_moves_without_inversions:
                        gens = rearrangements.all_transposition_instances(framework, single_segment_limit = 1, with_inversion = False, canonical_reps_only = True)
                    case MODEL.three_region_transpositions:
                        gens = rearrangements.all_transposition_instances(framework, scope_limit = 3, canonical_reps_only = True)
                    case MODEL.two_region_adjacent_transpositions:
                        gens = rearrangements.all_adjacent_transpositions_representatives(framework, num_regions=2)
                for generator in gens: 
                    model[generator] = relative_prob/len(gens)

        model = cls(framework, model)
        model.names += list(named_model_dictionary.keys())
        return model

    def reg_rep_of_zs(self):
        """Return the regular representation of zs as comptued by PositionParadigmFramework.reg_rep_zs, but store the sparse result"""
        if self._reg_rep_of_zs is None:
            self._reg_rep_of_zs = self.framework.reg_rep_of_zs(self, sparse=True)
        return self._reg_rep_of_zs
    
    def s_element(self, in_algebra=ALGEBRA.genome):
        if in_algebra not in {ALGEBRA.group, ALGEBRA.genome}:
            raise NotImplementedError(f"Model element for {str(in_algebra)} algebra not yet implemented")
        A = self.framework.group_algebra()
        s = A(0) # Zero element in algebra
        gens_dict = self.generating_dictionary
        gens = {x for x in self.generating_dictionary.keys()}
        while len(gens) > 0:
            gen = gens.pop()
            prob = gens_dict[gen]
            if in_algebra is ALGEBRA.group:
                conj_class = rearrangements.conjugacy_class(self.framework, gen)
            elif in_algebra is ALGEBRA.genome:
                conj_class = {gen}
            for perm in conj_class:
                s += QQ(prob/len(conj_class)) * A(perm)
            gens -= conj_class
        return s