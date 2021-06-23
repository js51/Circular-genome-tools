"""
"""

from enum import Enum, auto
from sage.all_cmdline import *
from .enums import *
from . import rearrangements
from . import structures

class Model:
	"""Defines a model. A model consists of some collection of permutations and a map from these permutations to probabilities [0,1]"""
	def __init__(self, framework, generating_dictionary):
		"""Define a model from a dictionary of single permutations, with their probabilities as the values."""
		self.framework = framework
		self.generating_dictionary = generating_dictionary
		# TODO: Check model properties

	@classmethod
	def named_model_with_relative_probs(cls, framework, named_model_dictionary):
		if not framework.oriented:
			raise NotImplementedError(f"not yet implemented for {str(framework)}")
		# look through all specified model names and handle them appropriately
		model = {}
		if abs(sum(named_model_dictionary.values()) - 1) > 0.00001: 
			raise ValueError("supplied probabilities do not sum to 1")
		for model_name, relative_prob in named_model_dictionary.items():
			if model_name is MODEL.one_region_inversions:
				gens = rearrangements.all_inversions_representatives(framework, num_regions=1)
				for generator in gens:
					model[generator] = relative_prob/len(gens)
			if model_name is MODEL.two_region_inversions:
				gens = rearrangements.all_inversions_representatives(framework, num_regions=2)
				for generator in gens:
					model[generator] = relative_prob/len(gens)
			if model_name is MODEL.all_inversions:
				gens = rearrangements.all_inversions_representatives(framework)
				for generator in gens:
					model[generator] = relative_prob/len(gens)
		return cls(framework, model)