from .. import example, distances, position_paradigm, Model
from .. import MODEL
import numpy as np
import pytest

def test_likelihood_function():
    for n in (3,4,5,6):
        framework = position_paradigm.GenomeFramework(n)
        converge_point = framework.symmetry_group().order()/framework.genome_group().order()
        model = Model.named_model_with_relative_probs(framework, {
            MODEL.all_inversions: 1
        })
        for _ in range(10):
            instance = framework.cycles(framework.canonical_instance(framework.random_instance()))
            identity_instance = framework.identity_instance()
            L = distances.likelihood_function(framework, model, instance)
            print(L(0))
            print(instance)
            if instance == identity_instance:
                assert np.isclose(L(0), 1)
            else:
                assert np.isclose(L(0), 0)
            assert np.isclose(float(L(500)), float(converge_point))


def test_small_inversions_min_dist():
    for n in (3,4):
        framework, model = example(n)
        min_distances_graph_based = distances.min_distance(framework, model)
        min_distance_irrep_based = distances.min_distance_using_irreps(framework, model)
        for instance, distance in min_distances_graph_based.items():
            assert distance == min_distance_irrep_based[instance]
