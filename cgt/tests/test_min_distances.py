from .. import example
from .. import distances
import pytest

def test_small_inversions_min_dist():
    for n in (2,3,4):
        framework, model = example(n)
        min_distances_graph_based = distances.min_distance(framework, model)
        min_distance_irrep_based = distances.min_distance_using_irreps(framework, model)
        for instance, distance in min_distances_graph_based.items():
            assert distance == min_distance_irrep_based[instance]