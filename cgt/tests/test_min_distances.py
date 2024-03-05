from .. import example
from .. import distances
import pytest

@pytest.fixture
def example_framework():
    yield example(4)

def test_small_inversions_min_dist(example_framework):
    framework, model = example_framework
    genome = framework.random_instance()
    min_distances_graph_based = distances.min_distance(framework, model)
    min_distance_irrep_based = distances.min_distance_using_irreps(framework, model)
    for instance, distance in min_distances_graph_based.items():
        assert distance == min_distance_irrep_based[instance]