from .. import example
import pytest

@pytest.fixture
def example_framework():
    yield example(10)

def test_one_row_and_cycles_random(example_framework):
    fw, _ = example_framework
    instance = fw.random_instance()
    assert instance == fw.cycles(fw.one_row(instance))

def test_matrix_and_cycles_random(example_framework):
    fw, _ = example_framework
    instance = fw.random_instance()
    assert instance == fw.cycles(fw.matrix(instance))

def test_same_canonical_instance_random(example_framework):
    fw, _ = example_framework
    canonical_instance = fw.canonical_instance(fw.random_instance())
    genome = fw._genome_coset(fw.cycles(canonical_instance))
    for instance in genome:
        assert canonical_instance == fw.canonical_instance(instance)