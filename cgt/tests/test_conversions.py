from .. import example
import pytest

@pytest.fixture
def example_framework():
    yield example()

def test_one_row_and_cycles_random(example_framework):
    fw, model = example_framework
    instance = fw.random_instance()
    assert instance == fw.cycles(fw.one_row(instance))

def test_matrix_and_cycles_random(example_framework):
    fw, model = example_framework
    instance = fw.random_instance()
    assert instance == fw.cycles(fw.matrix(instance))