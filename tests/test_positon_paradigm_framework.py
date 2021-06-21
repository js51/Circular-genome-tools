import pytest

@pytest.fixture(scope='class')
def test_cases_conversions():
    cases = [
        { 
         #   'newick_string'  : "(0:0.1,1:0.4);",
         #   'number_of_taxa' : 2,
         #   'true_splits'    : [],
         #   'model'          : None,
        },
    ]
    return cases


def test_conversions(test_cases_conversions):
    """ Test a bunch of things """
    for case in test_cases_conversions:
        pass # assert some stuff
