from .. import example

def test_example():
    fw, model = example()
    assert model.framework is fw