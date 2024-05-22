from .. import position_paradigm
from .. import hyperoctahedral

def test_hyperoctahedral_group_irreps():
    for n in (2,3,4,5):
        framework = position_paradigm.GenomeFramework(n)
        for irrep in framework.irreps():
            r1 = framework.random_instance()
            r2 = framework.random_instance()
            assert (irrep(r1) * irrep(r2)) == irrep(r1*r2)

def test_hyperoctahedral_group_map():
    for n in (3, 4, 5, 10, 20, 30):
        framework = position_paradigm.GenomeFramework(n)
        H = hyperoctahedral.HyperoctahedralGroup(n)
        for _ in range(10):
            r1 = framework.random_instance()
            r2 = framework.random_instance()
            assert H.Phi(r1) * H.Phi(r2) == H.Phi(r1*r2)

