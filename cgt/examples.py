"""
"""
from sage.modular.modform.element import ModularFormElement
from .enums import MODEL, SYMMETRY
from .models import Model
from .position_paradigm import PositionParadigmFramework

def example():
    framework = PositionParadigmFramework(6, symmetry=SYMMETRY.circular)
    model = Model.named_model_with_relative_probs(framework, {
        MODEL.one_region_inversions: 0.75, 
        MODEL.two_region_inversions: 0.25
    })
    return framework, model