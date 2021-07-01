"""
"""
from .enums import MODEL, SYMMETRY
from .models import Model
from .position_paradigm import PositionParadigmFramework

def example(n):
    framework = PositionParadigmFramework(n, symmetry=SYMMETRY.circular)
    model = Model.named_model_with_relative_probs(framework, {
        MODEL.one_region_inversions: 0.5, 
        MODEL.two_region_inversions: 0.5
    })
    return framework, model