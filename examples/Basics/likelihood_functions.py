# %%
import cgt
from cgt.enums import *
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
import line_profiler 
from line_profiler import * 
from decorator import decorator
from line_profiler import LineProfiler 

framework = cgt.PositionParadigmFramework(5, oriented=True, symmetry=SYMMETRY.circular)

model = cgt.Model.named_model_with_relative_probs(framework, {
    cgt.MODEL.one_region_inversions: 2/3,
    cgt.MODEL.two_region_inversions: 1/3
})
s = model.s_element(in_algebra=ALGEBRA.genome)
z = framework.symmetry_element()
zs = s*z


#%% 
%%time

from cgt.distances import likelihood_function

instances = [
    [1, 2, 3, 4, 5],
    [-5, 4, -1, -3, 2],
    [-5, -3, 2, 4, -1],
    [-1, -3, 5, -4, -2],
    [5, -4, -3, 2, 1],
    [-4, 1, 5, 3, 2],
    [-4, -3, -5, 2, -1],
    [-4, -2, 1, -5, 3],
    [4, -3, 2, -5, -1],
    [-2, -3, -5, 4, 1],
    [-3, 4, 5, -2, -1],
    [-4, 1, -2, -5, -3],
    [-4, -5, 2, -1, 3],
    [1, 3, -4, 2, 5],
    [-3, -2, -1, 4, -5],
    [3, 2, 5, 4, 1]
]
instance_gen = map(framework.cycles, instances)

dim = int(np.ceil(np.sqrt(len(instances))))
fig, ax = plt.subplots(dim, dim, sharex=False, sharey=False)
fig.set_size_inches(10, 8)
fig.tight_layout(h_pad=2.5, w_pad=1.5)

for i in range(dim):
    for j in range(dim):
        try:
            instance = next(instance_gen)
        except StopIteration:
            break
        L = likelihood_function(framework, model, instance, attempt_exact=False)
        times = np.arange(0, 25, 0.5)
        ax[i, j].plot(times, [L(t) for t in times], color='black')
        ax[i, j].set_title(str(framework.one_row(instance)))

# %%
