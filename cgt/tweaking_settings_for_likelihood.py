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

framework = cgt.PositionParadigmFramework(6, oriented=True, symmetry=SYMMETRY.circular)

model = cgt.Model.named_model_with_relative_probs(framework, {
    cgt.MODEL.one_region_inversions: 2/3,
    cgt.MODEL.two_region_inversions: 1/3
})
s = model.s_element(in_algebra=ALGEBRA.genome)
z = framework.symmetry_element()
zs = s*z


#%% 
%%time

from threading import Thread

from cgt.distances import likelihood_function

instances = [
    [1,2,3,4,  5,-6],
    [1,2,-3,4, 5,-6],
    [1,2,3,4,  5,-6],
    [1,2,-3,4, 5,-6],
]
instance_gen = map(framework.cycles, instances)

dim = int(np.ceil(np.sqrt(len(instances))))
fig, ax = plt.subplots(dim, dim, sharex=False, sharey=False)
fig.tight_layout()
threads = []
for i in range(dim):
    for j in range(dim):
        try:
            instance = next(instance_gen)
        except StopIteration:
            break
        def do_L(i,j):
            L = likelihood_function(framework, model, instance, attempt_exact=False)
            times = np.arange(0, 2100, 0.1)
            ax[i, j].plot(times, [L(t) for t in times], color='black')
        thread = Thread(target=do_L, args=(i,j))
        thread.start()
        threads.append(thread)

for thread in threads:
    thread.join()

fig.show()

# %%
