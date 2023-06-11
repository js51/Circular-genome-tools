
# %% Imports
import cgt
from cgt import *
import cgt.simulations
import cgt.parsers
import splitp as sp
import numpy as np
from splitp import phylogenetics
from cgt import distances

# %% Define framework and model
framework = cgt.PositionParadigmFramework(5)
model = cgt.Model.named_model_with_relative_probs(framework, {
    MODEL.all_inversions: 1
})

# %% Define tree
newick_string = "(((Z:1.1,A:1.3):1.3,B:1.2):1.2,((X:1.4,Y:1.9):1.8,D:1.7):1.1);"
tree = cgt.simulations.newick_to_tree(newick_string)

#%% Draw tree with SplitP
import splitp as sp
sp.Phylogeny(newick_string).draw()

# %%
new_tree = cgt.simulations.evolve_on_tree(tree, framework, model)
cgt.simulations.draw_tree(new_tree)

# %% Get the genomes at the leaves of the tree
leaves = [n for n, d in new_tree.out_degree() if d==0]
genomes = [new_tree.nodes[leaf]["genome"] for leaf in leaves]
labels = [new_tree.nodes[leaf]["label"] for leaf in leaves]

# %% Compute MLE distance matrix
D_min = cgt.distances.distance_matrix(framework, model, genomes, cgt.DISTANCE.min)

# %% Compute MLE distance matrix
D_MLE = cgt.distances.distance_matrix(framework, model, genomes, cgt.DISTANCE.MLE)

# %%
reconstructed_tree = sp.phylogenetics.neighbour_joining(D_MLE, labels=labels)
cgt.simulations.draw_tree(reconstructed_tree)
# %%
