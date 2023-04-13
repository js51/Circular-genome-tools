
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

# %% Get the genomes g that we need dist(id -> g) for:
need_distances = {}
for i, g_i in enumerate(genomes):
    for j, g_j in enumerate(genomes):
        if i < j:
            canonical_i = framework.canonical_instance(framework.random_instance(g_i))
            canonical_j = framework.canonical_instance(framework.random_instance(g_j))
            need_distances[(canonical_i, canonical_j)] = canonical_i.inverse() * canonical_j

# %% Compute distance matrix
min_distances = cgt.distances.min_distance(framework, model, weighted=False)

# %% Likelihood distances
L_distances = cgt.distances.mles(framework, model, genome_instances=need_distances.values(), verbose=True)

# %%
MFPT_distances = cgt.distances.MFPT(framework, model)

#%%
D_min = cgt.distances.dict_to_distance_matrix(min_distances, framework, genomes)
D_MFPT = cgt.distances.dict_to_distance_matrix(MFPT_distances, framework, genomes)

# %%
tree1 = sp.phylogenetics.neighbour_joining(D_min, labels=labels)
cgt.simulations.draw_tree(tree1)

# %%
tree2 = sp.phylogenetics.neighbour_joining(D_MFPT, labels=labels)
cgt.simulations.draw_tree(tree2)
# %%
