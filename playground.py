
#%%
import cgt
from cgt import *
import cgt.simulations
import cgt.parsers
import splitp as sp
import numpy as np
from splitp import phylogenetics
from cgt import distances

# %%
framework = cgt.PositionParadigmFramework(5)
model = cgt.Model.named_model_with_relative_probs(framework, {
    MODEL.all_inversions: 1
})

# %%
newick_string = "(((Z:1.1,a:1.3):1.3,B:1.2):1.2,((X:1.4,Y:1.9):1.8,D:1.7):1.1);"
tree = cgt.simulations.newick_to_tree(newick_string)

# %%
new_tree = cgt.simulations.evolve_on_tree(tree, framework, model)
cgt.simulations.draw_tree(new_tree)

# %% Get the genomes at the leaves of the tree
leaves = [n for n, d in new_tree.out_degree() if d==0]
genomes = [new_tree.nodes[leaf]["genome"] for leaf in leaves]
labels = [new_tree.nodes[leaf]["label"] for leaf in leaves]
distances = cgt.distances.min_distance(framework, model, weighted=True)

#%%
# Convert list of distances to a distance matrix
D = np.zeros((len(genomes), len(genomes)))
for i in range(len(genomes)):
    canonical_i = framework.canonical_instance(framework.random_instance(genomes[i]))
    distances_copy = { canonical_i * k : v for k, v in distances.items() }
    for j in range(len(genomes)):
        canonical_j = framework.canonical_instance(framework.random_instance(genomes[j]))
        D[i,j] = distances_copy[canonical_j]

# %%
tree = sp.phylogenetics.neighbour_joining(D, labels=labels)
cgt.simulations.draw_tree(tree)
cgt.simulations.draw_tree(new_tree)

# %%
import numpy as np
import networkx as nx
import splitp as sp
from splitp import phylogenetics
import cgt
from cgt import simulations

matrix = np.array([[0,  5,  9,  9, 8],
                   [5,  0, 10, 10, 9],
                   [9, 10,  0,  8, 7],
                   [9, 10,  8,  0, 3],
                   [8,  9,  7,  3, 0]])

tree = sp.phylogenetics.neighbour_joining(matrix, labels=["A", "B", "C", "D", "E"])
sp.phylogenetics.midpoint_rooting(tree)
sp.parsers.newick.move_tree_edge_labels_to_nodes(tree)
json = nx.tree_data(tree, root=-1)
sp.parsers.newick.json_to_newick(json, lengthstring="weight")
# %%
