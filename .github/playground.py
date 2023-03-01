
#%%
import cgt
from cgt import *
import cgt.simulations
import cgt.parsers

# %%

framework = cgt.PositionParadigmFramework(6)
model = cgt.Model.named_model_with_relative_probs(framework, {
    MODEL.all_inversions: 1
})

# %%
newick_string = "(((Z:1.1,a:1.3):1.3,B:1.2):1.2,((X:1.4,Y:1.9):1.8,D:1.7):1.1);"
tree = cgt.simulations.newick_to_tree(newick_string)

# %%
new_tree = cgt.simulations.evolve_on_tree(tree, framework, model)
cgt.simulations.draw_tree(new_tree)
# %%
