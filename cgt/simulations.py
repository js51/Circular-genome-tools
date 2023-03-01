#%%
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import cgt
from cgt import *
from networkx.drawing.nx_agraph import graphviz_layout
from networkx import json_graph
import networkx as nx

def newick_to_tree(newick_string):
    json_tree = cgt.parsers.newick_to_json(
        newick_string, 
        generate_names = True, 
        lengthstring = "weight"
    )
    tree = json_graph.tree_graph(json_tree)
    return tree

def evolve_on_tree(tree, framework, model, root="random"):
    model_element = model.s_element()
    # The tree should be a networkx DiGraph
    # Choose a genome at the root of the tree
    if root == "random": 
        root = framework.random_genome()
    else: 
        root = framework.identity_genome()
    
    # Get the root and set the genome at the root
    root_id = [n for n, d in tree.in_degree() if d==0][0]
    tree.nodes[root_id]["genome"] = root

    # Traverse the tree
    root = [n for n,d in tree.in_degree() if d==0][0]
    for successor, node in dict(nx.bfs_predecessors(tree, root)).items():
        # Get the genome at the parent node
        parent = tree.nodes[node]["genome"]
        # Get the branch length of the in-edge
        branch_length = tree.nodes[successor]["weight"]
        tree[node][successor]['weight'] = branch_length
        tree[node][successor]['len'] = branch_length
        # Decide how many rearrangements to apply
        num_rearrangements = np.random.poisson(branch_length)
        # Apply them
        tree.nodes[successor]["genome"] = parent
        if num_rearrangements == 0: 
            tree.nodes[successor]["label"] = str(framework.collect_genome_terms(parent)).replace(" ", "")
            tree.nodes[successor]["genome"] = parent
        for r in range(num_rearrangements):
            # Draw a rearrangement
            rearrangement = draw_rearrangement(model)
            # Apply the rearrangement
            genome = tree.nodes[successor]["genome"] * rearrangement * framework.symmetry_element()
            selected_genome = cgt.simulations.draw_genome(framework.collect_genome_terms(genome))
            selected_genome = [int(x) for x in selected_genome[:-1].strip('][').split(',')]
            selected_genome = framework.cycles(selected_genome) * framework.symmetry_element()
            # Update the genome at the successor node
            tree.nodes[successor]["genome"] = selected_genome
            tree.nodes[successor]["label"] = str(framework.collect_genome_terms(selected_genome)).replace(" ", "")
    
    return tree

def grow_tree(framework, model, num_taxa):
    model_element = model.s_element()
    # Start with the identity genome
    genome = framework.genome_group().identity()
    print(framework.draw_instance(genome, shortened=True))
    # Instantiate the tree
    G = nx.DiGraph()
    G.add_node(0, label=str(genome), genome=genome) # initialize root
    for t in range(0, num_taxa+1, 2):
        # Select a genome from the tree
        node = np.random.choice([node for node in G.nodes() if G.out_degree(node)==0])
        parent = G.nodes[node]['genome']

        # Draw a branch length From a multinomial distribution
        # Draw the number of events from a poisson distribution
        # Choose the rearangements from the model and apply them
        
    # Grow two new leaves
        for t_1 in range(2):
            # Draw a rearrangement
            rearrangement = draw_rearrangement(model)
            # Apply the rearrangement
            genome = parent * rearrangement
            # Draw a branch length
            branch_length = draw_branch_length()
            # Add the new genome to the tree
            G.add_node(t+t_1+1, label=str(genome), genome=genome)
            G.add_edge(node, t+t_1+1, weight=round(branch_length,1))

    # Return the tree
    return G

def draw_tree(tree):
    # Draw the tree
    plt.figure(figsize=(9,5))
    G = tree
    pos = graphviz_layout(G, prog='dot')#, args='-Granksep=0.1 -Gnodesep=0.1')
    labels = nx.get_node_attributes(G, 'label')
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw(G, pos=pos, node_size=0)
    nx.draw_networkx_labels(G, pos, labels=labels, bbox={"facecolor":"white", "edgecolor":"white"}, font_size=10)
    nx.draw_networkx_edge_labels(G, pos=pos, rotate=False, edge_labels=edge_labels, font_size=10)
    plt.show()


def draw_branch_length():
    return np.random.exponential(1)


def draw_rearrangement(model):
    model_element = model.s_element()
    rearrangements, probs = tuple(zip(*[list(x)[0] for x in (model_element).terms()]))
    return np.random.choice(rearrangements, p=probs)

def draw_genome(formal_sum):
    probs, genomes = zip(*list(formal_sum))
    return np.random.choice(genomes, p=probs)

def set_node(node, framework, genome):
    node["genome"] = genome
    node["label"] = str(framework.canonical_instance(genome))
