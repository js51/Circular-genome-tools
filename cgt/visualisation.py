"""What I want to be able to write:
ppf = PositionParadigmFramework(2)
genome = ppf.random_genome()

gt.visualise_genome(ppf, genome, type=VISUALISE_AS.interactive_graph)
gt.visualise_genome(ppf, genome, type=VISUALISE_AS.matplotlib_figure)

ppf.random_genome() # Returnes a 'genome' object

gt.visualisation.interactive(genome)

class genome:
    def __init__(self, framework, instance) # Problem is what is an instance???

maybe the framework should create these?

class genome:
    def __init__(self, framework, instance) # Only the framework creates these, so it can just convert to the group element first
    # example: return genome(self, instance_group_element)

    def as_formal_sum():
        return self.framework.formal_sum(self.instance) # or something?

Then we can go:

genome.as_matrix()
genome.as_formal_sum()
genome.random_instance()
genome.canonical_instance()

All of these use the embedded framework?
what would as_matrix do for a genome? Lots of repeated code because of genome/instance difference, but perhaps not avoidable?
or I suppose the genome methods are repeated instance ones, instance ones appearing in PPF and the genome ones appearing for the genome itself??
LEt's test some ideas!
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from sage.all_cmdline import *
from .position_paradigm import *
from .enums import *
import warnings

def draw_genome_instance(framework, instance, show=False):
    """Produce a drawing of a circular genome with n regions"""
    if not framework.oriented:
        raise NotImplementedError("not implemented for linear genomes.")
    try:
        instance = framework(instance)
    except:
        TypeError(f"{instance} does not belong in {str(framework)}")

    matplotlib.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "pgf.texsystem": "pdflatex",
        'pgf.rcfonts': False,
        "font.size": 18,
    })
    plt.style.use('default')

    n = framework.n
    region_labels     = [None for _ in range(1, n+1)]
    orientation_list  = [False for _ in range(1,n+1)]

    for i, r in reversed(list(enumerate(framework.one_row(instance, as_list=True)))):
        region_labels[abs(r)-1] = i+1
        orientation_list[abs(r)-1] = r < 0
    
    segment_sizes = [100/n for _ in range(n)]
    segment_colors = ['white' for _ in range(n)]
    
    fig1, ax1 = plt.subplots(figsize=(2,2))
    wedges, texts = ax1.pie(segment_sizes, labeldistance=1.2, colors=segment_colors, 
                            labels=region_labels, radius=1, startangle=90, 
                            wedgeprops={"edgecolor":"0", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
				
    if framework.oriented and orientation_list:
        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1)/2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))
            arrow_length = 0.00001 # We only want the arrow head
            angle = np.deg2rad(ang) + (np.pi/2 if orientation_list[i] else 3*np.pi/2)
            if orientation_list:
            	ax1.arrow(x, y, arrow_length*np.cos(angle), arrow_length*np.sin(angle), 
                          overhang = 0.3, head_width=0.12, head_length=0.12, head_starts_at_zero=True, 
                          length_includes_head=True, linewidth=1, facecolor='black', fill=True)

	# Draw inner circle to hide most of the segment edges
    centre_circle = plt.Circle((0,0),0.85,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)
    ax1.axis('equal') # makes sure it's a circle
    plt.tight_layout()

    if not os.path.exists('_output'):
        os.makedirs('_output')
    plt.savefig(f'_output/{str(instance)}.png', transparent=True)
    plt.savefig(f'_output/{str(instance)}.pdf', transparent=True)
    if show: plt.show()