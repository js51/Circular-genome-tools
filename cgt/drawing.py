from sage.all_cmdline import *

def draw_genome(n=None, permutation=None, orientation_list=None, save_as_pgf=False, show=True):
	"""Produce a drawing of a circular genome with n regions"""
	# Imports
	import matplotlib
	import matplotlib.pyplot as plt
	if save_as_pgf:
		matplotlib.use("pgf")
		matplotlib.rcParams.update({
			"pgf.texsystem": "pdflatex",
			'font.family': 'serif',
			'text.usetex': True,
			'pgf.rcfonts': False
		})
	import numpy as np
	###
	
	region_labels = None
	orientation_list = None
	if permutation:
		n = len(permutation)
		region_labels = [None for _ in range(1, n+1)]
		orientation_list = [False for _ in range(1,n+1)]
		for i, r in enumerate(permutation):
			region_labels[abs(r)-1] = i+1
			orientation_list[abs(r)-1] = r < 0
		region_labels = list(reversed(region_labels))
		orientation_list = list(reversed(orientation_list))
	else:
		pass #region_labels = [k for k in range(1,n+1)]
		
	# Stuff for the pie chart
	segment_sizes = [100/n for _ in range(n)]
	segment_colors = ['white' for _ in range(n)]
	 
	fig1, ax1 = plt.subplots(figsize=(2,2))
	wedges, texts = ax1.pie(segment_sizes, labeldistance=1.2, colors=segment_colors, labels=region_labels, radius=1, startangle=90, wedgeprops={"edgecolor":"0", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
				
	if orientation_list:
		for i, p in enumerate(wedges):
			ang = (p.theta2 - p.theta1)/2. + p.theta1
			y = np.sin(np.deg2rad(ang))
			x = np.cos(np.deg2rad(ang))
			arrow_length = 0.00001 # We only want the arrow head
			angle = np.deg2rad(ang) + (np.pi/2 if orientation_list[i] else 3*np.pi/2)
			if orientation_list:
				ax1.arrow(x, y, arrow_length*np.cos(angle), arrow_length*np.sin(angle),overhang = 0.3, head_width=0.1, head_length=0.1, head_starts_at_zero=True, length_includes_head=True, linewidth=1, facecolor='black', fill=True)

	# Draw inner circle to hide most of the segment edges
	centre_circle = plt.Circle((0,0),0.85,fc='white')
	fig = plt.gcf()
	fig.gca().add_artist(centre_circle)
	ax1.axis('equal') # makes sure it's a circle
	plt.tight_layout()


	import os
	if not os.path.exists('output'):
		os.makedirs('output')
	plt.savefig('output/genome.png', facecolor='white', transparent=True)
	plt.savefig('output/genome' + ('.pgf' if save_as_pgf else '.pdf'), facecolor='white', transparent=True)
	if show and (not save_as_pgf):
		plt.show()