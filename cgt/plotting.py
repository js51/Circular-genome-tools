"""
This module contains functions for producing a selection of pre-styled matplotlib plots, with settings specific to the context of circular genomes with minimal other options available. For fully custom plots, use the matplotlib library directly along with the cgt.distances module.
"""

import cgt
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np


def latex_figure_setup():
    """
    Sets up the matplotlib environment to use LaTeX for all text.
    """
    nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        "axes.labelsize": 11,
        "font.size": 11,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    }
    plt.style.use("default")
    mpl.rcParams.update(nice_fonts)
    monochrome = (
        cycler("color", ["k"])
        * cycler("marker", ["", "."])
        * cycler("linestyle", ["-", "--", ":"])
    )
    plt.rc("axes", prop_cycle=monochrome)


def x_axis_limit_for_framework(framework):
    n = framework.n
    from_testing = {1: 10, 2: 10, 3: 10, 4: 14, 5: 22, 6: 26, 7: 34, 8: 40, 9: 98}
    if n in from_testing:
        return from_testing[n]
    else:
        return from_testing[2 * max(from_testing.keys())]


def plot_likelihood_and_probability(framework, model, instance, width=6, height=4, tmax = None, dpi=300, use_eigenvectors=True):
    a = cgt.distances.prob_to_reach_in_steps_func(framework, model, instance, use_eigenvectors=use_eigenvectors)
    L = cgt.distances.likelihood_function(
        framework, model, instance, use_eigenvectors=use_eigenvectors
    )
    if tmax is None:
        max_x = x_axis_limit_for_framework(framework)
    else:
        max_x = tmax
    tvec_L = list(np.arange(0, max_x + 1, 0.2))
    tvec = list(np.arange(0, max_x + 1, 1))
    latex_figure_setup()
    plt.figure(figsize=(width, height), dpi=dpi)
    plt.scatter(tvec, [a(t) for t in tvec], color="black", s=3)
    plt.plot(
        tvec,
        [a(t) for t in tvec],
        color="black",
        linestyle="dotted",
        label="Probability",
        linewidth=0.9,
    )
    plt.plot(
        tvec_L,
        [L(t) for t in tvec_L],
        color="black",
        linestyle="solid",
        label="Likelihood",
    )
    plt.title(f"Genome: z{str(framework.canonical_instance(instance))}")
    n = framework.n
    if n < 9:
        x_tick_list = list(range(0, 9, 1)) + list(range(10, max_x + 1, 2))
    else: 
        x_tick_list = list(range(0, max_x + 1, 5))
    plt.xticks(x_tick_list)
    # Legend without border
    # bbox is (x, y, width, height)
    plt.legend(loc="best", frameon=False, bbox_to_anchor=(0.7, 0.0, 0.3, 1.0))
    # As the x axis label, show the MLE
    mle = cgt.distances.maximise(framework, L)
    plt.xlabel(
        f"MLE: {round(mle, 1) if not np.isnan(mle) else 'N/A'}, Min dist: {cgt.distances.first_nonzero_value(framework, a)}"
    )
    plt.show()
    return a, L
