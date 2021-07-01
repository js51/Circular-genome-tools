"""
"""
from cgt.enums import ALGEBRA
import numpy as np
import networkx as nx
from sage.all import ComplexDoubleField, UniversalCyclotomicField, matrix, Matrix, real, exp, round
from scipy.optimize import minimize
from warnings import warn
import matplotlib.pyplot as plt

def mles(framework, model, genome_reps, plots=False):
    """
    Return dictionary of MLEs (if they exist) for time elapsed rearranging ref->genome for each genome
    """
    warn("This... might take a while!")
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    G = framework.genome_group()
    Z = framework.symmetry_group()
    z = framework.symmetry_element()
    s = model.s_element(in_algebra=ALGEBRA.genome)
    irreps_of_z, irreps_of_s = framework.irreps(z), framework.irreps(s)
    irreps_of_zs = [ Matrix(CDF, matrix(UCF, irrep_z*irrep_s)) for irrep_z, irrep_s in zip(irreps_of_z, irreps_of_s) ]
    for irrep in irreps_of_zs:
	    irrep.set_immutable()
    dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]
    eig_lists   = [_eigenvalues(irrep_zs, round_to=9, make_real=True, inc_repeated=False) for irrep_zs in irreps_of_zs]
    projections = [_projection_operators(*vals) for vals in zip(irreps_of_zs, eig_lists)]
    traces = {
        r: {
            genome_rep: { 
                eigenvalue: {} for eigenvalue in eig_lists[r]
            } for genome_rep in genome_reps
        } for r in range(len(irreps_of_zs))
    }
    inv_sigmas = { instance : framework.cycles(instance.inverse()) for instance in genome_reps } # pre-compute for efficiency
    irreps = framework.irreps()
    for r, irrep in enumerate(irreps): # Iterate over irreducible representations
        for instance in genome_reps:
            sigd = Matrix(CDF, matrix(UCF, irrep(inv_sigmas[instance]))) * irreps_of_z[r]
            for e, eigenvalue in enumerate(eig_lists[r]):
                traces[r][instance][eigenvalue] = round(real((sigd*projections[r][e]).trace()), 6)
    def likelihood(t, instance):
        ans = 0
        for r in range(len(irreps_of_zs)):
            ans += Z.order()*(exp(-t)/G.order())*dims[r]*sum(exp(eigenvalue*CDF(t)) * traces[r][instance][eigenvalue] for eigenvalue in eig_lists[r])
        return real(ans)
    mle_dict = {}
    bound = (0, 35)
    for instance in genome_reps:
        def L(time):
            return -likelihood(time, instance)
        mle = minimize(L, 0.1, method="TNC", bounds=[bound])
        mle_dict[instance] = mle.x[0]
        if plots:
            plt.figure()
            times = np.arange(bound[0], bound[1], 0.5)
            likelihood_functional_values = [likelihood(time, instance) for time in times]
            plt.plot(times, likelihood_functional_values)
            plt.savefig(f'./cgt_{framework.one_row(instance)}')

    return {framework.one_row(key) : value for key, value in mle_dict.items()}

def _projection_operators(mat, eigs):
        """Return projection operators for given matrix and its eigenvalues"""
        dim = mat.nrows()
        projections = [matrix.identity(dim) for _ in eigs]
        for e1, eig1 in enumerate(eigs):
            for eig2 in eigs:
                if eig1 != eig2:
                    projections[e1] *= (mat-(eig2*np.eye(dim)))*(1/(eig1-eig2))
        return projections

def _eigenvalues(mat, round_to=9, make_real=True, inc_repeated=False):
    col = list if inc_repeated else set
    collection = col(round(real(eig) if make_real else eig, round_to) for eig in mat.eigenvalues())
    return sorted(collection)

def likelihood_function():
    """Return the likelihood function for a given genome"""
    pass

def min_distance(framework, model, weighted=False):
    """Return dictionary of minimum distances ref->genome, using inverse of model probabilities as weights (or not)"""
    matrix = model.reg_rep_of_zs().toarray()
    if weighted:
        matrix = np.reciprocal(matrix, where=matrix!=0)
    graph = nx.Graph(matrix)
    genome_reps = framework.genomes().keys()
    return { 
        rep : nx.shortest_path_length(graph, source=0, target=r, weight='weighted' if weighted else None) 
        for r, rep in enumerate(genome_reps) 
    }

def MFPT(framework, model, scale_by=1):
    """Return the mean time elapsed rearranging ref->genome where the target is an absorbing state"""
    genomes = framework.genomes()
    P = model.reg_rep_of_zs().toarray()
    P_0 = P.copy()
    for i in range(len(P)):
        P_0[0,i] = 0 
    Z = np.linalg.inv(np.eye(len(P)) - P_0)
    def MFTP_dists():
        return (Z@P_0@Z)[:,0]
    reps = list(genomes.keys())
    MFTP_distances = list(MFTP_dists())
    MFTP_distances = { rep : scale_by * MFTP_distances[r] for r, rep in enumerate(reps) }
    return MFTP_distances