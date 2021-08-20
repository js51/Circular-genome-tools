"""
Implements a number of distance measures for genomes under the position paradigm.
"""

from cgt.enums import ALGEBRA
import numpy as np
import networkx as nx
from sage.all import ComplexDoubleField, UniversalCyclotomicField, matrix, Matrix, real, exp, round, CC
from scipy.optimize import minimize_scalar

def mles(framework, model, genome_instances):
    """Return maximum likelihood estimates for a set of genome instances under the given model and framework"""
    mles = {}
    for instance in genome_instances:
        mles[instance] = mle(framework, model, instance)
    return mles

def mle(framework, model, genome_instance):
    """Return maximum likelihood estimates for a genome instance under the given model and framework"""
    return maximise(framework, likelihood_function(framework, model, genome_instance, attempt_exact=False))

def maximise(framework, L, max_time=100):
    """Return the time that maximises likelihood function L, using additional information from the framework"""
    limit = 1/framework.num_genomes()
    t_max = minimize_scalar(lambda t: -1*L(t), 10, bounds=[(0, max_time)])['x'][0]
    mle = t_max if L(t_max)>limit else np.nan
    return mle

def _projection_operators(mat, eigs):
    """Return projection operators for given matrix and its eigenvalues"""
    dim = mat.nrows()
    projections = [matrix.identity(dim) for _ in eigs]
    for e1, eig1 in enumerate(eigs):
        for eig2 in eigs:
            if eig1 != eig2:
                projections[e1] *= (mat-(eig2*matrix.identity(dim)))*(1/(eig1-eig2))
    return projections

def _irreps_of_zs(framework, model, attempt_exact=False):
    """Return a set of matrices---images of zs under each irrep"""
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    z = framework.symmetry_element()
    s = model.s_element(in_algebra=ALGEBRA.genome)
    irreps_of_z, irreps_of_s = framework.irreps(z), framework.irreps(s)
    irreps_of_zs = (matrix(UCF, irrep_z*irrep_s) for irrep_z, irrep_s in zip(irreps_of_z, irreps_of_s))
    irreps_of_zs = [Matrix(UCF if attempt_exact else CDF, irrep_zs) for irrep_zs in irreps_of_zs]
    for irrep in irreps_of_zs:
	    irrep.set_immutable()
    return irreps_of_zs

def _eigenvalues(mat, round_to=8, make_real=True, inc_repeated=False, attempt_exact=False, use_numpy=True, bin_eigs=False, tol=10^(-8)):
    """Return all the eigenvalues for a given matrix mat"""
    col = list if inc_repeated else set
    if use_numpy:
        new_mat = np.array(matrix(CC, mat))
        all_eigs = np.linalg.eigvals(new_mat)
    else:
        all_eigs = (eig for eig in mat.eigenvalues())
    if attempt_exact: 
        return sorted(col(all_eigs))
    elif bin_eigs:
        return _bin(sorted(all_eigs), return_bin_size=False, tol=tol)
    else:
        return sorted(col(round(real(eig) if make_real else eig, round_to) for eig in all_eigs))

def _bin(eigenvalues, tol=10^(-8), return_bin_size=False):
    binned_eigenvals = []
    num_eigs = len(eigenvalues) # ??
    e = 0
    eigs = [real(eig) for eig in eigenvalues]
    while e < num_eigs:
        bin=[eigs[e]] # this is the eigenvalue
        while len(bin) < num_eigs-e and abs(eigs[e]-eigs[e+len(bin)]) < tol:
            bin.append(eigs[e+len(bin)])
        average = real(sum(bin))/len(bin)
        binned_eigenvals.append((average, len(bin)) if return_bin_size else average)
        e=e+len(bin)
    return binned_eigenvals

def _eigenvectors(mat, tol=10^(-8)):
    eigen_tuples = sorted(mat.eigenvectors_right())
    binned_eigenvals = _bin([et[0] for et in eigen_tuples], tol=tol, return_bin_size=True)
    # Orthogonalise the eigenvectors                      
    q=0
    eigenvectors = []
    for v in range(len(binned_eigenvals)):
        c=binned_eigenvals[v][1] 
        vec, _ = matrix([eigen_tuples[q+l][1][0].list() for l in range(c)]).gram_schmidt(orthonormal=True)
        eigenvectors.append(vec)
        q=q+c
    return binned_eigenvals, eigenvectors

def _partial_traces_for_genome(framework, instance, irreps, irreps_of_zs, projections, eig_lists):
    """Return dictionary of partial traces, indexed first by irrep index and then by eigenvalaue"""
    irreps_of_z = [irrep(framework.symmetry_element()) for irrep in irreps]
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    traces = {
        r: {
            eigenvalue: {} for eigenvalue in eig_lists[r]
        } for r in range(len(irreps_of_zs))
    }
    for r, irrep in enumerate(irreps): # Iterate over irreducible representations
        sigd = irreps_of_z[r] * Matrix(CDF, matrix(UCF, irrep(framework.cycles(instance.inverse()))))
        for e, eigenvalue in enumerate(eig_lists[r]):
            traces[r][eigenvalue] = round(real((sigd*projections[r][e]).trace()), 6)
    return traces

def _partial_traces_for_genome_using_eigenvectors(framework, instance, irreps, irreps_of_zs):
    """Return dictionary of partial traces, indexed first by irrep index and then by eigenvalaue"""
    irreps_of_z = [irrep(framework.symmetry_element()) for irrep in irreps]
    eigenvectors_list = [_eigenvectors(irrep_zs) for irrep_zs in irreps_of_zs]
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    traces = {
        r: {
            eigenvalue[0]: {} for eigenvalue in eigenvectors_list[r][0]
        } for r in range(len(irreps_of_zs))
    }
    for r, irrep in enumerate(irreps): # Iterate over irreducible representations
        sigd = irreps_of_z[r] * Matrix(CDF, matrix(UCF, irrep(framework.cycles(instance.inverse()))))
        eigenvalues = eigenvectors_list[r][0]
        eigenvector_list = eigenvectors_list[r][1]
        for e, eigenvalue in enumerate(eigenvalues):
            eigenvectors = eigenvector_list[e]
            traces[r][eigenvalue[0]]=round(real((sum([(eigenvectors[m,:].H)*eigenvectors[m,:] for m in range(eigenvalue[1])])*sigd).trace()),6)
    return traces


def likelihood_function(framework, model, genome, attempt_exact=False, use_projections=True):
    """Return the likelihood function for a given genome"""
    instance = genome
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    G, Z = framework.genome_group(), framework.symmetry_group()
    irreps = framework.irreps()
    irreps_of_zs = _irreps_of_zs(framework, model, attempt_exact=attempt_exact)
    if use_projections:
        eig_lists = [_eigenvalues(irrep_zs, round_to=8, make_real=True, inc_repeated=False, attempt_exact=attempt_exact) for irrep_zs in irreps_of_zs]
        projections = [_projection_operators(*vals) for vals in zip(irreps_of_zs, eig_lists)]
        traces = _partial_traces_for_genome(framework, instance, irreps, irreps_of_zs, projections, eig_lists)
    else:
        eig_lists = [[x[0] for x in _eigenvectors(irrep_zs)[0]] for irrep_zs in irreps_of_zs]
        traces = _partial_traces_for_genome_using_eigenvectors(framework, instance, irreps, irreps_of_zs)
    dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]
    def likelihood(t):
        ans = 0
        for r, dim in enumerate(dims):
            ans += Z.order()*(exp(-t)/G.order())*dim*sum(exp(eigenvalue*CDF(t)) * traces[r][eigenvalue] for eigenvalue in eig_lists[r])
        return real(ans)
    return likelihood

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