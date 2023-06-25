"""
Implements a number of distance measures for genomes under the position paradigm.

The primary function of this module is the distance_matrix function, which returns a distance matrix for a given set of genomes and distance measure.

The distance measures implemented are:
    - min: the minimum number of rearrangements required to transform one genome into another
    - min_weighted: the minimum number of rearrangements required to transform one genome into another, weighted by the inverse of the probability of the rearrangement
    - MFPT: the mean first passage time from the identity to a given genome, where the target is an absorbing state
    - MLE: the maximum likelihood estimate of the time elapsed between the identity and a given genome

Individual likelihood functions for the time elapsed between the identity and a given genome can also be obtained
"""

from cgt.enums import ALGEBRA, DISTANCE
import numpy as np
import networkx as nx
from sage.all import ComplexDoubleField, UniversalCyclotomicField, matrix, Matrix, real, exp, round, CC, log
from scipy.optimize import minimize_scalar

def mles(framework, model, genome_instances=None, verbose=False):
    """
    Returns a dictionary of maximum likelihood estimates for each genome instance under the given model and framework.

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        genome_instances (list): a list of genome instances to compute MLEs for. If None, all genomes in the framework are used.
        verbose (bool): whether to print progress
    
    Returns:
        dict: a dictionary of maximum likelihood estimates, indexed by genome instance
    """
    mles = {}
    if genome_instances is None:
        genome_instances = [framework.canonical_instance(g) for g in framework.genomes()]
    for instance in genome_instances:
        if verbose: print(f"Computing MLE for {instance}")
        mles[instance] = mle(framework, model, instance)
        return mles

def mle(framework, model, genome_instance):
    """
    Returns the maximum likelihood estimate for a given genome instance under the given model and framework.

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        genome_instance (group element): the genome instance to compute the MLE for

    Returns:
        float: the maximum likelihood estimate
    """
    return maximise(framework, likelihood_function(framework, model, genome_instance))

def maximise(framework, L, max_time=100):
    """
    Return the time that maximises likelihood function L, using additional information from the framework
    
    Args:
        framework (PositionParadigmFramework): the framework
        L (function): the likelihood function to maximise
        max_time (float): the maximum time to consider (default: 100)
    
    Returns:
        float: the time that maximises the likelihood function
    """
    limit = 1/framework.num_genomes()
    t_max = minimize_scalar(lambda t: -1*L(t), method='bounded', bounds=(0, max_time))['x']
    mle = t_max if L(t_max)>limit else np.nan
    return mle

def _projection_operators(mat, eigs):
    """
    Return projection operators for given matrix and its eigenvalues
    
    Args:
        mat (matrix): a representation of zs
        eigs (list): a list of eigenvalues of mat

    Returns:
        list: a list of projection operators
    """
    dim = mat.nrows()
    projections = [matrix.identity(dim) for _ in eigs]
    for e1, eig1 in enumerate(eigs):
        for eig2 in eigs:
            if eig1 != eig2:
                projections[e1] *= (mat-(eig2*matrix.identity(dim)))*(1/(eig1-eig2))
    return projections

def _irreps_of_zs(framework, model, force_recompute=False):
    """
    Return a set of matrices---images of zs under each irrep
    
    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model 
        force_recompute (bool): whether to force recomputation and invalidate cache (default: False)

    Returns:
        list: a list of irredicuble representations of the group applied to zs
    """
    key = "irreps_of_zs"
    if key in model.data_bundle and not force_recompute:
        irreps_of_zs =  model.data_bundle[key]
    else:
        CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
        z = framework.symmetry_element()
        s = model.s_element(in_algebra=ALGEBRA.genome)
        irreps_of_z, irreps_of_s = framework.irreps(z), framework.irreps(s)
        irreps_of_zs = (matrix(UCF, irrep_z*irrep_s) for irrep_z, irrep_s in zip(irreps_of_z, irreps_of_s))
        irreps_of_zs = [Matrix(CDF, irrep_zs) for irrep_zs in irreps_of_zs]
        for irrep in irreps_of_zs:
            irrep.set_immutable()
        model.data_bundle[key] = irreps_of_zs
    return irreps_of_zs

def _irreps_of_z(framework, model, force_recompute=False):
    """
    Return a set of matrices---images of z (the symmetry element) under each irrep
    
    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        force_recompute (bool): whether to force recomputation and invalidate cache (default: False)

    Returns:
        list: a list of irredicuble representations of the group applied to z
    """
    key = "irreps_of_z"
    if key in model.data_bundle and not force_recompute:
        irreps_of_z =  model.data_bundle[key]
    else:
        CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
        z = framework.symmetry_element()
        irreps_of_z = framework.irreps(z)
        irreps_of_z = (matrix(UCF, irrep_z) for irrep_z in irreps_of_z)
        irreps_of_z = [Matrix(CDF, irrep_z) for irrep_z in irreps_of_z]
        for irrep in irreps_of_z:
            irrep.set_immutable()
        model.data_bundle[key] = irreps_of_z
    return irreps_of_z


def _eigenvalues(mat, round_to=7, make_real=True, inc_repeated=False, use_numpy=True, bin_eigs=False, tol=10**(-8)):
    """Return all the eigenvalues for a given matrix mat"""
    col = list if inc_repeated else set
    if use_numpy:
        new_mat = np.array(matrix(CC, mat))
        all_eigs = np.linalg.eigvals(new_mat)
    else:
        all_eigs = (eig for eig in mat.eigenvalues())
    if bin_eigs:
        return _bin(sorted(all_eigs), return_bin_size=False, tol=tol)
    else:
        return sorted(col(round(real(eig) if make_real else eig, round_to) for eig in all_eigs))

def _bin(eigenvalues, tol=10**(-8), return_bin_size=False):
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

def _eigenvectors(mat, tol=10**(-8)):
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

def _partial_traces_for_genome(framework, instance, irreps, irreps_of_zs, projections, eig_lists, irreps_of_z=None):
    """Return dictionary of partial traces, indexed first by irrep index and then by eigenvalaue"""
    if irreps_of_z is None:
        irreps_of_z = [irrep(framework.symmetry_element()) for irrep in irreps]
    traces = {
        r: {
            eigenvalue: {} for eigenvalue in eig_lists[r]
        } for r in range(len(irreps_of_zs))
    }
    for r, irrep in enumerate(irreps): # Iterate over irreducible representations
        sigd = irreps_of_z[r] * irrep(instance.inverse())
        for e, eigenvalue in enumerate(eig_lists[r]):
            traces[r][eigenvalue] = real((sigd*projections[r][e]).trace())
    return traces

def likelihood_function(framework, model, genome):
    """Return the likelihood function for a given genome"""
    instance = genome
    CDF, UCF = ComplexDoubleField(), UniversalCyclotomicField()
    G, Z = framework.genome_group(), framework.symmetry_group()
    irreps = framework.irreps()
    irreps_of_zs = _irreps_of_zs(framework, model)
    irreps_of_z = _irreps_of_z(framework, model)
    if "eig_lists" in model.data_bundle:
        eig_lists = model.data_bundle["eig_lists"]
    else:
        eig_lists = [_eigenvalues(irrep_zs) for irrep_zs in irreps_of_zs]
        model.data_bundle["eig_lists"] = eig_lists
    if "projections" in model.data_bundle:
        projections = model.data_bundle["projections"]
    else:
        projections = [_projection_operators(*vals) for vals in zip(irreps_of_zs, eig_lists)]
        model.data_bundle["projections"] = projections
    traces = _partial_traces_for_genome(framework, instance, irreps, irreps_of_zs, projections, eig_lists, irreps_of_z)
    dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]
    def likelihood(t):
        ans = 0
        for r, dim in enumerate(dims):
            ans += Z.order()*(exp(-t)/G.order())*dim*sum(exp(eigenvalue*CDF(t)) * traces[r][eigenvalue] for eigenvalue in eig_lists[r])
        return real(ans)
    return likelihood

def min_distance(framework, model, genome_reps=None, weighted=False):
    """Return dictionary of minimum distances ref->genome, using inverse of model probabilities as weights (or not)"""
    matrix = model.reg_rep_of_zs().toarray()
    if weighted:
        matrix[matrix != 0] = 1 / matrix[matrix != 0]
    graph = nx.Graph(matrix)
    genomes = list(framework.genomes().keys())
    if genome_reps is None:
        genome_reps = genomes
    return { 
        rep : nx.shortest_path_length(graph, source=0, target=genomes.index(rep), weight='weight' if weighted else None) 
        for r, rep in enumerate(genome_reps) 
    }

def MFPT(framework, model, genome_reps=None, scale_by=1):
    """Return the mean time elapsed rearranging ref->genome where the target is an absorbing state"""
    genomes = framework.genomes()
    reps = list(genomes.keys())
    P = model.reg_rep_of_zs().toarray()
    P_0 = P.copy()
    for i in range(len(P)):
        P_0[0,i] = 0 
    Z = np.linalg.inv(np.eye(len(P)) - P_0)
    def MFTP_dists():
        return (Z@P_0@Z)[:,0]
    MFTP_distances = list(MFTP_dists())
    MFTP_distances = { rep : scale_by * MFTP_distances[r] for r, rep in enumerate(reps) }
    if genome_reps is not None:
        MFTP_distances = { rep : MFTP_distances[framework.canonical_instance(rep)] for rep in genome_reps }
    return MFTP_distances

def dict_to_distance_matrix(distances, framework, genomes=None):
    """If need to convert to pairwise distances, supply a list of genomes."""
    if genomes is not None:
        D = np.zeros((len(genomes), len(genomes)))
        for i in range(len(genomes)):
            canonical_i = framework.canonical_instance(framework.random_instance(genomes[i]))
            distances_copy = { canonical_i * k : v for k, v in distances.items() }
            for j in range(len(genomes)):
                canonical_j = framework.canonical_instance(framework.random_instance(genomes[j]))
                D[i,j] = distances_copy[canonical_j]
    else:
        D = np.zeros((len(distances), len(distances)))
        for (i, j), distance in distances.items():
            D[i,j] = distance
    return D

def distance_matrix(framework, model, genomes, distance):
    """
    Compute a distance matrix for a given set of genomes and distance measure.

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        genomes (list): a list of genomes to compute distances for
        distance (DISTANCE): the distance measure to use
    """
    instances = [framework.canonical_instance(g) for g in genomes]
    # Get the genomes g that we need dist(id -> g) for:
    need_distances = {}
    for i, _ in enumerate(genomes):
        for j, _ in enumerate(genomes):
            if i < j:
                canonical_i, canonical_j = instances[i], instances[j]
                need_distances[(i, j)] = canonical_i.inverse() * canonical_j
    pairs, need_distances = map(list, zip(*need_distances.items()))

    # Compute the distances:
    if distance == DISTANCE.MFPT:
        distances = MFPT(framework, model, genome_reps=need_distances)
    elif distance == DISTANCE.min:
        distances = min_distance(framework, model, genome_reps=need_distances)
    elif distance == DISTANCE.min_weighted:
        distances = min_distance(framework, model, genome_reps=need_distances, weighted=True)
    elif distance == DISTANCE.MLE:
        distances = mles(framework, model, genome_instances=need_distances, verbose=True)

    # Construct the distance matrix
    D = np.zeros((len(genomes), len(genomes)))
    for p, (i, j) in enumerate(pairs):
        D[i,j] = distances[need_distances[p]]
    
    return D