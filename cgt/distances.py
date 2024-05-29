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

from cgt.enums import ALGEBRA, DISTANCE, DATA, SYMMETRY
from cgt.constants import VALUES_OF_N_WITH_SAVED_IRREPS_OF_Z
import numpy as np
import networkx as nx
from sage.all import matrix, real, exp, round
import scipy
from scipy.optimize import minimize_scalar
from cgt import pickle_manager
from warnings import warn


def mles(framework, model, genome_instances=None, verbose=False, show_work=False):
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
        genome_instances = [
            framework.canonical_instance(g) for g in framework.genomes()
        ]
    for instance in genome_instances:
        if verbose:
            print(f"Computing MLE for {instance}")
        mles[instance] = mle(framework, model, instance, show_work=show_work)
    return mles


def distance(framework, model, genome_instance, distance_measure):
    distance_func = get_distance_function(distance_measure)
    return distance_func(framework, model, [genome_instance])


def distance_between(framework, model, genome_instance_1, genome_instance_2, distance_measure):
    return distance(framework, model, genome_instance_1.inverse() * genome_instance_2, distance_measure)


def mle(framework, model, genome_instance, show_work=False):
    """
    Returns the maximum likelihood estimate for a given genome instance under the given model and framework.

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        genome_instance (group element): the genome instance to compute the MLE for

    Returns:
        float: the maximum likelihood estimate
    """
    L = likelihood_function(framework, model, genome_instance)
    max_t = maximise(framework, L)
    if show_work:
        return max_t, L
    else:
        return max_t


def maximise(framework, L, max_time=None):
    """
    Return the time that maximises likelihood function L, using additional information from the framework

    Args:
        framework (PositionParadigmFramework): the framework
        L (function): the likelihood function to maximise
        max_time (float): the maximum time to consider (default: 100)

    Returns:
        float: the time that maximises the likelihood function
    """
    if max_time is None: 
        max_time = framework.n * 25
    limit = 1 / framework.num_genomes()
    t_max = minimize_scalar(
        lambda t: -1 * L(t), method="bounded", bounds=(0, max_time)
    )["x"]
    mle = t_max if L(t_max) > limit * (1 + (1/3)*0.01) else np.nan

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
                projections[e1] *= (mat - (eig2 * matrix.identity(dim))) * (
                    1 / (eig1 - eig2)
                )
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
    key = DATA.irreps_zs
    if key in model.data_bundle and not force_recompute:
        irreps_of_zs = model.data_bundle[key]
    else:
        z = framework.symmetry_element()
        s = model.s_element(in_algebra=ALGEBRA.genome)
        irreps_of_z = _irreps_of_z(framework, model)
        irreps_of_s = framework.irreps(s)
        irreps_of_zs = [
            irrep_z * irrep_s for irrep_z, irrep_s in zip(irreps_of_z, irreps_of_s)
        ]
        model.data_bundle[key] = irreps_of_zs
    return irreps_of_zs


def _irreps_of_z(framework, model=None, force_recompute=False):
    """
    Return a set of matrices---images of z (the symmetry element) under each irrep

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        force_recompute (bool): whether to force recomputation and invalidate cache (default: False)

    Returns:
        list: a list of irredicuble representations of the group applied to z
    """
    key = DATA.irreps_z
    if model is not None and key in model.data_bundle and not force_recompute:
        return model.data_bundle[key]
    if framework.n in VALUES_OF_N_WITH_SAVED_IRREPS_OF_Z and framework.symmetry is SYMMETRY.circular and not force_recompute:
        try:
            irreps_of_z = pickle_manager.retrieve_irrep_of_z(n=framework.n)
            if model is not None: 
                model.data_bundle[key] = irreps_of_z
                model.data_bundle[DATA.irreps_z_np] = [irrep.numpy() for irrep in irreps_of_z] # TODO: This is a hack
            return irreps_of_z
        except:
            pass
    z = framework.symmetry_element()
    irreps_of_z = framework.irreps(z)
    if model is not None: 
        model.data_bundle[key] = irreps_of_z
        model.data_bundle[DATA.irreps_z_np] = [irrep.numpy() for irrep in irreps_of_z] # TODO: This is a hack
    return irreps_of_z


def _eigenvalues(mat, round_to=7, make_real=True, inc_repeated=False):
    """Return all the eigenvalues for a given matrix mat"""
    col = list if inc_repeated else set
    all_eigs = np.linalg.eigvals(mat)
    return sorted(
        col(round(real(eig) if make_real else eig, round_to) for eig in all_eigs)
    )


def _eigen_data(framework, model, irreps_of_zs, round_vecs_to=10, round_vals_to=7):
    """
    Return an eigenvector matrix for each irrep of zs. Also return eigenvalue and eigenvector lists

    Args:
        irreps_of_zs (list): a list of irreps of zs
        round_vecs_to (int): the number of decimal places to round eigenvectors to (default: 10)
        round_vals_to (int): the number of decimal places to round eigenvalues to (default: 7)

    Returns:
        list: a list of eigenvector matrices
        list: a list of eigenvalue lists
        list: a list of eigenvector lists
    """
    eigen_data = {
        DATA.eigval_lists: [],
        DATA.eigval_sets: [],
        DATA.eigvec_lists: [],
        DATA.eigvec_mat_inv: [],
    }

    # Fix for weird issue with SVD convergence (n=9)
    if framework.n == 9 and round_vecs_to == 10:
        round_vecs_to = 12

    # If the data is already stored in the model, return it
    if DATA.eig_data in model.data_bundle:
        return model.data_bundle[DATA.eig_data]

    # Otherwise compute the eigen_data dictionary
    for irrep_of_zs in irreps_of_zs:
        # Prepare irrep matrices
        irrep_zs_np = irrep_of_zs.numpy()
        # Compute eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eig(irrep_zs_np)
        # Rounding
        eigenvectors = np.around(eigenvectors, round_vecs_to)  # a matrix
        #print(f"Condition number is: {np.linalg.cond(eigenvectors)}")
        eigenvalues = [round(val, round_vals_to) for val in eigenvalues]
        # Store data
        eigen_data[DATA.eigval_lists].append(eigenvalues)
        eigen_data[DATA.eigval_sets].append(set(eigenvalues))
        eigen_data[DATA.eigvec_lists].append(eigenvectors.T)
        inverse_of_eigenvector_mat = np.linalg.pinv(eigenvectors)
        #print(f"Elements of inverted eigenvector matrix range from ({np.amin(inverse_of_eigenvector_mat)} to {np.amax(inverse_of_eigenvector_mat)})")
        eigen_data[DATA.eigvec_mat_inv].append(inverse_of_eigenvector_mat)


    # Store the data in the model
    model.data_bundle[DATA.eig_data] = eigen_data
    return eigen_data


def probability_to_reach_in_steps(framework, model, sigma, k):
    alpha = prob_to_reach_in_steps_func(framework, model, sigma)
    return alpha(k)


def prob_to_reach_in_steps_func(framework, model, sigma, use_eigenvectors=True):
    G, Z = framework.genome_group(), framework.symmetry_group()
    instance = framework.cycles(sigma)
    irreps = framework.irreps()
    irreps_of_zs = _irreps_of_zs(framework, model)
    irreps_of_z = _irreps_of_z(framework, model)
    if DATA.eigval_sets_old in model.data_bundle:
        eig_lists = model.data_bundle[DATA.eigval_sets_old]
    else:
        eig_lists = [_eigenvalues(irrep_zs) for irrep_zs in irreps_of_zs]
        model.data_bundle[DATA.eigval_sets_old] = eig_lists
    if not use_eigenvectors:
        if DATA.projections in model.data_bundle:
            projections = model.data_bundle[DATA.projections]
        else:
            projections = [
                _projection_operators(*vals) for vals in zip(irreps_of_zs, eig_lists)
            ]
            model.data_bundle[DATA.projections] = projections
        traces = _partial_traces_for_genome(
            framework, model, instance, irreps, irreps_of_zs, projections, eig_lists, irreps_of_z
        )
    else:
        eig_lists, traces = _partial_traces_for_genome_using_eigenvectors(
            framework, model, instance, irreps, irreps_of_zs, irreps_of_z
        )
    dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]

    def prob_in_steps(k):
        alpha = 0
        for ptrace, dim, eig_list in zip(traces.values(), dims, eig_lists):
            term = 0
            for e, eig in enumerate(eig_list):
                term += (eig**k) * ptrace[eig]
            alpha += dim * term
        return (Z.order() / G.order()) * alpha

    return prob_in_steps


def discrete_MFPT(framework, model, genome_reps=None, verbose=False):
    if genome_reps is None:
        genomes = list(framework.genomes().keys())
        genome_reps = genomes
    dists = {}
    a_star = framework.symmetry_group().order() / framework.genome_group().order()
    identity = framework.one_row(framework.genome_group().identity())
    for instance in genome_reps:
        if instance == identity:
            dists[instance] = 0
        else:
            a = prob_to_reach_in_steps_func(framework, model, instance)
            summands = [0]
            a_values = [0]
            k = 1
            prod = 1
            while True:
                a_k = a(k)
                a_values.append(a_k)
                prod *= 1 - a_values[k - 1]
                summands.append(k * a_k * prod)
                if (abs(a_values[k] - a_values[k - 1]) < 10 ** (-3)) and (
                    abs(a_values[k] - a_star) < 10 ** (-3)
                ):
                    if verbose:
                        print(
                            f"Terms {a_values[k]} and {a_values[k-1]} are close enough to {a_star}"
                        )
                    break
                k += 1
            c = k
            if verbose:
                print(f"Converged at {k=}")
            adjustment_term = sum(
                j * a_star * (1 - a_star) ** j for j in range(0, c + 1)
            )
            remaining_sum = ((prod * a_star) / ((1 - a_star) ** (c + 1))) * (
                ((1 - a_star) / (a_star**2)) - adjustment_term
            )
            dists[instance] = sum(summands) + remaining_sum
    return dists


def _partial_traces_for_genome(
    framework,
    model,
    instance,
    irreps,
    irreps_of_zs,
    projections,
    eig_lists,
    irreps_of_z=None,
):
    """Return dictionary of partial traces, indexed first by irrep index and then by eigenvalaue"""
    instance_inverse = instance.inverse()
    if (
        irreps_of_z is None
    ):  # It almost never is. Can't used the cached version below because it's stored in the model (which we don't have)
        irreps_of_z = framework.irreps(framework.symmetry_element())
    traces = {
        r: {eigenvalue: {} for eigenvalue in eig_lists[r]}
        for r in range(len(irreps_of_zs))
    }
    for r, irrep in enumerate(irreps):  # Iterate over irreducible representations
        zero_irrep = not irreps_of_z[r] 
        if zero_irrep:  # matrix of zeros
            sigd = irreps_of_z[r]
        else:
            sigd = irreps_of_z[r] * irrep(instance_inverse)
        for e, eigenvalue in enumerate(eig_lists[r]):
            if zero_irrep:  # matrix of zeros
                traces[r][eigenvalue] = 0
            else:
                traces[r][eigenvalue] = real((sigd.numpy() @ projections[r][e]).trace())
    return traces


def _partial_traces_for_genome_using_eigenvectors(
    framework, model, instance, irreps, irreps_of_zs, irreps_of_z
):
    """Return dictionary of partial traces, indexed first by irrep index and then by eigenvalaue"""
    instance_inverse = instance.inverse()
    eigen_data = _eigen_data(framework, model, irreps_of_zs)
    if DATA.partial_traces in model.data_bundle and instance_inverse in model.data_bundle[DATA.partial_traces]:
        return eigen_data[DATA.eigval_sets], model.data_bundle[DATA.partial_traces][instance_inverse]
    traces = {}
    irreps_of_z_np = model.data_bundle[DATA.irreps_z_np]
    for r, irrep in enumerate(irreps):  # Iterate over irreducible representations
        # Faster convert to numpy
        traces[r] = {}
        irrep_instance = irrep(instance_inverse)
        irrep_instance_np = np.zeros(irrep_instance.dimensions())
        for (i, j), val in irrep_instance.items():
            irrep_instance_np[i,j] = val
        eigenvalue_list = eigen_data[DATA.eigval_lists][r]
        if not irreps_of_z[r]:  # matrix of zeros
            irrep_of_g_inverse_z_np = irreps_of_z[r]
            traces[r] = {eigenvalue: 0 for eigenvalue in eigenvalue_list}
        else:
            irrep_of_g_inverse_z_np = irreps_of_z_np[r] @ irrep_instance_np
            eigenvectors = eigen_data[DATA.eigvec_lists][r]
            P_inv = eigen_data[DATA.eigvec_mat_inv][r]
            mat = ((P_inv @ irrep_of_g_inverse_z_np) * eigenvectors).sum(-1)
            for v, eig in enumerate(eigenvalue_list):
                if eig in traces[r]:
                    traces[r][eig] += mat[v]
                else:
                    traces[r][eig] = mat[v]
        full_trace = irrep_of_g_inverse_z_np.trace()
        sum_of_partial_traces = sum(traces[r].values())
        if not sum_of_partial_traces - full_trace < 1e-04:
            print(f"Sum of partial traces ({round(sum_of_partial_traces, 5)}) is very different from the full trace ({round(full_trace, 5)}) for irrep {r}.")
    
    if DATA.partial_traces not in model.data_bundle:
        model.data_bundle[DATA.partial_traces] = {}
    model.data_bundle[DATA.partial_traces][instance_inverse] = traces

    return eigen_data[DATA.eigval_sets], traces


def _eig_lists(model, irreps_of_zs):
    if DATA.eigval_sets_old in model.data_bundle:
        eig_lists = model.data_bundle[DATA.eigval_sets_old]
    else:
        eig_lists = [_eigenvalues(irrep_zs) for irrep_zs in irreps_of_zs]
        model.data_bundle[DATA.eigval_sets_old] = eig_lists
    return eig_lists


def _projections(model, irreps_of_zs, eig_lists):
    if DATA.projections in model.data_bundle:
        projections = model.data_bundle[DATA.projections]
    else:
        projections = [
            _projection_operators(*vals) for vals in zip(irreps_of_zs, eig_lists)
        ]
        model.data_bundle[DATA.projections] = projections
    return projections


def likelihood_function(framework, model, genome, use_eigenvectors=True):
    """Return the likelihood function for a given genome"""
    instance = framework.cycles(genome)
    G, Z = framework.genome_group(), framework.symmetry_group()
    irreps = framework.irreps()
    irreps_of_zs = _irreps_of_zs(framework, model)
    irreps_of_z = _irreps_of_z(framework, model)
    if not use_eigenvectors:
        eig_lists = _eig_lists(model, irreps_of_zs)
        projections = _projections(model, irreps_of_zs, eig_lists)
        traces = _partial_traces_for_genome(
            framework,
            model,
            instance,
            irreps,
            irreps_of_zs,
            projections,
            eig_lists,
            irreps_of_z,
        )
    else:
        eig_lists, traces = _partial_traces_for_genome_using_eigenvectors(
            framework, model, instance, irreps, irreps_of_zs, irreps_of_z
        )
    dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]

    def likelihood(t):
        ans = 0
        for r, dim in enumerate(dims):
                ans += (
                    Z.order()
                    * (exp(-t) / G.order())
                    * dim
                    * sum(
                        exp(eigenvalue * t) * traces[r][eigenvalue]
                        for eigenvalue in eig_lists[r]
                    )
                )
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
        rep: nx.shortest_path_length(
            graph,
            source=0,
            target=genomes.index(rep),
            weight="weight" if weighted else None,
        )
        for r, rep in enumerate(genome_reps)
    }


def min_distance_using_irreps(framework, model, genome_reps=None):
    num_genomes = int(
        framework.genome_group().order() / framework.symmetry_group().order()
    )
    if genome_reps is None:
        genomes = list(framework.genomes().keys())
        genome_reps = genomes
    dists = {}
    for rep in genome_reps:
        a = prob_to_reach_in_steps_func(framework, model, rep)
        for k in range(0, num_genomes):
            if a(k) > 10 ** (-7):
                dists[rep] = k
                break
    return dists

def first_nonzero_value(framework, function, limit=None):
    if limit is None:
        limit = framework.genome_group().order() / framework.symmetry_group().order()
    first_k = limit
    for k in range(1, limit):
        if function(k) > 10 ** (-8):
            first_k = k
            break
    return first_k

def MFPT(framework, model, genome_reps=None, scale_by=1):
    """Return the mean time elapsed rearranging ref->genome where the target is an absorbing state"""
    genomes = framework.genomes()
    reps = list(genomes.keys())
    P = model.reg_rep_of_zs().toarray()
    P_0 = P.copy() 
    for i in range(len(P)):
        P_0[0, i] = 0
    Z = np.linalg.pinv(np.eye(len(P)) - P_0)

    def MFTP_dists():
        return (Z @ P_0 @ Z)[:, 0]

    MFTP_distances = list(MFTP_dists())
    MFTP_distances = {rep: scale_by * MFTP_distances[r] for r, rep in enumerate(reps)}
    if genome_reps is not None:
        MFTP_distances = {
            rep: MFTP_distances[framework.canonical_instance(rep)]
            for rep in genome_reps
        }
    return MFTP_distances

def fast_MFPT(framework, model):
    reg_rep, genomes = framework.fast_reg_rep_of_zs(model)
    Q = reg_rep[1:,1:] # Remove the absorbing state
    m = Q.shape[0] # n - 1
    A = scipy.sparse.identity(m) - Q
    scipy.sparse.linalg.cg(A, np.ones(m))

def dict_to_distance_matrix(distances, framework, genomes=None):
    """If need to convert to pairwise distances, supply a list of genomes."""
    if genomes is not None:
        D = np.zeros((len(genomes), len(genomes)))
        for i in range(len(genomes)):
            canonical_i = framework.canonical_instance(
                framework.random_instance(genomes[i])
            )
            distances_copy = {canonical_i * k: v for k, v in distances.items()}
            for j in range(len(genomes)):
                canonical_j = framework.canonical_instance(
                    framework.random_instance(genomes[j])
                )
                D[i, j] = distances_copy[canonical_j]
    else:
        D = np.zeros((len(distances), len(distances)))
        for (i, j), distance in distances.items():
            D[i, j] = distance
    return D


def genomes_for_dist_matrix(framework, genomes):
    instances = [framework.canonical_instance(g) for g in genomes]
    # Get the genomes g that we need dist(id -> g) for:
    need_distances = {}
    for i, _ in enumerate(genomes):
        for j, _ in enumerate(genomes):
            if i < j:
                canonical_i, canonical_j = instances[i], instances[j]
                need_distances[(i, j)] = canonical_i.inverse() * canonical_j
    return need_distances


def get_distance_function(distance_type):
    match distance_type:
        case DISTANCE.MFPT: return MFPT
        case DISTANCE.discrete_MFPT: return discrete_MFPT
        case DISTANCE.min: return min_distance_using_irreps
        case DISTANCE.MLE: return mles
        case DISTANCE.min_weighted: return min_distance


def distance_matrix(
    framework,
    model,
    genomes,
    distance,
    replace_nan_with=np.nan
):
    """
    Compute a distance matrix for a given set of genomes and distance measure.

    Args:
        framework (PositionParadigmFramework): the framework
        model (Model): the model
        genomes (list): a list of genomes to compute distances for
        distance (DISTANCE): the distance measure to use
    """
    need_distances = genomes_for_dist_matrix(framework, genomes)
    pairs, need_distances = map(list, zip(*need_distances.items()))

    # Parameters
    params = { "framework" : framework, "model" : model, "genome_reps" : need_distances }

    # Compute the distances:
    distance_func = get_distance_function(distance)
    distances = distance_func(**params)

    # Construct the distance matrix
    D = np.zeros((len(genomes), len(genomes)))
    for p, (i, j) in enumerate(pairs):
        D[i, j] = distances[need_distances[p]]

    if replace_nan_with != np.nan:
        D[np.isnan(D)] = replace_nan_with

    return D
