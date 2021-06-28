"""
"""
import numpy as np
import networkx as nx
from sage.calculus.predefined import W
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

def mle(framework, model, genomes):
    """Return dictionary of MLEs (if they exist) for time elapsed rearranging ref->genome for each genome"""
    reps = framework.irreps()
    G = framework.genome_group()
    Z = framework.symmetry_group()
    num_reps = len(reps)
    z_coefficient = 1/Z.order()
    UCF = UniversalCyclotomicField()
    #irreps_of_s = [sum([s_coefficients[sigma]*rep(sigma, as_gap_matrix=True) for sigma in M]) for rep in reps ]
    #irreps_of_z = [sum([z_coefficient*rep(d, as_gap_matrix=True) for d in D]) for rep in reps ]
    #irreps_of_zs = [ Matrix(F, matrix(UCF irrep_d*irrep_s)) for irrep_d, irrep_s in zip(irreps_of_z, irreps_of_s)]
    #irreps_of_s = [ Matrix(F, matrix(UCF, irrep)) for irrep in irreps_of_s]
    #irreps_of_z = [ Matrix(F, matrix(UCF, irrep)) for irrep in irreps_of_z]
    #for irrep in irreps_of_zs:
	#    irrep.set_immutable()

    #eigs = [sorted(set([round(real(eig),9) for eig in irrep_of_zs.eigenvalues()]).list()) for irrep_of_zs in irreps_of_zs]
    #numEs = [len(eigs) for eigs in eigs] 

    #print("COMPUTE PROJECTION OPERATORS")
    #dims = [irrep_of_zs.nrows() for irrep_of_zs in irreps_of_zs]
    #projections=[[Matrix.identity(dim) for _ in range(numEs[d])] for d, dim in enumerate(dims)] # Identity matrices for now
    #for rep_index, proj in enumerate(projections):
    #    for j in range(0, numEs[rep_index]):
    #        for k in range(0, numEs[rep_index]):
    #            if k!=j:
    #                dim = irreps_of_zs[rep_index].nrows()
    #                proj[j]=proj[j]*(irreps_of_zs[rep_index]-eigs[rep_index][k]*Matrix.identity(dim))*(1/(eigs[rep_index][j]-eigs[rep_index][k]))


    #%% Compute equivalence classes for genomes, ASSUMING DIHEDRAL SYMMETRY AND TIME REVERSIBILITY ###
    #!%%time
    #print("COMPUTE GENOME EQUIVALENCE CLASSES")
    #classlist = hyperoctahedral_groups.EquivalenceClasses(H, n)
    #num_classes = len(classlist) # Number of equivalence classes
    
    ### 'traces' is indexed by irreps, and each matrix is indexed by equivalence class (row) and eigenvalue (column) ###
    #traces = [matrix(CDF, num_classes, numE) for numE in numEs]
    #traces = [ [ [ 0 for e in range(numE) ] for c in range(num_classes) ] for numE in numEs]
    #inverted_sigmas = [ H(conversions.signed_permutation_to_cycles(n, each_class[0].inverse(), signed=signed)) for each_class in classlist ]
    #for r in range(num_reps):
    #    print(f'Representation {r+1} of the Hyperoctahedral Group, of dimension {dims[r]} with eigenvalues: {[round(real(eig),4) for eig in eigs[r]]}')
    #    for c in range(num_classes):
    #        sigd = Matrix(F, matrix(UniversalCyclotomicField(), reps[r](inverted_sigmas[c], as_gap_matrix=True))) * irreps_of_z[r]
    #        for e in range(len(eigs[r])):
    #            if eigs[r] == [0]: # then the projection matrix is the identity, so...
    #                traces[r][c][e] = round(real(sigd.trace()),4)
    #            else:
    #                traces[r][c][e] = round(real((sigd*projections[r][e]).trace()),4)

    #%% SET UP LIKELIHOOD FUNCTION ###
    #!%%time
    # print("SET UP LIKELIHOOD FUNCTION")
    # genome_class_reps = [equiv_class[0] for equiv_class in classlist]
    # def likelihood(t, genome):
    #     ans = 0
    #     genome_index = genome_class_reps.index(genome)
    #     for r in range(len(irreps_of_zs)):
    #         ans += D.order()*(exp(-t)/H.order())*dims[r]*sum(exp(eigenvalue*CDF(t)) * traces[r][genome_index][e] for e, eigenvalue in enumerate(eigs[r]))
    #     return real(ans)

    # #%% PLOT LIKELIHOOD FUNCTIONS AND COMPUTE MLES
    # print("Print likelihood functions and compute MLEs")
    # mle_list = dict()
    # bound = (0, 35)
    # path = f'output/figures_test/n{n}_{settings["model"]}_{bound[0]}to{bound[1]}'
    # step = 0.5
    # os.makedirs(path)
    # for genome in genome_class_reps:
    #     print(f'Genome {genome}')
    #     times = np.arange(bound[0], bound[1], step)
    #     likelihood_functional_values = [likelihood(time, genome) for time in times]
    #     #likelihood_functional_values_interpolated = interp1d(times, likelihood_functional_values, kind='cubic')
    #     plt.plot(times, likelihood_functional_values)
    #     #plt.plot(times, likelihood_functional_values_interpolated(times), color='green')
    #     plt.savefig(path + f'/{genome}')
    #     plt.show()
    #     def L(time):
    #         return -likelihood(time, genome)
    #     mle = minimize(L, 0.1, method="TNC", bounds=[bound])
    #     mle_list[genome] = "NA" if mle.x[0] >= bound[1] and mle.x[0] > 0.0001 else mle.x[0] # CHECK THIS!!!
    #     plt.clf()




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