from sage.all_cmdline import *
from cgt import *

##calculates partial traces for all equivalence classes (reps to be loaded in from file) for up to N=8 regions
n = 4
# S = SymmetricGroup(n)
H = SymmetricGroup(n) #hyperoctahedral_groups.HyperoctahedralGroup(n)
D = hyperoctahedral_groups.DihedralSubgroup(H, n)

#make a list of the irreps
reps = representations.irreducible_representations(n, signed=False) 
num_reps = len(reps)

# Model 1
# M = [H(f'({i},{-1*(i+1)})({-1*i},{i+1})') for i in range(1,n)] + [H(f'({n},-1)(-{n},1)')]
M = [H("(1,2)"), H("(2,3)"), H("(3,4)"), H("(1,4)")]
ss = tuple(m for m in M)
dd = tuple(d for d in D)

#make a list of the irreps of s
irreps_of_s = [ sum([rep(sigma) for sigma in ss]) for rep in reps ]

#make a list of the dimensions of the irreps
dims=[irrep_of_s.nrows() for irrep_of_s in irreps_of_s]

#make these matrices immutable
for irrep in irreps_of_s:
	irrep.set_immutable()

#save(irrepss,'irrepss'+str(N))
#save(dims,'dims'+str(N))


#find the eigenvalues of the irreps as a sensible sorted list
Elist=[sorted(Set(irrep_of_s.eigenvalues()).list()) for irrep_of_s in irreps_of_s]

#how many eigenvalues does each irrep have?
numEs=[len(eigs) for eigs in Elist]

#projections.... go!
projections=[[matrix.identity(dim) for _ in range(numEs[d])] for d, dim in enumerate(dims)] ##this is the list skeleton, filled with identity matrices

for rep_index, proj in enumerate(projections):
	for j in range(0, numEs[rep_index]):
		for k in range(0, numEs[rep_index]):
			if k!=j:
				proj[j]=proj[j]*(irreps_of_s[rep_index]-(Elist[rep_index][k]))*(1/(UniversalCyclotomicField()(Elist[rep_index][j]-Elist[rep_index][k])))

##we also need the irreps of d
irreps_of_d = [sum([rep(sigma) for sigma in dd]) for rep in reps ]


## Equivalence classes
classlist = hyperoctahedral_groups.EquivalenceClasses(H, n)
num_classes = len(classlist)

##now calculate the partial traces
traces=[ matrix(RDF, num_classes, numE) for numE in numEs ]
tracef=open('partial_traces_S'+str(n)+'.txt','w')
i = 0
while i < num_reps:
	reps[i]
	print("dimension", dims[i])
	print("e-values", Elist[i])
	print("partial-traces")
	tracef.write(str(reps[i]))
	tracef.write(str('\n'))
	tracef.write(str("dimension ")+str(dims[i]))
	tracef.write(str('\n'))
	tracef.write(str("e-values ")+str(Elist[i]))
	tracef.write(str('\n'))
	tracef.write(str("partial-traces "))
	tracef.write(str('\n'))
	k=0
	while k<num_classes:
		sigd = reps[i](H(conversions.signed_permutation_to_cycles(n, classlist[k][0], signed=False))) * irreps_of_d[i]
		j=0
		while j<len(Elist[i]):
			traces[i][k,j]=real((sigd*projections[i][j]).trace())
			j=j+1
		tracef.write(str(traces[i][k,:]))
		tracef.write(str('\n'))
		k=k+1
	print(i)
	i=i+1
tracef.close()