def regular_representation(G, g): 
	""" regular representation """
	C = QQbar # Complex field
	CG =  GroupAlgebra(G,C)
	return matrix(QQbar, [(CG(g)*G(h)).to_vector(QQbar) for h in G]).transpose()