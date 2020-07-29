def reflection_permutation(n):
	string = "(1,-1)"
	string += "".join(("(%s, -%s)(%s, -%s)" % (i+2, n-i, n-i, i+2) if i+2 != n-i else "(%s,-%s)" % (i+2, i+2)) for i in range(0,int(n/2)))			
	return string
	
def rotation_permutation(n):
	positive = tuple(i for i in range(1,n+1))
	negative = tuple(-i for i in positive)
	return str(positive) + str(negative)
	
def signed_inversion(of_positions = (1,2)):
	return "not done yet, will return a signed inversion in cycles form"