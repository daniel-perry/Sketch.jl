# Frequent Directions
# from: Liberty, "Simple and Deterministic Matrix Sketching"
#
# Danny Perry (dperry@cs.utah.edu)
# Oct 2014

using Debug

function reduce(S,V,p)
	delta = S[end-p+1]^2
	SS = (S .^2) - delta
	SS[find(SS .< 0)] = 0
	diagm( sqrt(SS) ) * V'
end

# Frequent Directions
# @param A - original matrix
# @param l - size of sketch
# @param p - the number of components to zero out at a time
# 
# @return B - the sketch
function FD(A, l, p)
  n = size(A,1)
  d = size(A,2)
  if l < p
    error("l must be >= p")
  end
  if p > d
    error ("p must be <= d")
  end
	B = A[1:l,:]
	for i=l+1:p:n
    numz = p
    if n-i < p
      numz = n-i+1
    end
		U,S,V = svd(B)
		corankB = length(find(S .<= eps(typeof(S[1]))))
		if corankB < numz # d-rank(B) < numz
    	B = reduce(S,V,numz)
		else
			B = diagm(S)*V'
		end
    B[end-numz+1:end,:] = A[i:i+numz-1,:]
	end
	B
end

# Frequent Directions - bit by bit
# @param A - next p or less rows
# @param B - the current B matrix
# @param l - size of sketch
# @param p - the number of components to zero out at a time
# 
# @return B - the sketch
function FD(A, B, l, p)
	numz = size(A,1)
  d = size(A,2)
  if l < p
    error("l must be >= p")
  end
  if p > d
    error ("p must be <= d")
  end
	U,S,V = svd(B)
	if length(find(S .<= eps(typeof(S[1])))) < numz # d-rank(B) < numz
		B = reduce(S,V,numz)
	else
		B = diagm(S)*V'
	end
	B[end-numz+1:end,:] = A
	B
end
