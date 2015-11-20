# Column sampling sketch of a matrix
# LinearTimeSVD - Drineas, et  al. "Fast monte carlo algorithms for matrices ii: computing a low rank approximation to a matrix", SIAM Journal of Computing, 2006.
#
# Danny Perry (dperry@cs.utah.edu)
# April 2015

using Debug

@debug function LinearTimeSVD(A, num_cols)
	n = size(A,1)
	d = size(A,2)
	c = num_cols

	C = zeros(c,d)
	norm_est = 0
	for i=1:n
		normi = norm(A[i,:])^2
		norm_est += normi
		pc = normi/norm_est
		scale = 1/sqrt(c*pc)
		for j=1:c
			if rand() <= pc
				C[j,:] = scale * A[i,:]
			end
		end
	end
	U,S2,V = svd(C*C')
	S = sqrt(S2) # S are the approximate singular values of A
	H = (C' * V' * pinv(diagm(S)))  # H is an approximate right singular vector of A (left singular vector of A').
	return H,S
end

# same but using leverage scores
@debug function LeverageLinearTimeSVD(A, Uk,k, num_cols)
	n = size(A,1)
	d = size(A,2)
	c = num_cols

	row_lev = sum(Uk.^2,2) ./ k

	C = zeros(c,d)
	for i=1:n
		pc = row_lev[i]
		scale = 1/sqrt(c*pc)
		for j=1:c
			if rand() <= pc
				C[j,:] = scale * A[i,:]
			end
		end
	end
	U,S2,V = svd(C*C')
	S = sqrt(S2) # S are the approximate singular values of A
	H = (C' * V' * pinv(diagm(S)))  # H is an approximate right singular vector of A (left singular vector of A').
	return H,S
end


