# CUR decomposition sketch
# LinearTimeCUR - Drineas, et  al. "Fast monte carlo algorithms for matrices iii: computing a compressed approximate matrix decomposition", SIAM Journal of Computing, 2006.
#
# Danny Perry (dperry@cs.utah.edu)
# April 2015

using Debug

# CUR decomposition (rank-k)
# A - the matrix to approximate
# fro_norm - frobenius norm of A
# num_cols - number of colums to sample for C (and U)
# num_rows - number of rows to sample for R (and U)
# k - rank k to approximate
#
# returns C,U,R matrices
@debug function LinearTimeCUR(A,fro_norm, num_cols,num_rows,k)
	n = size(A,1)
	d = size(A,2)

	c = num_cols
	r = num_rows

	C = zeros(n,c)
	U = zeros(c,r)
	R = zeros(r,d)

	for i=1:d
		pc = norm(A[:,i])^2 / fro_norm^2
		scale = 1/sqrt(c*pc)
		for j=1:c
			if rand() <= pc
				C[:,j] = scale * A[:,i]
			end
		end
	end
	CU,CS,CV = svd(C'*C)
	k = min(k, rank(diagm(CS)))
	for i=1:n
		pc = norm(A[i,:])^2 / fro_norm^2
		scale = 1/sqrt(r*pc)
		for j=1:r
			if rand() <= pc
				R[j,:] = scale * A[i,:]
				U[:,j] = scale * C[i,:]
			end
		end
	end
	U = CV[:,1:k] * pinv(diagm(CS[1:k])) * CU[:,1:k]' * U

	C,U,R
end

# This uses a decomposition to choose samples based on the leverage scores
# Uk,Vk - the first k columns of the matrix U and V, such that  A = USV'
@debug function LeverageCUR(A,fro_norm,Uk,Vk, num_cols,num_rows,k)
	n = size(A,1)
	d = size(A,2)

	c = num_cols
	r = num_rows

	C = zeros(n,c)
	U = zeros(c,r)
	R = zeros(r,d)

	col_lev = sum(Vk.^2,2) ./ k
	row_lev = sum(Uk.^2,2) ./ k

	for i=1:d
		pc = col_lev[i]
		scale = 1/sqrt(c*pc)
		for j=1:c
			if rand() <= pc
				C[:,j] = scale * A[:,i]
			end
		end
	end
	CU,CS,CV = svd(C'*C)
	k = min(k, rank(diagm(CS)))
	for i=1:n
		pc = row_lev[i]
		scale = 1/sqrt(r*pc)
		for j=1:r
			if rand() <= pc
				R[j,:] = scale * A[i,:]
				U[:,j] = scale * C[i,:]
			end
		end
	end
	U = CV[:,1:k] * pinv(diagm(CS[1:k])) * CU[:,1:k]' * U

	C,U,R
end


