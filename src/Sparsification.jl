# Sparsification sketch of a matrix
#
# Danny Perry (dperry@cs.utah.edu)
# April 2015
#

using Debug

# Random sparsification that preserves eigen spectrum of A
# Theorem 1 from Achlioptas, "Fast computation of low rank matrix approximations", STOC 2001
# http://web.stanford.edu/class/ee378b/papers/achlioptas.pdf
#
# A - matrix to sparsify
# sparsity - controls amount sparsity, 1 <= s <= (m+n)/(log(m+n)^6)
#            sparse terms will be chosen w.p. (1-1/s)
#
function RandomSparsify(A, sparsity)
	n = size(A,1)
	d = size(A,2)
	s = sparsity

	entries = n*d
	I = zeros(Int32,int(.5*entries))
	J = zeros(Int32,int(.5*entries))
	V = zeros(typeof(A[1,1]), int(.5*entries))
	k = int(0)

	for i=1:n
		for j=1:d
			if rand() < 1.0/s
				#B[i,j] = s * A[i,j]
				k += 1
				if k > size(I,1)
					append!(I, zeros(typeof(I[1]),int(.25*size(I,1))))
					append!(J, zeros(typeof(J[1]),int(.25*size(J,1))))
					append!(V, zeros(typeof(V[1]),int(.25*size(V,1))))
				end
				I[k] = i
				J[k] = j
				V[k] = s*A[i,j]
			#else
			# B[i,j] = 0  # implied via sparse matrix
			end
		end
	end
	println("sparse k: ", k)

	I = I[1:k]
	J = J[1:k]
	V = V[1:k]

	B = sparse(I,J,V,n,d)
end


# Random sparsification that preserves eigen spectrum of A
# Theorem 2 from Achlioptas, "Fast computation of low rank matrix approximations", STOC 2001
# http://web.stanford.edu/class/ee378b/papers/achlioptas.pdf
#
# A - matrix to sparsify
#
# returns sparse binary matrix B with non-zero entries = 2, 
# however estimates should use matrix  b*(B .- 1), so that each 
# entry is either +b or -b.
#
# returns B,b
function RandomRoundSparsify(A)
	n = size(A,1)
	d = size(A,2)

	entries = n*d
	I = zeros(Int32,int(.5*entries))
	J = zeros(Int32,int(.5*entries))
	V = zeros(typeof(A[1,1]), int(.5*entries))
	k = int(0)

	b = maximum(abs(A))
	for i=1:n
		for j=1:d
			p = .5 + A[i,j] / (2*b)
			if rand() < p
				#B[i,j] = 2
				k += 1
				if k > size(I,1)
					append!(I, zeros(typeof(I[1]),int(.25*size(I,1))))
					append!(J, zeros(typeof(J[1]),int(.25*size(J,1))))
					append!(V, zeros(typeof(V[1]),int(.25*size(V,1))))
				end
				I[k] = i
				J[k] = j
				V[k] = 2

			#else
			# B[i,j] = 0  # implied via sparse matrix
			end
		end
	end

	I = I[1:k]
	J = J[1:k]
	V = V[1:k]

	B = sparse(I,J,V,n,d)
	B,b
end


