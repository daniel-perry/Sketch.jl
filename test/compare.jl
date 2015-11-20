# compare FD to a traditional rank-k approximation
# Danny Perry (dperry@cs.utah.edu)
# Oct 2014
#

using PyPlot
using Debug
using RDatasets

using Sketch

X = []
Labels = []

#name = "iris"
name = "baseball"

if name == "iris"
	iris = dataset("datasets","iris")

	X = [ convert(Array{Float64,1}, iris[1]) convert(Array{Float64,1}, iris[2]) convert(Array{Float64,1}, iris[3]) convert(Array{Float64,1}, iris[4]) ]
	n = size(X,1)
	println("X: ", size(X))
	Labels = ASCIIString[]
	for i=1:n
		push!(Labels, iris[5][i])
	end                                                                                                                                            
else

	baseball = dataset("plyr","baseball")

	numerical_cols = [1,4,8,9,10,11,12,13,17]
	X = zeros(size(baseball,1),length(numerical_cols))

	for i=1:length(numerical_cols)
		println(i,": ",numerical_cols[i])
		X[:,i] = convert(Array{Float64,1}, baseball[numerical_cols[i]])
	end

end

n = size(X,1)
X_centered = X .- (ones(n,1)*mean(X,1))


@debug function compare(X)

	n = size(X,1)
	d = size(X,2)

	errs = Float64[]
	minerrs = Float64[]
	maxerrs = Float64[]

	figure()
	for l=2:20
		B = FD(X,l)
		BB = B'*B
		U,S,V = svd(BB)
		
		projerr = sum((X .- (X*V*V')).^2,2)
		err = mean(projerr)
		minerr = minimum(projerr)
		maxerr = maximum(projerr)

		push!(errs,err)
		push!(minerrs,minerr)
		push!(maxerrs,maxerr)

		clf()
		plot(1:length(errs),errs,"r-")
		hold(:on)
		plot(1:length(errs),minerrs,"r--")
		plot(1:length(errs),maxerrs,"r--")
		title("FD projection error")
		xlabel("size of sketch")
		ylabel("proj err squared")
		sleep(0.0001)
	end

	@bp
end

compare(X_centered)

