# Sketch algorithms
#
# Danny Perry (dperry@cs.utah.edu)
# Oct 2014
#

module Sketch

export 
	# Frequent Directions
	FD, 
	# Column Sampling
	LinearTimeSVD, LeverageLinearTimeSVD, 
	# CUR
	LinearTimeCUR, LeverageCUR,
	# Sparsification
	RandomSparsify, RandomRoundSparsify


# Frequent Directions
include("FD.jl")

# Column sampling
include("ColumnSampling.jl")

# CUR
include("CUR.jl")

# Sparsify
include("Sparsification.jl")

end # module

