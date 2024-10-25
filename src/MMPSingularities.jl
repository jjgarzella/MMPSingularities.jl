module MMPSingularities

using GPUPolynomials
using Oscar
using Memoize
using Combinatorics
using StaticArrays
using CUDA
using Serialization
using InteractiveUtils

#include("../DeRham.jl/src/Utils.jl")



include("delta1/delta1.jl")

include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")


# using GPUFiniteFieldMatrices

include("QFSCalabiYau.jl")
include("QFSGeneralCase.jl")

#include("RandomPolynomials.jl")
#include("PolyData.jl")

# exports here

end
