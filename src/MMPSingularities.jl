module MMPSingularities

using Oscar
using Memoize
using Combinatorics
using StaticArrays

using DeRham
#include("../DeRham.jl/src/Utils.jl")

include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")

#include("../GPUPolynomials.jl/benchmarks/Benchmarks.jl")
#include("../GPUPolynomials.jl/src/Delta1.jl")
#using .Delta1
using GPUPolynomials
using GPUFiniteFieldMatrices

include("QFSCalabiYau.jl")
include("QFSGeneralCase.jl")

include("RandomPolynomials.jl")
include("PolyData.jl")

# exports here

end
