module MMPSingularities

using Oscar
using Memoize
using Combinatorics
using StaticArrays
using BitIntegers
using InteractiveUtils
using CUDA
using Primes
import Adapt

using GradedRingUtilities
using GPUPolynomials
#include("../DeRham.jl/src/Utils.jl")

include("utils/gpuhashmap.jl")
include("utils/int128.jl")

include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")
include("FPureThresholds.jl")

# using GPUFiniteFieldMatrices

include("QFSCalabiYau.jl")
include("QFSGeneralCase.jl")

#include("RandomPolynomials.jl")
#include("PolyData.jl")

# exports here

end
