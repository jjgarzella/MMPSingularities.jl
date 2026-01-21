module MMPSingularities

using Oscar
using Memoize
using Combinatorics
using StaticArrays
using BitIntegers
using InteractiveUtils
using CUDA
using AcceleratedKernels
using Primes
using SparseArrays
import Adapt

using GradedRingUtilities
using GPUPolynomials
using CudaNTTs
#include("../DeRham.jl/src/Utils.jl")

include("utils/GPUHashMap.jl")
include("utils/int128.jl")
include("utils/Delta1Helpers.jl")
include("utils/goldilocks.jl")

include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")
include("FPureThresholds.jl")

# using GPUFiniteFieldMatrices

include("QFSCalabiYau.jl")
include("QFSGeneralCase.jl")

include("RandomPolynomials.jl")
#include("PolyData.jl")

include("HarveyTrace/TraceFormula.jl")
include("HarveyTrace/GenericMultiplyThenSplit.jl")
include("HarveyTrace/diag_momts_naive_little.jl")

include("naivepointcounts/ProjectiveSpace.jl")
include("naivepointcounts/CountPoints.jl")

# exports here

end
