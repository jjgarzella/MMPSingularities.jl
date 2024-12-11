module MMPSingularities

using Oscar
using CUDA
using Memoize
using Serialization
using Primes
using StaticArrays
using BitIntegers
import Adapt

include("../GPUPolynomials.jl/src/GPUPolynomials.jl")
using .GPUPolynomials

include("RandomPolynomials.jl")
include("delta1/delta1.jl")
include("FrobSplittingInfra.jl")
include("gpu_hashmap.jl")
include("MatricesOfSplittings.jl")
include("FPureThresholds.jl")

#include("../GPUPolynomials.jl/benchmarks/Benchmarks.jl")
#include("../GPUPolynomials.jl/src/Delta1.jl")
#using .Delta1
#using GPUPolynomials
#using GPUFiniteFieldMatrices

include("QFSCalabiYau.jl")

end
