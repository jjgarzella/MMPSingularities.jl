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
include("QFSCalabiYau.jl")

end
