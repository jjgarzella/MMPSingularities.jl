module MMPSingularities

using Oscar
using CUDA
using Memoize
using Serialization

include("../GPUPolynomials.jl/src/GPUPolynomials.jl")
using .GPUPolynomials

include("RandomPolynomials.jl")
include("delta1/delta1.jl")
include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")
include("QFSCalabiYau.jl")

end
