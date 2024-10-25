module MMPSingularities

include("../GPUPolynomials.jl/src/GPUPolynomials.jl")
using .GPUPolynomials
using Oscar
using CUDA
using Memoize
using Serialization

include("RandomPolynomials.jl")
include("delta1/delta1.jl")
include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")
include("QFSCalabiYau.jl")

end
