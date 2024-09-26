module MMPSingularities

using Oscar
using GPUPolynomials
using CUDA

include("RandomPolynomials.jl")
include("GPUDelta1.jl")
include("FrobSplittingInfra.jl")
include("MatricesOfSplittings.jl")
include("QFSCalabiYau.jl")

end
