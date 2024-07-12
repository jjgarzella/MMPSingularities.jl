module MMPSingularities

using Oscar
using Memoize
using Combinatorics

include("../DeRham.jl/src/Utils.jl")

include("FrobSplittingInfra.jl")

include("../GPUPolynomials.jl/benchmarks/Benchmarks.jl")
include("../GPUPolynomials.jl/src/Delta1.jl")
using .Delta1

include("QFSCalabiYau.jl")
include("QFSGeneralCase.jl")

include("RandomPolynomials.jl")
include("PolyData.jl")

# exports here

end
