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

using GPUFiniteFieldMatrices

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

global fpminus1_time = Float64[]
global delta1_time = Float64[]
global move_matrix_time = Float64[]
global move_vector_time = Float64[]
global stripe_mul_time = Float64[]
global if_time = Float64[]

function reset_times()
    global fpminus1_time = Float64[]
    global delta1_time = Float64[]
    global move_matrix_time = Float64[]
    global move_vector_time = Float64[]
    global stripe_mul_time = Float64[]
    global if_time = Float64[]
end

function display_times()
  # s(vec) = begin
  #   if length(vec) == 0
  #     return 0
  #   else
  #     return sum(vec)
  #   end
  # end

  println("fpminus1 avg time: $(sum(fpminus1_time) / length(fpminus1_time))")
  println("delta1 avg time: $(sum(delta1_time) / length(delta1_time))")
  println("move_matrix avg time: $(sum(move_matrix_time) / length(move_matrix_time))")
  println("move_vector avg time: $(sum(move_vector_time) / length(move_vector_time))")
  println("stripe_mul avg time: $(sum(stripe_mul_time) / length(stripe_mul_time))")
  println("if avg time: $(sum(if_time) / length(if_time))")
  println("total time: $(sum(fpminus1_time) + sum(delta1_time) + sum(move_matrix_time) + sum(move_vector_time) + sum(stripe_mul_time) + sum(if_time))")
  println("total itrs: $(length(fpminus1_time) + length(delta1_time) + length(move_matrix_time) + length(move_vector_time) + length(stripe_mul_time) + length(if_time))")
  println("avg time: $((sum(fpminus1_time) + sum(delta1_time) + sum(move_matrix_time) + sum(move_vector_time) + sum(stripe_mul_time) + sum(if_time)) / (length(fpminus1_time) + length(delta1_time) + length(move_matrix_time) + length(move_vector_time) + length(stripe_mul_time) + length(if_time)))")
  println()
  
  println("fpminus times: $(fpminus1_time)")
  println("delta1 times: $(delta1_time)")
  println("move_matrix times: $(move_matrix_time)")
  println("move_vector times: $(move_vector_time)")
  println("stripe_mul times: $(stripe_mul_time)")
  println("if times: $(if_time)")
end

end
