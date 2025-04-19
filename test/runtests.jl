

using Oscar
using CUDA

using Test
using BenchmarkTools

using Profile
# using MMPSingularities
#include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")
using MMPSingularities

#include("TestCases.jl")
#include("CalabiYauHeights.jl")
#include("QuasiFSplitMatrices.jl")
include("delta1_tests.jl")
include("height_tests.jl")

@testset "K3 surfaces" begin
  #test_heights_all()
  #test_matrices_all()
  height_run_tests()
  delta1_run_tests()
end

