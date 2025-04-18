

using Oscar
using CUDA

using Test
using Profile
using BenchmarkTools

# using MMPSingularities
#include("../src/MMPSingularities.jl")
using MMPSingularities

#include("TestCases.jl")
#include("CalabiYauHeights.jl")
#include("QuasiFSplitMatrices.jl")
include("delta1_tests.jl")
include("height_tests.jl")

#@testset "K3 surfaces" begin
#  test_heights_all()
#  test_matrices_all()
#end

@testset "new_algorithms" begin
    delta1_tests()
    height_tests()
end
