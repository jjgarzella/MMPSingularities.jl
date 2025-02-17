

using Oscar
using CUDA

using Test

using Profile
# using MMPSingularities
include("../src/MMPSingularities.jl")
using .MMPSingularities

using Revise
includet("TestCases.jl")
includet("CalabiYauHeights.jl")
includet("QuasiFSplitMatrices.jl")

@testset "K3 surfaces" begin
  test_heights_all()
  test_matrices_all()
end
