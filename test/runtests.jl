

using Oscar
using CUDA

using Test

using Profile
#using JET

#include("../src/MMPSingularities.jl")
using MMPSingularities

#include("TestCases.jl")
#include("CalabiYauHeights.jl")
#include("QuasiFSplitMatrices.jl")

using Revise
includet("TestCases.jl")
includet("CalabiYauHeights.jl")
includet("QuasiFSplitMatrices.jl")



@testset "K3 surfaces" begin
  #test_heights_all()
  #test_matrices_all()
end
