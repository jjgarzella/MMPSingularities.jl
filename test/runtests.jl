

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
# include("QuasiFSplitMatrices.jl")
include("delta1_tests.jl")
include("height_tests.jl")

@testset "K3 surfaces" begin
  # test_heights_all()
  # test_matrices_all()
  height_run_tests()
  delta1_run_tests()
end

# TODO: fix this test and put it in its own file.
# function test_qfs_general_case()
#     R, (w,x,y,z,u) = polynomial_ring(GF(7),["w","x","y","z","u"])
#     f = w^4 + x^4 + y^4 + z^4 + u^4
#     h = MMPSingularities.quasiFSplitHeight_lift_mingens_wics(p,f,3)
#
#     @test h == 2
# end
