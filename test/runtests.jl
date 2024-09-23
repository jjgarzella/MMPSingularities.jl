

using Oscar
using CUDA

using Test

using Profile
#using JET

#include("../src/MMPSingularities.jl")
using MMPSingularities

include("QuasiFSplitMatrices.jl")
include("CalabiYauHeights.jl")

@testset "K3 surfaces" begin
  test_matrix_qfs_cy_char_2()
  test_high_height_CY3_char_2()
  
  test_matrix_qfs_cy_char_3()
  test_delta1_char_3()
  #test_K3_3()
  #test_K3_3_gpu()

  test_K3_5()
  test_K3_5_matrix()
  test_time_K3_5()
end
