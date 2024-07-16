
include("../src/MMPSingularities.jl")

using Oscar
using CUDA

using Test
#using JET

function test_matrix_qfs_cy_char_2()
  R, (w,x,y,z) = polynomial_ring(GF(2),4)

  f = w^3*x + w^3*y + w^3*z + w^2*x^2 + w^2*y^2 + w^2*z^2 + w*x^3 + w*x^2*y + w*x*y^2 + w*y^2*z + x^3*z + x*y^3 + x*y^2*z + x*z^3 + y^4 + y^3*z + y^2*z^2 + y*z^3

  @test MMPSingularities.quasiFSplitHeight_CY_lift(2,f,10) == 3

  # Code from naive way of getting matrix:
  #
  #N = 4
  #p = 2
  #
  #fpminus1 = f^(p-1)
  #
  #Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  #
  #m = N*(p-1)
  #critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = vector(fpminus1,m)
  #println("creating matrix...")
  #@time M = MMPSingularities.matrix_of_lin_op(θFstar,m,parent(f))
  #println("matrix finished:")
  #display(M)
 
  # Code from new way of getting matrix:
  N = 4
  p = 2

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = MMPSingularities.Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)

  #println("creating matrix...")
  #=@time=# M = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(p,coefs,degs,m)

  #println("doing it again after jitter is primed")
  #@time M = MMPSingularities.matrix_of_multiply_then_split_sortmodp(p,coefs,degs,m)

  M_true = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 1  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  0  1  1  1  0  0  0  0  0  0  0  0  1  0  1  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  1  0  0  1  0  1  1  1  0  0  0  0  0  0  1  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  1  0  0  0  1  1  0  1  0  0  1  1  1  0  1  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  1  1  1  0  0  1  1  1  1  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  1  0  0  1  1  0  1  0  0  0  0  0  0  0  0  0  0
 1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  0  0  1  0  0  1  0  0  0  1  0  1  0  0  0  0  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
 1  1  0  1  0  0  1  1  0  0  0  0  0  0  0  1  0  0  1  0  1  1  0  0  0  1  1  0  1  0  0  0  1  1  0
 0  1  1  1  0  1  1  0  1  0  0  0  0  0  0  1  0  1  0  1  0  0  0  0  0  0  0  1  1  1  0  0  1  0  0
 0  0  0  1  1  0  0  1  1  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  1  0  0  1  0  1  0  1  0
 1  0  0  0  0  1  1  0  0  0  1  0  1  0  0  0  0  0  0  1  0  0  0  1  1  0  0  0  0  0  0  1  0  0  0
 0  0  1  0  0  0  0  1  1  1  1  0  0  1  0  0  1  0  0  0  1  1  1  0  0  0  0  0  1  0  1  1  1  1  0
 0  0  0  0  1  0  0  0  0  0  1  1  1  1  0  0  0  0  1  0  1  0  1  0  1  0  0  0  0  0  0  0  1  1  0
 0  1  0  0  0  1  1  0  0  1  0  0  1  1  0  1  1  0  0  0  0  0  0  0  1  1  0  0  0  1  0  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  1  0  0  1  1  0  1  0  0  1  0  0  0  1  1  0  0  0  0  1  0
 0  0  0  0  0  0  0  0  0  0  1  0  1  1  1  0  0  0  0  0  0  0  1  1  0  0  0  0  0  1  1  0  0  0  0
 0  1  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  1  1  0
 0  1  0  1  0  1  1  1  1  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0  1  1  0  0  1  0  1  0  0  1
 0  0  0  1  0  0  0  1  1  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  1  1  1  1  0  1  0  1  0
 0  1  0  0  0  1  1  0  0  0  1  0  1  1  0  0  1  0  0  1  1  0  0  0  0  1  0  0  1  1  0  0  0  1  0
 0  0  0  1  0  0  0  1  1  0  1  0  1  1  0  0  0  0  1  0  1  1  1  0  1  0  0  1  0  0  1  0  1  0  0
 0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  0  0  0  0  1  0  0  1  1  0  1  0  1  1  1  1  1  0  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  0  0  0  0  0  1  0  1  1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  0  0  0  1  0  1  1  0  1  0  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  0  1  1  0  0  1  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0]


  #display(M)
  #println("True result:")
  #display(M_true)

  @test M == M_true
end

function test_K3_5()

  best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(5,f,10)

  R, (x1,x2,x3,x4) = polynomial_ring(GF(5),4)

  ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4
  ftwo2 = x1^4 + x1^3*x2 + 2*x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + x1^2*x4^2 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 3*x1*x2^2*x4 + 2*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x4^3 + x2^3*x3 + x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 4*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + 2*x2*x4^3 + x3^2*x4^2 + x3*x4^3
  fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4
  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3


  @test best_qfs_height(ftwo1) == 2
  @test best_qfs_height(ftwo2) == 2
  @test best_qfs_height(fthree1) == 3
  @test best_qfs_height(ffour1) == 4
end

using Profile

function test_time_K3_5()
  qfs_height_fn = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu

  pregen = MMPSingularities.pregen_delta1(4,5)


  R, (x1,x2,x3,x4) = polynomial_ring(GF(5),4)

  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

  # prime the jitter
  println("PRIMING THE JITTER:")
  CUDA.@time qfs_height_fn(5,ffour1,10,pregen)

  println()
  #CUDA.@time qfs_height_fn(5,ffour1,10,pregen)
  #@report_opt qfs_height_fn(5,ffour1,10,pregen)


  println()
  Profile.Allocs.clear()
  Profile.Allocs.@profile qfs_height_fn(5,ffour1,10,pregen)
end


@testset "K3 surfaces" begin
  test_matrix_qfs_cy_char_2()
  test_K3_5()
  test_time_K3_5()
end
