
include("../src/MMPSingularities.jl")

using Oscar
using CUDA

using Test

using Profile
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

function test_matrix_qfs_cy_char_3()
  p = 3 
  N = 4

  R, (w,x,y,z) = polynomial_ring(GF(p),4)

  f = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w

  @test MMPSingularities.quasiFSplitHeight_CY_lift(p,f,10) == 3

  # Code from naive way of getting matrix:
  
  
  fpminus1 = f^(p-1)
  
  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  
  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = vector(fpminus1,m)
  println("creating matrix...")
  @time M_true = MMPSingularities.matrix_of_lin_op(θFstar,m,parent(f))
  #println(M)
 
  # Code from new way of getting matrix:

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = MMPSingularities.Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)

  #println("creating matrix...")
  #=@time=# M = MMPSingularities.matrix_of_multiply_then_split_sortmodp_dict(p,coefs,degs,m)

  #println("doing it again after jitter is primed")
  #@time M = MMPSingularities.matrix_of_multiply_then_split_sortmodp(p,coefs,degs,m)

#  M_true = 


  println("matrix finished:")
  display(M)
  println("True result:")
  display(M_true)

  @test M == M_true
end

function test_high_height_CY3_char_2()
  p = 2
  
  R, (u,w,x,y,z) = polynomial_ring(GF(p),["u","w","x","y","z"])
  
  f60 = x^5 + y^5 + z^5 + w^5 + u^5 + x*z^3*w + y*z*w^3 + x^2*z*u^2 + y^2*z^2*w + x*y^2*w*u + y*z*w*u^2

  @test best_qfs_height(f60) == 60
end

function test_K3_3()

  p = 3
  best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift_matrix(3,f,10)

  R, (w,x,y,z) = polynomial_ring(GF(3),["w","x","y","z"])
 
  f1 = x^4 + y^4 + z^4 + 2w^4 + x^2*y*w + y*z^2*w
 2*w^4 + w*x^2*y + w*y*z^2 + x^4 + y^4 + z^4
  f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
 2*w^4 + x^4 + x*y*z^2 + 2*y^4 + 2*z^4
  f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w
 w^4 + w*z^3 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
  f4 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2
 w^4 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
  f5 = x^4 + y^4 + z^4 + w^4 + x^3*z + z^3*w + y*z^2*w + y*z*w^2
 w^4 + w^2*y*z + w*y*z^2 + w*z^3 + x^4 + x^3*z + y^4 + z^4
  f6 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3
 w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4
  f7 = x^4 + y^4 + z^4 + w^4 + x*y^2*z + x*z^2*w + y*z^2*w + y*z*w^2
 w^4 + w^2*y*z + w*x*z^2 + w*y*z^2 + x^4 + x*y^2*z + y^4 + z^4
  f8 = x^4 + x^2*y*z + x^2*y*w + 2x^2*z^2 + x*y*w^2 + 2*y^4 + y^3*w + z^4 + w^4
 w^4 + w^2*x*y + w*x^2*y + w*y^3 + x^4 + x^2*y*z + 2*x^2*z^2 + 2*y^4 + z^4
  f9 = x^4 + y^4 + z^4 + w^4 + x*y^3 + y^3*w + z^2*w^2 + 2*x*y*z^2 + y*z*w^2
 w^4 + w^2*y*z + w^2*z^2 + w*y^3 + x^4 + x*y^3 + 2*x*y*z^2 + y^4 + z^4
  f10 = x^4 + 2*x^2*y*z + x^2*y*w + x*y^2*w + y^4 + y^3*w + y^2*z^2 + 2*y^2*z*w + y^2*w^2 + y*z^3 + y*z^2*w + y*z*w^2 + z^4 + z*w^3
 w^3*z + w^2*y^2 + w^2*y*z + w*x^2*y + w*x*y^2 + w*y^3 + 2*w*y^2*z + w*y*z^2 + x^4 + 2*x^2*y*z + y^4 + y^2*z^2 + y*z^3 + z^4
  finfty = x^4 + y^4 + z^4 + w^4
 w^4 + x^4 + y^4 + z^4
  fm10 = 2w^4 + x*w^3 + w^2*x^2 + w^2*x*y + w^2*y*z + 2w*x^3 + 2w*x^2*y + w*x*y^2 + w*x*y*z + w*x*z^2 + w*y^3 + 2w*y^2*z + w*y*z^2 + x^4 + 2*x^3*y + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + 2x*y*z^2 + x*z^3 + 2*y^4 + 2*y^3*z + y^2*z^2 + y*z^3 + 2*z^4

  (x1,x2,x3,x4) = (w,x,y,z)
  fr8 = x1^3*x2 + 2*x1^2*x2^2 + 2*x1^2*x2*x3 + x1^2*x3^2 + x1^2*x3*x4 + 2*x1*x2^2*x3 + 2*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + x2^4 + x2^3*x3 + 2*x2^3*x4 + x2^2*x3^2 + 2*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + x2*x4^3 + x3^4 + x3^2*x4^2 + x3*x4^3


  @test best_qfs_height(f1) == 1
  @test best_qfs_height(f2) == 2
  @test best_qfs_height(f3) == 3
  @test best_qfs_height(f4) == 4
  @test best_qfs_height(f5) == 5
  @test best_qfs_height(f6) == 6
  @test best_qfs_height(f7) == 7
  @test best_qfs_height(f8) == 8
  @test best_qfs_height(f9) == 9
  @test best_qfs_height(f10) == 10
  infty_substitute = best_qfs_height(finfty)
  @test infty_substitute == 11 || infty_substitute == 12
  @test best_qfs_height(fm10) == 10 
  @test best_qfs_height(fr8) == 8 
end

function test_K3_5()

  p = 5
  best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p,f,10)

  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),4)

  ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4
  ftwo2 = x1^4 + x1^3*x2 + 2*x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + x1^2*x4^2 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 3*x1*x2^2*x4 + 2*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x4^3 + x2^3*x3 + x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 4*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + 2*x2*x4^3 + x3^2*x4^2 + x3*x4^3
  fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4
  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3


  @test best_qfs_height(ftwo1) == 2
  @test best_qfs_height(ftwo2) == 2
  @test best_qfs_height(fthree1) == 3
  @test best_qfs_height(ffour1) == 4


end

function test_K3_5_matrix()

  p = 5 
  N = 4

  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),4)


  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

  @test MMPSingularities.quasiFSplitHeight_CY_lift(p,ffour1,10) == 4

  f = ffour1
  fpminus1 = f^(p-1)

  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = MMPSingularities.Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)

  #println("creating matrix...")
  #=@time=# 
  M_true = MMPSingularities.matrix_of_multiply_then_split_sortmodp_dict(p,coefs,degs,m)

  M = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker_noslice(p,coefs,degs,m)

  @test M == M_true

end

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
  #Profile.Allocs.clear()
  #Profile.Allocs.@profile qfs_height_fn(5,ffour1,10,pregen)

  Profile.clear()
  Profile.@profile qfs_height_fn(5,ffour1,10,pregen)
end


@testset "K3 surfaces" begin
  #test_matrix_qfs_cy_char_2()
  #test_high_height_CY3_char_2()
  #test_matrix_qfs_cy_char_3()
  #test_K3_3()
  #test_K3_5()
  #test_K3_5_matrix()
  test_time_K3_5()
end
