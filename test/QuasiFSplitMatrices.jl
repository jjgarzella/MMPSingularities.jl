
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


