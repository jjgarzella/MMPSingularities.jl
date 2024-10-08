
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
  (coefs,degs) = MMPSingularities.convert_to_gpu_representation(Δ₁fpminus1)

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


function test_matrix_K3(p,f,sortmodp_function)
  N = 4

  fpminus1 = f^(p-1)
  
  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  println("Delta1 number of terms: $(length(terms(Δ₁fpminus1)))")
  θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  
  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = vector(fpminus1,m)
  (coefs,degs) = MMPSingularities.convert_to_gpu_representation(Δ₁fpminus1)
  display(degs)

  println("creating matrix...")
  @time M_true = MMPSingularities.matrix_of_lin_op(θFstar,m,parent(f))
  #println(M)
 
  M_true = Int64.(map(x -> lift(ZZ,x),M_true))
  # Code from new way of getting matrix:

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
  println("Delta1 number of terms: $(length(terms(Δ₁fpminus1)))")
  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = MMPSingularities.convert_to_gpu_representation(Δ₁fpminus1)
  display(degs)

  #println("creating matrix...")
  #=@time=# M = sortmodp_function(p,coefs,degs,m)

  #println("doing it again after jitter is primed")
  #@time M = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(p,coefs,degs,m)

#  M_true = 


  println("matrix finished:")
  display(M)
  println("True result:")
  display(M_true)

  println(findall(==(0),(M .== M_true)))
  @test M == M_true

end

function test_k3_3_matrix()
  #p = 3 

  #R, (w,x,y,z) = polynomial_ring(GF(p),4)

  #f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w

  #test_matrix_K3(p,f3)

  #f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
  cases_3, _ = test_cases_k3surf(3)
  f2 = cases_3[2]
  f3 = cases_3[3]
  
  test_matrix_K3(3,f2,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker)
  test_matrix_K3(3,f3,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker)
end

function test_k3_5_matrix()

  cases_5, _ = test_cases_k3surf(5)
  ffour1 = cases_5[4]

  test_matrix_K3(5,ffour1,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker)

#  p = 5 
#  N = 4
#
#  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),4)
#
#
#  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3
#
#  @test MMPSingularities.quasiFSplitHeight_CY_lift(p,ffour1,10) == 4
#
#  f = ffour1
#  fpminus1 = f^(p-1)
#
#  Δ₁fpminus1 = MMPSingularities.Δ₁l(p,fpminus1)
#  #θFstar(a) = MMPSingularities.polynomial_frobenius_generator(p,Δ₁fpminus1*a)
#
#  m = N*(p-1)
#  critical_ind = MMPSingularities.index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)  start_vector = lift_to_Int64(vector(fpminus1,m))
#  (coefs,degs) = MMPSingularities.GPUPolynomials.convert_to_gpu_representation(Δ₁fpminus1)
#
#  #println("creating matrix...")
#  #=@time=# 
#  M_true = MMPSingularities.matrix_of_multiply_then_split_sortmodp_dict(p,coefs,degs,m)
#
#  M = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(p,coefs,degs,m)
#
#  @test M == M_true

end

function test_delta1_char_3()
  p = 3
  N = 4

  R, (x, y, z, w) = polynomial_ring(GF(3), 4)
  f = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2

  # CPU method
  
  fpminus1 = f^(p-1)
  DD1 = MMPSingularities.Δ₁l(p,fpminus1)
  gpu_rep = MMPSingularities.convert_to_gpu_representation(DD1)

  # GPU method
  
  isfsplit, fpminus1 = MMPSingularities.isFSplit2(p, f)
  fpminus1_gpu = MMPSingularities.convert_to_gpu_representation(fpminus1)
  fpminus1_homog = MMPSingularities.GPUPolynomials.HomogeneousPolynomial(fpminus1_gpu...)
  pregen = MMPSingularities.pregen_delta1(size(fpminus1_homog.degrees, 2),p)
  MMPSingularities.GPUPolynomials.sort_to_kronecker_order(fpminus1_homog, pregen.key1)
  Δ₁fpminus1 = MMPSingularities.delta1(fpminus1_homog,p;pregen)

  gpuDelta1 = zero(R)

  for (i, coeff) in enumerate(Δ₁fpminus1.coeffs)
    exp_row = Δ₁fpminus1.degrees[i, :]
    term = coeff * x^exp_row[1] * y^exp_row[2] * z^exp_row[3] * w^exp_row[4]
    gpuDelta1 += term
  end

  @test DD1 == gpuDelta1
#   cpu_degs = gpu_rep[2]
#   gpu_degs = Δ₁fpminus1.degrees

#   cpu_coeffs = gpu_rep[1]
#   gpu_coeffs = Δ₁fpminus1.coeffs

#   @test size(cpu_degs) == size(gpu_degs)
#   @test length(cpu_coeffs) == length(gpu_coeffs)

#   @test cpu_degs == gpu_degs
#   @test cpu_coeffs == gpu_coeffs 

end

function test_matrices_all()

  #test_k3_5_matrix()
  #TODO: known to fail
  #test_k3_3_matrix()
  cases_3, _ = test_cases_k3surf(3)
  f2 = cases_3[2]
  f3 = cases_3[3]
  #
  test_matrix_K3(3,f2,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker2_correct)
  test_matrix_K3(3,f3,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker2_correct)

  cases_2, _ = test_cases_k3surf(2)
  ff = cases_2[10]
  
  test_matrix_K3(2,ff,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker2_correct)

  
  #test_matrix_qfs_cy_char_2()
  
  #test_delta1_char_3()

  #test_time_K3_5()
end
