include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")

using Test
using CUDA
using Oscar

function run_tests()
    test_height()
    # test_matrix()
end

function test_height()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)
    # (x, y, z, w) = vars
    # f = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3 + w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4

    (x1, x2, x3, x4) = vars

    fnine1 = 3*x1^4 + 3*x1^3*x2 + 3*x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + 3*x1^2*x2*x4 + 3*x1^2*x3^2 + 2*x1^2*x3*x4 + 2*x1^2*x4^2 + 4*x1*x2^3 + 2*x1*x2^2*x3 + 4*x1*x2*x3^2 + 2*x1*x2*x3*x4 + 4*x1*x2*x4^2 + x1*x3^3 + 3*x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + 3*x2^3*x3 + 4*x2^3*x4 + 3*x2^2*x3*x4 + x2^2*x4^2 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 4*x2*x4^3 + 3*x3*x4^3 + 4*x4^4

    pregen = MMPSingularities.pregen_delta1(n, p)

    # @test 9 == MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, f, 10)
    # @test 9 == MMPSingularities.quasiFSplitHeight_CY_lift_matrix(p, fnine1, 10)
    @test 9 == MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, fnine1, 10, pregen)
end

function convert_to_gpu_representation(p)
    coeffs = coefficients(p)

  # julia by default doesn't realize that "ZZ" is not
  # an array, so insert it as a one-element tuple "(ZZ,)"
  # so that julia will know not to broadcast along it.
  coeffs_as_int64arr = UInt.(lift.((ZZ,),coeffs))

  exp_vecs = leading_exponent_vector.(terms(p))

  # shamelessly taken from 
  # https://discourse.julialang.org/t/how-to-convert-vector-of-vectors-to-matrix/72609/2 
  exponent_mat = reduce(hcat,exp_vecs)

  (coeffs_as_int64arr,exponent_mat)
end

function test_matrix()
    n = 4
    p = 5
    
    R, vars = polynomial_ring(GF(p), n)

    (x, y, z, w) = vars
    f = random_homog_poly_mod(p, vars, n)

    pregen = MMPSingularities.pregen_delta1(n, p)

    fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
    Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly

    @time mat2 = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
    @time mat3 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
    @test mat2 == mat3
    # @test mat1 == mat2
end

function istwo(x)
    return x == typeof(x)(2)
end

run_tests()