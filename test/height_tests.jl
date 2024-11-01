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
    p = 7

    R, vars = polynomial_ring(GF(p), n)
    # (x, y, z, w) = vars
    # f = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3 + w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4

    (x1, x2, x3, x4) = vars
    f = 2*x1^3*x2 + 2*x1^3*x3 + 3*x1^2*x2^2 + 5*x1^2*x2*x4 + 5*x1^2*x3^2 + 2*x1^2*x3*x4 + 5*x1^2*x4^2 + 2*x1*x2^3 + x1*x2^2*x3 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 5*x1*x3^2*x4 + 5*x1*x3*x4^2 + x1*x4^3 + 5*x2^3*x3 + 4*x2^3*x4 + 5*x2^2*x3^2 + 4*x2^2*x3*x4 + 4*x2^2*x4^2 + x2*x3^3 + 4*x2*x3^2*x4 + 5*x2*x3*x4^2 + 6*x2*x4^3 + x3*x4^3
    pregen = MMPSingularities.pregen_qfsheight(n, p, true)

    height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, f, 10, pregen)
    # realheight = MMPSingularities.quasiFSplitHeight_CY_lift(p, f, 10)
    println("height: $height")
    # println("realheight: $realheight")
    # @test height == realheight
end

function test_matrix()
    n = 4
    p = 7
    
    R, vars = polynomial_ring(GF(p), n)

    pregen = MMPSingularities.pregen_delta1(n, p)

    # fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
    # Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly
    # for i in 1:10
    f = random_homog_poly_mod(p, vars, n)
    # fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
    # Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly
    # Δ₁fpminus1 = MMPSingularities.Δ₁l(p, f ^ (p - 1))
    momtspregen = MMPSingularities.pregen_MOMTS(n, p)
    for i in 1:10
        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
        Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly

        mat3 = MMPSingularities.matrix_of_multiply_then_split_alex_gpu(Δ₁fpminus1, momtspregen)
        mat2 = MMPSingularities.matrix_of_multiply_then_split_alex(Δ₁fpminus1)
        mat1 = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
        println()
        @test Array(mat3) == mat2
        @test mat1 == mat2
    end
    # end
end

function istwo(x)
    return x == typeof(x)(2)
end

run_tests()