include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")

using Test
using CUDA
using Oscar

function run_tests()
    # test_height()
    test_K3_5()
    # time_K3_7()
    # test_matrix()
end

function test_height()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)
    # (x, y, z, w) = vars
    # f = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3 + w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4

    (x, y, z, w) = vars
    f = x^4 + y^4 + z^4 + w^4
    pregen = MMPSingularities.pregen_qfsheight(n, p)

    height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, f, 10, pregen)
    # realheight = MMPSingularities.quasiFSplitHeight_CY_lift(p, f, 10)
    println("height: $height")
    # println("realheight: $realheight")
    # @test height == realheight
end

function test_K3_5()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    fone1 = x1^4 + x2^4 + x3^4 + x4^4

    ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4

    fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4

    ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

    ffive1 = 4*x1^2*x2^2 + 4*x1^2*x2*x3 + 2*x1^2*x2*x4 + 3*x1^2*x3*x4 + 4*x1*x2^2*x3 + 3*x1*x2*x3*x4 + x1*x3^3 + x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + x2^4 + 2*x2*x3^2*x4 + 4*x2*x3*x4^2

    fsix1 = 4*x1^3*x2 + 4*x1^3*x4 + 4*x2^4 + 2*x2*x3^2*x4 + x2*x3*x4^2 + x3^3*x4

    fseven1 = 4*x1^3*x4 + x1*x2^3 + 2*x1*x3^3 + 2*x2^3*x4 + x2^2*x4^2 + x3^3*x4 + 3*x3^2*x4^2 + 2*x3*x4^3

    feight1 = 2*x1^3*x3 + x1^2*x2^2 + x1^2*x4^2 + x1*x2^2*x4 + 3*x1*x3^2*x4 + x2^4 + 2*x2^3*x3 + 3*x2^2*x4^2 + 4*x2*x3^2*x4 + 3*x3*x4^3

    fnine1 = 3*x1^4 + 3*x1^3*x2 + 3*x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + 3*x1^2*x2*x4 + 3*x1^2*x3^2 + 2*x1^2*x3*x4 + 2*x1^2*x4^2 + 4*x1*x2^3 + 2*x1*x2^2*x3 + 4*x1*x2*x3^2 + 2*x1*x2*x3*x4 + 4*x1*x2*x4^2 + x1*x3^3 + 3*x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + 3*x2^3*x3 + 4*x2^3*x4 + 3*x2^2*x3*x4 + x2^2*x4^2 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 4*x2*x4^3 + 3*x3*x4^3 + 4*x4^4

    ften1 = 2*x1^4 + 4*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x2*x4 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + 2*x1^2*x4^2 + x1*x2^3 + 4*x1*x2^2*x4 + 3*x1*x2*x3^2 + 3*x1*x2*x4^2 + 2*x1*x3^3 + 3*x1*x3^2*x4 + 2*x1*x3*x4^2 + x1*x4^3 + 3*x2^4 + 2*x2^3*x3 + 2*x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 3*x2^2*x4^2 + x2*x3^3 + 2*x2*x3*x4^2 + 2*x2*x4^3 + 4*x3^4 + x3^3*x4 + 3*x3^2*x4^2 + 4*x3*x4^3 + 3*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, x, 10, pregen)
    
    @test qfs_height_fn(fone1) == 1
    @test qfs_height_fn(ftwo1) == 2
    @test qfs_height_fn(fthree1) == 3
    @test qfs_height_fn(ffour1) == 4
    @test qfs_height_fn(ffive1) == 5
    @test qfs_height_fn(fsix1) == 6
    @test qfs_height_fn(fseven1) == 7
    @test qfs_height_fn(feight1) == 8
    @test qfs_height_fn(fnine1) == 9
    @test qfs_height_fn(ften1) == 10
end

function time_K3_7()
    n = 4
    p = 7

    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    f5 = 5*x1^4 + 6*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x4 + 2*x1^2*x3^2 + 3*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + 5*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 6*x2^4 + 5*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 3*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 5*x3^4 + 6*x3^2*x4^2 + 6*x3*x4^3 + 3*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)
    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, x, 10, pregen)

    @assert qfs_height_fn(f5) == 5

    for i in 1:10
        qfs_height_fn(f5)
    end
end

function test_matrix()
    n = 4
    p = 5
    
    R, vars = polynomial_ring(GF(p), n)
    (x1, x2, x3, x4) = vars

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
        fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
        Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly

        # @time mat0 = MMPSingularities.matrix_of_multiply_then_split_correct(Δ₁fpminus1)
        @time mat1 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
        CUDA.@time mat2 = MMPSingularities.matrix_of_multiply_then_split_gpu(Δ₁fpminus1, momtspregen)
        @time mat3 = MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
        @time mat4 = MMPSingularities.matrix_of_multiply_then_split_wics(Δ₁fpminus1)
        CUDA.@time mat5 = MMPSingularities.matrix_of_multiply_then_split_wics_gpu(Δ₁fpminus1, momtspregen)
        
        println()

        # @assert mat0 == mat1
        @assert mat2 == mat5 string(f)
        @assert mat1 == mat3 string(f)
        @assert mat1 == mat4 string(f)
        @assert mat1 == Array(mat2) string(f)

    end
    # end
end

run_tests()