#include("../src/MMPSingularities.jl")
#include("../src/RandomPolynomials.jl")
#
#using Test
#using CUDA
#using Oscar

function height_run_tests()
    # height_test_K3_3()
    # height_test_K3_5()
    # height_test_K3_7()
    # height_test_K3_11()
    height_test_K3_11()
    MMPSingularities.reset_times()
    height_test_K3_11()
    MMPSingularities.display_times()
    # test_matrix()
end

function height_test_K3_3()
    n = 4
    p = 3

    R, (w,x,y,z) = polynomial_ring(GF(3),["w","x","y","z"])

    f1 = x^4 + y^4 + z^4 + 2w^4 + x^2*y*w + y*z^2*w
    f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
    f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w
    f4 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2
    f5 = x^4 + y^4 + z^4 + w^4 + x^3*z + z^3*w + y*z^2*w + y*z*w^2
    f6 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3
    f7 = x^4 + y^4 + z^4 + w^4 + x*y^2*z + x*z^2*w + y*z*w^2 + y^2*z*w
    f8 = x^4 + x^2*y*z + x^2*y*w + 2*x^2*z^2 + x*y*w^2 + 2*y^4 + y^3*w + z^4 + w^4
    f9 = x^4 + y^4 + z^4 + w^4 + x*y^3 + y^3*w + z^2*w^2 + 2*x*y*z^2 + y*z*w^2
    f10 = x^4 + 2*x^2*y*z + x^2*y*w + x*y^2*w + y^4 + y^3*w + y^2*z^2 + 2*y^2*z*w + y^2*w^2 + y*z^3 + y*z^2*w + y*z*w^2 + z^4 + z*w^3
    finfty = x^4 + y^4 + z^4 + w^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    println("Running K3_3 tests...")
    @test qfs_height_fn(f1) == 1
    @test qfs_height_fn(f2) == 2
    @test qfs_height_fn(f3) == 3
    @test qfs_height_fn(f4) == 4
    @test qfs_height_fn(f5) == 5
    @test qfs_height_fn(f6) == 6
    @test qfs_height_fn(f7) == 7
    @test qfs_height_fn(f8) == 8
    @test qfs_height_fn(f9) == 9
    @test qfs_height_fn(f10) == 10
    @test qfs_height_fn(finfty) > 10
end

function height_test_K3_5()
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

    finfty = x1^4 + x2^4 + x3^4 + x4^4 + x1*x2*x3*x4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)
    
    println("Running K3_5 tests...")
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
    @test qfs_height_fn(finfty) > 10
end

function height_test_K3_7()
    n = 4
    p = 7

    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    fone1 = 5*x1^4 + 5*x1^3*x3 + 2*x1^3*x4 + 3*x1^2*x2*x3 + x1^2*x2*x4 + 6*x1^2*x3^2 + 3*x1^2*x3*x4 +
    3*x1^2*x4^2 + 4*x1*x2^3 + 6*x1*x2^2*x3 + 2*x1*x2^2*x4 + 4*x1*x2*x3^2 + 5*x1*x2*x3*x4 + 4*x1*x2*x4^2
    + 5*x1*x3^3 + 4*x1*x3^2*x4 + 5*x1*x4^3 + 5*x2^4 + x2^3*x3 + 4*x2^3*x4 + 5*x2^2*x3^2 + x2^2*x3*x4 +
    x2*x3^3 + 2*x2*x3^2*x4 + 2*x2*x3*x4^2 + x2*x4^3 + 3*x3^4 + 5*x3^3*x4 + 3*x3^2*x4^2 + x4^4

    ftwo1 = 3*x1^4 + 4*x1^3*x2 + x1^3*x3 + x1^3*x4 + x1^2*x2^2 + 5*x1^2*x2*x3 + 5*x1^2*x2*x4 + 
    3*x1^2*x3^2 + 5*x1^2*x3*x4 + 6*x1^2*x4^2 + 2*x1*x2^3 + x1*x2^2*x3 + 5*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 
    x1*x2*x4^2 + 2*x1*x3^3 + 3*x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 4*x2^3*x3 + 4*x2^3*x4 + 
    6*x2^2*x3^2 + 3*x2^2*x3*x4 + 3*x2*x3^3 + 4*x2*x3*x4^2 + 2*x3^4 + 4*x3^3*x4 + 4*x3^2*x4^2 + 2*x3*x4^3 + 6*x4^4

    fthree1 = 4*x1^4 + x1^3*x2 + 2*x1^3*x3 + 6*x1^3*x4 + 6*x1^2*x2^2 + 3*x1^2*x2*x3 +
    3*x1^2*x2*x4 + 2*x1^2*x3*x4 + 4*x1^2*x4^2 + 2*x1*x2^3 + 5*x1*x2^2*x4 + 5*x1*x2*x3^2 +
    4*x1*x2*x3*x4 + 4*x1*x2*x4^2 + 6*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 2*x1*x4^3 +
    3*x2^4 + 2*x2^3*x3 + 5*x2^2*x3^2 + 5*x2^2*x3*x4 + 3*x2^2*x4^2 + 4*x2*x3^3 +
    6*x2*x3^2*x4 + 5*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3 + 5*x4^4

    ffour1 = 2*x1^4 + 6*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x3 +
    3*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x3*x4 + 2*x1^2*x4^2 + 3*x1*x2^3 + 6*x1*x2^2*x4 +
    x1*x2*x3^2 + 6*x1*x2*x3*x4 + x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x3^2*x4 + 5*x1*x3*x4^2 +
    2*x1*x4^3 + 6*x2^4 + 3*x2^3*x3 + 5*x2^2*x3^2 + x2^2*x3*x4 + 5*x2^2*x4^2 + 4*x2*x3^3 +
    3*x2*x3^2*x4 + x2*x4^3 + 6*x3^4 + 2*x3^3*x4 + x3^2*x4^2 + 3*x3*x4^3 + 2*x4^4

    ffive1 = 5*x1^4 + 6*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x4 +
    2*x1^2*x3^2 + 3*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 2*x1*x2*x3^2 +
    3*x1*x2*x3*x4 + 5*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 6*x2^4 +
    5*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 3*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 +
    3*x2*x4^3 + 5*x3^4 + 6*x3^2*x4^2 + 6*x3*x4^3 + 3*x4^4

    fsix1 = x1^4 + x1^3*x2 + 4*x1^3*x3 + 6*x1^3*x4 + 6*x1^2*x2^2 + 2*x1^2*x2*x4 +
    6*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + 2*x1*x2^2*x4 + 2*x1*x2*x3^2 +
    5*x1*x2*x3*x4 + 6*x1*x2*x4^2 + 6*x1*x3^2*x4 + 3*x1*x3*x4^2 + 6*x2^4 + 2*x2^3*x3 +
    3*x2^3*x4 + 5*x2^2*x3^2 + 4*x2^2*x3*x4 + 6*x2^2*x4^2 + 5*x2*x3^2*x4 + x2*x3*x4^2 +
    3*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 5*x3^2*x4^2 + 2*x3*x4^3 + 4*x4^4

    fseven1 = 2*x1^3*x2 + 2*x1^3*x3 + 2*x1^3*x4 + x1^2*x2^2 + 2*x1^2*x2*x3 + 3*x1^2*x2*x4 +
    5*x1^2*x3^2 + 6*x1^2*x3*x4 + x1^2*x4^2 + 2*x1*x2^3 + 5*x1*x2^2*x3 + x1*x2*x3^2 +
    2*x1*x2*x3*x4 + 6*x1*x2*x4^2 + 4*x1*x3^3 + 6*x1*x3^2*x4 + 5*x1*x3*x4^2 + 2*x1*x4^3 +
    2*x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x4^2 + 3*x2*x3^3 + x2*x3^2*x4 +
    5*x2*x3*x4^2 + 5*x2*x4^3 + 5*x3^3*x4 + x3^2*x4^2 + 6*x3*x4^3 + 6*x4^4

    feight1 = 2*x1^3*x2 + 2*x1^3*x4 + 4*x1^2*x2^2 + 6*x1^2*x2*x3 + 5*x1^2*x2*x4 + 4*x1^2*x3^2 +
    3*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + 4*x1*x2*x3^2 +
    5*x1*x2*x3*x4 + x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 +
    5*x2^3*x3 + 5*x2^3*x4 + 6*x2^2*x3*x4 + 6*x2^2*x4^2 + 4*x2*x3^2*x4 + 3*x2*x3*x4^2 +
    2*x2*x4^3 + 6*x3^3*x4 + 6*x3^2*x4^2 + 4*x3*x4^3

    fnine1 = 2*x1^3*x2 + x1^3*x3 + 6*x1^3*x4 + 6*x1^2*x2^2 + 4*x1^2*x2*x3 + 2*x1^2*x2*x4 +
    3*x1^2*x3*x4 + x1^2*x4^2 + x1*x2^3 + x1*x2^2*x3 + 6*x1*x2^2*x4 + 6*x1*x2*x3^2 +
    6*x1*x2*x3*x4 + 6*x1*x2*x4^2 + 2*x1*x3^3 + 4*x1*x3*x4^2 + 6*x1*x4^3 + 6*x2^3*x3 +
    4*x2^3*x4 + 3*x2^2*x3^2 + 4*x2*x3^3 + 5*x2*x3^2*x4 + 4*x2*x3*x4^2 + 5*x2*x4^3 +
    3*x3^3*x4 + 4*x3^2*x4^2 + 2*x3*x4^3 + 3*x4^4

    ften1 = 3*x1^4 + 2*x1^3*x2 + x1^3*x3 + x1^3*x4 + 4*x1^2*x2*x3 + 2*x1^2*x2*x4 +
    5*x1^2*x3*x4 + 6*x1^2*x4^2 + x1*x2^3 + 2*x1*x2^2*x4 + 5*x1*x2*x3^2 + 3*x1*x2*x3*x4 +
    4*x1*x2*x4^2 + 5*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 6*x2^4 + x2^3*x4 +
    6*x2^2*x3^2 + x2^2*x3*x4 + 4*x2^2*x4^2 + x2*x3^3 + 5*x2*x4^3 + 2*x3^4 + 5*x3^3*x4 +
    5*x3^2*x4^2 + x3*x4^3 + 6*x4^4

    finfty = 3*x1^4 + 3*x1^3*x2 + 3*x1^3*x3 + 6*x1^2*x2^2 + 3*x1^2*x2*x4 + 2*x1^2*x3^2 + 2*x1^2*x3*x4 + 3*x1^2*x4^2 + 6*x1*x2^3 + 5*x1*x2^2*x3 + x1*x2*x3*x4 + 5*x1*x2*x4^2 + 5*x1*x3^3 + 4*x1*x3^2*x4 + 3*x1*x3*x4^2 + 6*x1*x4^3 + x2^4 + 4*x2^3*x4 + 3*x2^2*x3^2 + 5*x2^2*x3*x4 + 5*x2*x3^3 + x2*x3^2*x4 + 6*x2*x3*x4^2 + x3^3*x4 + x3^2*x4^2 + 3*x3*x4^3 + 4*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)
    
    println("Running K3_7 tests...")
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
    @test qfs_height_fn(finfty) > 10
end

function height_test_K3_11()
    n = 4
    p = 11

    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    fone1 = 4*x1^4 + 6*x1^3*x2 + x1^3*x3 + 2*x1^3*x4 + 3*x1^2*x2^2 + x1^2*x2*x3 + 3*x1^2*x2*x4 + 6*x1^2*x3^2 + 6*x1^2*x3*x4 + 8*x1^2*x4^2 + 7*x1*x2^3 + 2*x1*x2^2*x3 + 8*x1*x2^2*x4 + 8*x1*x2*x3*x4 + 10*x1*x2*x4^2 + 10*x1*x3^3 + 9*x1*x3^2*x4 + 6*x1*x3*x4^2 + 3*x1*x4^3 + 6*x2^4 + 7*x2^3*x3 + 4*x2^3*x4 + 10*x2^2*x3^2 + 3*x2^2*x3*x4 + 5*x2^2*x4^2 + 4*x2*x3^2*x4 + 6*x2*x4^3 + 3*x3^4 + 4*x3^3*x4 + 7*x3^2*x4^2 + 9*x3*x4^3 + 5*x4^4

    ftwo1 = 4*x1^4 + 5*x1^3*x2 + 9*x1^3*x3 + 2*x1^3*x4 + 8*x1^2*x2^2 + x1^2*x2*x3 + 9*x1^2*x2*x4 + x1^2*x3^2 + 8*x1^2*x3*x4 + 6*x1*x2^3 + 10*x1*x2^2*x3 + 2*x1*x2^2*x4 + 10*x1*x2*x3^2 + 9*x1*x2*x3*x4 + 6*x1*x2*x4^2 + 8*x1*x3^3 + 4*x1*x3^2*x4 + 7*x1*x3*x4^2 + 9*x1*x4^3 + 3*x2^4 + 7*x2^3*x3 + 6*x2^3*x4 + 10*x2^2*x3^2 + 8*x2^2*x3*x4 + x2^2*x4^2 + 9*x2*x3^3 + 6*x2*x3^2*x4 + x2*x3*x4^2 + 9*x3^4 + 10*x3^3*x4 + x3^2*x4^2 + x3*x4^3 + 4*x4^4

    # This example is nondegenerate in the sense of Costa, Harvey, and Kedlaya
    ftwo2 = 2*x1^4 + x1^3*x2 + x1^3*x3 + 7*x1^3*x4 + 10*x1^2*x2^2 + 7*x1^2*x2*x3 + 6*x1^2*x2*x4 + 5*x1^2*x3^2 + 9*x1^2*x3*x4 + 9*x1*x2^3 + 3*x1*x2^2*x3 + 3*x1*x2^2*x4 + 6*x1*x2*x3^2 + 5*x1*x2*x3*x4 + 7*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 7*x1*x3*x4^2 + 6*x1*x4^3 + 7*x2^4 + 2*x2^3*x4 + 3*x2^2*x3^2 + 10*x2^2*x3*x4 + x2^2*x4^2 + 4*x2*x3^3 + 3*x2*x3^2*x4 + 10*x2*x3*x4^2 + 2*x2*x4^3 + 4*x3^4 + 8*x3^3*x4 + 9*x3*x4^3 + 8*x4^4

    fthree1 = 10*x1^4 + 9*x1^3*x2 + 5*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 9*x1^2*x2*x3 + 4*x1^2*x2*x4 + 10*x1^2*x3^2 + 4*x1^2*x3*x4 + 8*x1^2*x4^2 + 8*x1*x2^3 + 9*x1*x2^2*x3 + 3*x1*x2^2*x4 + 7*x1*x2*x3^2 + 3*x1*x2*x4^2 + 8*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 7*x1*x4^3 + 2*x2^4 + 3*x2^3*x4 + x2^2*x3^2 + x2^2*x3*x4 + x2^2*x4^2 + 5*x2*x3^3 + 9*x2*x3^2*x4 + 9*x2*x3*x4^2 + 4*x2*x4^3 + 5*x3^4 + 10*x3^3*x4 + 10*x3*x4^3 + 10*x4^4

    ffour1 = 2*x1^4 + 4*x1^3*x2 + 9*x1^3*x3 + 10*x1^3*x4 + 2*x1^2*x2^2 + 4*x1^2*x2*x3 + 4*x1^2*x2*x4 + 4*x1^2*x3^2 + 10*x1^2*x3*x4 + 9*x1^2*x4^2 + 5*x1*x2^3 + 5*x1*x2^2*x3 + x1*x2^2*x4 + 8*x1*x2*x3^2 + 2*x1*x2*x3*x4 + 10*x1*x2*x4^2 + 8*x1*x3^3 + 7*x1*x3^2*x4 + 5*x1*x3*x4^2 + 4*x1*x4^3 + 3*x2^4 + 6*x2^3*x3 + 4*x2^3*x4 + 10*x2^2*x3^2 + 5*x2^2*x3*x4 + 5*x2^2*x4^2 + x2*x3^3 + 5*x2*x4^3 + 5*x3^4 + 7*x3^2*x4^2 + 5*x3*x4^3 + 9*x4^4

    ffive1 = 10*x1^4 + x1^3*x2 + 6*x1^3*x3 + 3*x1^3*x4 + x1^2*x2^2 + 9*x1^2*x2*x3 + 6*x1^2*x2*x4 + 6*x1^2*x3^2 + 8*x1^2*x3*x4 + 4*x1^2*x4^2 + 3*x1*x2^3 + 7*x1*x2^2*x3 + 3*x1*x2^2*x4 + 7*x1*x2*x3^2 + 9*x1*x2*x3*x4 + 8*x1*x2*x4^2 + 7*x1*x3^3 + x1*x3*x4^2 + 7*x1*x4^3 + x2^4 + 3*x2^3*x3 + 7*x2^3*x4 + 5*x2^2*x3^2 + 7*x2^2*x3*x4 + 8*x2^2*x4^2 + 8*x2*x3^3 + 5*x2*x3^2*x4 + x2*x3*x4^2 + 9*x2*x4^3 + 7*x3^4 + 4*x3^3*x4 + 4*x3^2*x4^2 + 3*x3*x4^3

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)
    
    println("Running K3_11 tests...")
    @test qfs_height_fn(fone1) == 1
    @test qfs_height_fn(ftwo1) == 2
    @test qfs_height_fn(ftwo2) == 2
    @test qfs_height_fn(fthree1) == 3
    @test qfs_height_fn(ffour1) == 4
    @test qfs_height_fn(ffive1) == 5
end

function height_test_K3_13()

# fthree1 = 8*x1^4 + 2*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + 6*x1^2*x2^2 + 7*x1^2*x2*x3 + 5*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x4^2 + 11*x1*x2^3 + 10*x1*x2^2*x3 + 3*x1*x2^2*x4 + 5*x1*x2*x3^2 + 10*x1*x2*x3*x4 + 7*x1*x2*x4^2 + 12*x1*x3^3 + 12*x1*x3^2*x4 + 5*x1*x3*x4^2 + 7*x1*x4^3 + 7*x2^4 + 6*x2^3*x3 + 3*x2^3*x4 + 10*x2^2*x3^2 + 5*x2^2*x3*x4 + 12*x2^2*x4^2 + x2*x3^3 + 3*x2*x3^2*x4 + 12*x2*x3*x4^2 + 8*x2*x4^3 + 10*x3^4 + 7*x3^3*x4 + 4*x3^2*x4^2 + 8*x3*x4^3 + 2*x4^4function height_test_K3_13()
    n = 4
    p = 13

    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    fone1 = 6*x1^4 + 7*x1^3*x3 + 4*x1^3*x4 + 6*x1^2*x2^2 + 7*x1^2*x2*x3 + 9*x1^2*x2*x4 + 2*x1^2*x3^2 + 3*x1^2*x3*x4 + 12*x1^2*x4^2 + 8*x1*x2^3 + 4*x1*x2^2*x3 + x1*x2^2*x4 + 9*x1*x2*x3^2 + 8*x1*x2*x3*x4 + 10*x1*x2*x4^2 + 8*x1*x3^3 + 2*x1*x3^2*x4 + 9*x1*x3*x4^2 + 4*x1*x4^3 + 5*x2^4 + 4*x2^3*x3 + 2*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 10*x2*x3^3 + 2*x2*x3^2*x4 + 2*x2*x3*x4^2 + 5*x2*x4^3 + 4*x3^4 + 3*x3^2*x4^2 + 2*x4^4

    ftwo1 = 5*x1^4 + 6*x1^3*x2 + 12*x1^3*x3 + 11*x1^3*x4 + 7*x1^2*x2^2 + 8*x1^2*x2*x3 + 3*x1^2*x2*x4 + 11*x1^2*x3^2 + 9*x1^2*x3*x4 + 6*x1^2*x4^2 + 3*x1*x2^3 + 2*x1*x2^2*x3 + x1*x2^2*x4 + 9*x1*x2*x4^2 + 6*x1*x3^3 + 3*x1*x3^2*x4 + 8*x1*x3*x4^2 + 7*x1*x4^3 + 6*x2^4 + 4*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 7*x2^2*x4^2 + 12*x2*x3^3 + 3*x2*x3^2*x4 + 2*x2*x3*x4^2 + x2*x4^3 + 6*x3^4 + 4*x3^3*x4 + 9*x3^2*x4^2 + 5*x3*x4^3 + 8*x4^4 

    fthree1 = 8*x1^4 + 2*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + 6*x1^2*x2^2 + 7*x1^2*x2*x3 + 5*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x4^2 + 11*x1*x2^3 + 10*x1*x2^2*x3 + 3*x1*x2^2*x4 + 5*x1*x2*x3^2 + 10*x1*x2*x3*x4 + 7*x1*x2*x4^2 + 12*x1*x3^3 + 12*x1*x3^2*x4 + 5*x1*x3*x4^2 + 7*x1*x4^3 + 7*x2^4 + 6*x2^3*x3 + 3*x2^3*x4 + 10*x2^2*x3^2 + 5*x2^2*x3*x4 + 12*x2^2*x4^2 + x2*x3^3 + 3*x2*x3^2*x4 + 12*x2*x3*x4^2 + 8*x2*x4^3 + 10*x3^4 + 7*x3^3*x4 + 4*x3^2*x4^2 + 8*x3*x4^3 + 2*x4^4

    ffour1 = 4*x1^4 + 4*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 9*x1^2*x2^2 + 6*x1^2*x2*x3 + 7*x1^2*x2*x4 + 10*x1^2*x3^2 + x1^2*x3*x4 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 12*x1*x2*x3^2 + 7*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 11*x1*x3^3 + 9*x1*x3^2*x4 + 10*x1*x3*x4^2 + 11*x1*x4^3 + 3*x2^4 + 5*x2^3*x3 + 8*x2^3*x4 + 5*x2^2*x3*x4 + 5*x2^2*x4^2 + 5*x2*x3^3 + 10*x2*x3^2*x4 + 2*x2*x3*x4^2 + 10*x2*x4^3 + 4*x3^4 + 5*x3^2*x4^2 + 4*x3*x4^3 + 6*x4^4

    ffive1 = 11*x1^4 + 4*x1^3*x2 + 12*x1^3*x3 + 4*x1^3*x4 + 6*x1^2*x2^2 + 10*x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 7*x1^2*x3*x4 + 4*x1^2*x4^2 + 6*x1*x2^3 + 11*x1*x2^2*x3 + 7*x1*x2^2*x4 + 8*x1*x2*x3^2 + 10*x1*x2*x4^2 + x1*x3^3 + 9*x1*x3^2*x4 + 8*x1*x3*x4^2 + 11*x1*x4^3 + 4*x2^4 + 8*x2^3*x3 + 5*x2^2*x3*x4 + 7*x2^2*x4^2 + 8*x2*x3^3 + 6*x2*x3^2*x4 + 5*x2*x4^3 + 2*x3^4 + 10*x3^3*x4 + 8*x3^2*x4^2 + 10*x3*x4^3 + 6*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)
    
    println("Running K3_13 tests...")
    @test qfs_height_fn(fone1) == 1
    @test qfs_height_fn(ftwo1) == 2
    @test qfs_height_fn(fthree1) == 3
    @test qfs_height_fn(ffour1) == 4
    @test qfs_height_fn(ffive1) == 5
end


function test_matrix()
    n = 4
    p = 5
    
    R, vars = polynomial_ring(GF(p), n)
    (x1, x2, x3, x4) = vars

    Δ₁plan = MMPSingularities.plan_Δ₁(n, p)

    f = random_homog_poly_mod(p, vars, n)

    momtspregen = MMPSingularities.pregen_MOMTS(n, p)

    fpminus1 = MMPSingularities.CufpMPolyRingElem((f ^ (p - 1)).data)
    fpminus1.opPlan = Δ₁plan
    Δ₁fpminus1 = MMPSingularities.Δ₁l(fpminus1)

    # display(@benchmark mat1 = MMPSingularities.matrix_of_multiply_then_split($Δ₁fpminus1; plan = $momtspregen, alg = 1))
    # display(@benchmark CUDA.@sync mat2 = MMPSingularities.matrix_of_multiply_then_split($Δ₁fpminus1; plan = $momtspregen, alg = 2))
    # display(@benchmark mat3 = MMPSingularities.matrix_of_multiply_then_split($Δ₁fpminus1; plan = $momtspregen, alg = 3))
    # display(@benchmark mat4 = MMPSingularities.matrix_of_multiply_then_split($Δ₁fpminus1; plan = $momtspregen, alg = 4))
    # display(@benchmark CUDA.@sync mat5 = MMPSingularities.matrix_of_multiply_then_split($Δ₁fpminus1; plan = $momtspregen, alg = 5))

    mat1 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1; plan = momtspregen, alg = 1)
    CUDA.@sync mat2 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1; plan = momtspregen, alg = 2)
    mat3 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1; plan = momtspregen, alg = 3)
    mat4 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1; plan = momtspregen, alg = 4)
    CUDA.@sync mat5 = MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1; plan = momtspregen, alg = 5)


    @assert mat2 == mat5 string(f)
    @assert mat1 == mat3 string(f)
    @assert mat1 == mat4 string(f)
    @assert mat1 == Array(mat2) string(f)

    # end
end

