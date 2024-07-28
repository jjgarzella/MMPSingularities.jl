include("MMPSingularities.jl")

using Oscar
using CUDA
using Test
using Profile
using Dates

function run_test(time = Second(30))
    println("Creating pregen...")
    pregen = MMPSingularities.pregen_delta1(4, 5)
    # R, (x1, x2, x3, x4) = polynomial_ring(GF(5),4)
    R, (vars) = polynomial_ring(GF(5), 4)
    # poly = 4*x1^2*x2^2 + 4*x1*x2^2*x3 + 2*x1*x4^3 + 3*x2^3*x4
    startTime = now()
    testspassed = 0
    println("Starting test")
    while now() - startTime < time
        poly = MMPSingularities.random_homog_poly_mod(5, vars, 4)
        # poly = 2*x1^3*x2 + x1^2*x2*x3 + 2*x1^2*x4^2 + x2^3*x4
        realHeight = MMPSingularities.quasiFSplitHeight_CY_gpu(5, poly, 10, pregen)
        fakeHeight = -1
        try
            fakeHeight = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(5, poly, 10, pregen) 
        catch (e)
            println("$poly")
            throw(e)
        end
        @test realHeight == fakeHeight
        println("realHeight: $realHeight, fakeHeight: $fakeHeight")
        testspassed += 1
    end
    println("tests passed: $testspassed")
end

function run_specific_test()
    println("Creating pregen...")
    pregen = MMPSingularities.pregen_delta1(4, 5)
    R, (x1, x2, x3, x4) = polynomial_ring(GF(5), 4)
    poly = 3*x1^4 + 4*x1^3*x2 + 4*x1^3*x3 + x1^3*x4 + 2*x1^2*x2^2 + x1^2*x2*x4 + 2*x1^2*x3^2 + 4*x1^2*x4^2 + x1*x2^2*x3 + x1*x2^2*x4 + 3*x1*x2*x3^2 + 3*x1*x2*x3*x4 + 2*x1*x2*x4^2 + x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 4*x2*x3*x4^2 + 2*x2*x4^3 + x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3 + 3*x4^4

    poly2 = 4*x1^2*x2^2 + 4*x1*x2^2*x3 + 2*x1*x4^3 + 3*x2^3*x4

    realHeight = MMPSingularities.quasiFSplitHeight_CY_gpu(5, poly, 10, pregen)
    fakeHeight = -1
    try
        fakeHeight = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(5, poly, 10, pregen) 
    catch (e)
        println("$poly")
        throw(e)
    end
    @assert realHeight == fakeHeight

    realHeight = MMPSingularities.quasiFSplitHeight_CY_gpu(5, poly2, 10, pregen)
    fakeHeight = -1
    try
        fakeHeight = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(5, poly2, 10, pregen) 
    catch (e)
        println("$poly2")
        throw(e)
    end
    @assert realHeight == fakeHeight
end