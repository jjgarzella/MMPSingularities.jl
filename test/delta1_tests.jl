include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")

using Test
using CUDA
using Oscar

function oscar_delta1(poly, p)

    R = parent(poly)
  
    originallift = map_coefficients(x -> lift(ZZ,x),poly)
  
    ZR = parent(originallift)
  
    nocrossterms = sum(terms(originallift) .^p)
    withcrossterms = originallift^p
    crossterms = withcrossterms - nocrossterms
    Δlift = map_coefficients(x -> divexact(x,p),crossterms)
  
    change_coefficient_ring(coefficient_ring(R),Δlift,parent=R)
end

function run_tests()
    test_K3_5()
    test_K3_7()
    time_K3_5()
    time_K3_7()
    # test_memory_safe()
end

function test_K3_5()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = oscar_delta1(fpminus1, p)

    pregen = MMPSingularities.pregen_delta1(n, p)

    fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)

    gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)

    gpu_result = MMPSingularities.convert(FqMPolyRingElem, gpud1)

    if gpu_result != oscar_result
        println("$f failed for test_K3_5!")
    end
    @test gpu_result == oscar_result
end

function time_K3_5()
    n = 4
    p = 5
    pregen = MMPSingularities.pregen_delta1(n, p)

    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)
    fpminus1 = f ^ (p - 1)
    fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)
    gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)
    println("Quartic K3_5 times: ")
    for i in 1:10
        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = f ^ (p - 1)
        fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)
        CUDA.@time gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)
    end
end

function test_K3_7()
    n = 4
    p = 7
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = oscar_delta1(fpminus1, p)

    pregen = MMPSingularities.pregen_delta1(n, p)

    fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)

    gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)

    gpu_result = MMPSingularities.convert(FqMPolyRingElem, gpud1)

    if gpu_result != oscar_result
        println("$f failed for test_K3_7!")
    end
    @test gpu_result == oscar_result
end

function time_K3_7()
    n = 4
    p = 7
    pregen = MMPSingularities.pregen_delta1(n, p)

    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)
    fpminus1 = f ^ (p - 1)
    fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)
    gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)
    println("Quartic K3_7 times: ")
    for i in 1:10
        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = f ^ (p - 1)
        fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)
        CUDA.@time gpud1 = MMPSingularities.delta1(fpminus1_gpu, p; pregen = pregen)
    end
end

function test_memory_safe()
    n = 4
    p = 7
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)
    fpminus1 = f ^ (p - 1)

    pregen = MMPSingularities.pregen_delta1(n, p)

    fpminus1_gpu = MMPSingularities.HomogeneousPolynomial(fpminus1)

    unsafed1 = MMPSingularities.memoryunsafe_delta1(fpminus1_gpu, p; pregen = pregen)
    pregen.gpupregen.nttpregen.butterfly = Array(MMPSingularities.generate_butterfly_permutations(length(pregen.gpupregen.nttpregen.butterfly)))
    safed1 = MMPSingularities.memorysafe_delta1(fpminus1_gpu, p; pregen = pregen)

    # @test safed1 == unsafed1
    if safed1.poly != unsafed1.poly
        println("$f failed for test_memory_safe!")
    end
    @test safed1.poly == unsafed1.poly
end

run_tests()
