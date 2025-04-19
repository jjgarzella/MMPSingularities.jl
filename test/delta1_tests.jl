#include("../src/MMPSingularities.jl")
#include("../src/RandomPolynomials.jl")

#using Test
#using CUDA
#using Oscar

function delta1_run_tests()
    delta1_test_K3_2()
    delta1_test_K3_3()
    delta1_test_K3_5()
    delta1_test_K3_7()
    # delta1_time_K3_5()
    # delta1_time_K3_7()
    # delta1_time_K3_11()
    # delta1_time_K3_13()
end

function delta1_test_K3_2()
    n = 4
    p = 2
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = MMPSingularities.Δ₁l(fpminus1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    gpu_result = fpMPolyRingElem(gpud1)

    # Need to compare strings because == doesn't work, think it's because I have
    # zero terms maybe, or different allocation length, idk
    @test string(gpu_result) == string(oscar_result)
end

function delta1_test_K3_3()
    n = 4
    p = 3
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = MMPSingularities.Δ₁l(fpminus1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    gpu_result = fpMPolyRingElem(gpud1)

    @test string(gpu_result) == string(oscar_result)
end

function delta1_test_K3_5()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = MMPSingularities.Δ₁l(fpminus1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    gpu_result = fpMPolyRingElem(gpud1)

    @test string(gpu_result) == string(oscar_result)
    # @test gpu_result == oscar_result
end

function delta1_time_K3_5()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan
    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)
    println("Quartic K3_5 times: ")
    for i in 1:100
        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = f ^ (p - 1)
        fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
        fpminus1_gpu.opPlan = plan
        CUDA.@time gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)
        # gpud1 = MMPSingularities.Δ₁(fpminus1_gpu)
    end
end

function delta1_test_K3_7()
    n = 4
    p = 7
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    fpminus1 = f ^ (p - 1)
    oscar_result = MMPSingularities.Δ₁l(fpminus1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    gpu_result = fpMPolyRingElem(gpud1)

    @test string(gpu_result) == string(oscar_result)
end

function delta1_time_K3_7()
    n = 4
    p = 7

    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁(n, p)
    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan
    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)
    println("Quartic K3_7 times: ")
    for i in 1:50
        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = f ^ (p - 1)
        fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
        fpminus1_gpu.opPlan = plan
        CUDA.@time gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)
        # gpud1 = MMPSingularities.Δ₁(fpminus1_gpu)
    end
end

function delta1_time_K3_11()
    n = 4
    p = 11
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)
    # f = x^4 + y^4

    fpminus1 = f ^ (p - 1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    CUDA.@time gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    return
end

function delta1_time_K3_13()
    n = 4
    p = 13
    R, vars = polynomial_ring(GF(p), n)
    f = random_homog_poly_mod(p, vars, n)
    # f = x^4 + y^4

    fpminus1 = f ^ (p - 1)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    gpud1 = MMPSingularities.Δ₁l(fpminus1_gpu)

    return
end


