using MMPSingularities
using CUDA
using BenchmarkTools
using Oscar

function run()
    # delta1_time_K3_5()
    # delta1_time_K3_7()
    # delta1_p²_time_K3_5()
    # delta1_p²_time_K3_7()
    # delta1_p²_time_K3_11()
    delta1_p²_time_K3_13()
end

function delta1_time_K3_5()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_5 time: ")
    @btime CUDA.@sync MMPSingularities.Δ₁l($fpminus1_gpu)

    GC.gc()
    return
end

function delta1_time_K3_7()
    n = 4
    p = 7

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_7 time: ")
    @btime CUDA.@sync MMPSingularities.Δ₁l($fpminus1_gpu)

    GC.gc()
    return
end

function delta1_p²_time_K3_5()
    n = 4
    p = 5

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁lp²(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_5 time with lifting to ZZ/p^2: ")
    @btime CUDA.@sync MMPSingularities.Δ₁lp²($fpminus1_gpu)

    GC.gc()
    return
end

function delta1_p²_time_K3_7()
    n = 4
    p = 7

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁lp²(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_7 time with lifting to ZZ/p^2: ")
    @btime CUDA.@sync MMPSingularities.Δ₁lp²($fpminus1_gpu)

    GC.gc()
    return
end

function delta1_p²_time_K3_11()
    n = 4
    p = 11

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁lp²(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_11 time with lifting to ZZ/p^2: ")
    @btime CUDA.@sync MMPSingularities.Δ₁lp²($fpminus1_gpu)

    GC.gc()
    return
end

function delta1_p²_time_K3_13()
    n = 4
    p = 13

    R, vars = polynomial_ring(GF(p), n)
    f = MMPSingularities.random_homog_poly_mod(p, vars, n)

    plan = MMPSingularities.plan_Δ₁lp²(n, p)

    fpminus1 = f ^ (p - 1)

    fpminus1_gpu = MMPSingularities.CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = plan

    println("Quartic K3_13 time with lifting to ZZ/p^2: ")
    @btime CUDA.@sync MMPSingularities.Δ₁lp²($fpminus1_gpu)

    # for i in 1:10
    #     CUDA.@time MMPSingularities.Δ₁lp²(fpminus1_gpu)
    #     GC.gc()
    # end 
    GC.gc()
    return
end


run()