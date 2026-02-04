using MMPSingularities
using CUDA
using Oscar

function run_benchmarks()
    bench_oscar_K3_3()
    bench_oscar_K3_5()
    bench_oscar_K3_7()

    bench_cpu_K3_3()
    bench_cpu_K3_5()
    bench_cpu_K3_7()

    bench_gpu_K3_3()
    bench_gpu_K3_5()
    bench_gpu_K3_7()
    bench_gpu_K3_11()
    bench_gpu_K3_13()
end

hrule() = println("-------------------------------")

function bench_oscar_K3_3()
    GC.gc()
    n = 4
    p = 3
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 2*x1^4 + 2*x1^3*x2 + x1^3*x3 + 2*x1^2*x2^2 + x1^2*x2*x3 + 2*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1*x2^3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 2*x1*x2*x3*x4 + x1*x2*x4^2 + x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 2*x2^4 + x2^3*x3 + x2^2*x3^2 + x2*x3^2*x4 + 2*x3^4 + x3^3*x4 + x3^2*x4^2 + 2*x4^4

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift(p, x, 10)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_3 Oscar Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()

    GC.gc()
    return
end

function bench_oscar_K3_5()
    GC.gc()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 3*x1^4 + x1^3*x4 + x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x4^2 + 2*x1*x2^3 + x1*x2^2*x3 + 3*x1*x2^2*x4 + 4*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 2*x1*x2*x4^2 + 4*x1*x3^3 + 3*x1*x3^2*x4 + x1*x3*x4^2 + 4*x1*x4^3 + 4*x2^3*x3 + 2*x2^3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 2*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3^2*x4^2 + 4*x3*x4^3 + 2*x4^4

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift(p, x, 10)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_5 Oscar Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()

    GC.gc()
    return
end

function bench_oscar_K3_7()
    GC.gc()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 5*x1^4 + 6*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x4 +
    2*x1^2*x3^2 + 3*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 2*x1*x2*x3^2 +
    3*x1*x2*x3*x4 + 5*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 6*x2^4 +
    5*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 3*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 +
    3*x2*x4^3 + 5*x3^4 + 6*x3^2*x4^2 + 6*x3*x4^3 + 3*x4^4

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift(p, x, 10)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:5
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_7 Oscar Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()

    GC.gc()
    return
end

function bench_cpu_K3_3()
    GC.gc()
    n = 4
    p = 3
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 2*x1^4 + 2*x1^3*x2 + x1^3*x3 + 2*x1^2*x2^2 + x1^2*x2*x3 + 2*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1*x2^3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 2*x1*x2*x3*x4 + x1*x2*x4^2 + x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 2*x2^4 + x2^3*x3 + x2^2*x3^2 + x2*x3^2*x4 + 2*x3^4 + x3^3*x4 + x3^2*x4^2 + 2*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_cpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_3 CPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()

    GC.gc()
    return
end

function bench_cpu_K3_5()
    GC.gc()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 3*x1^4 + x1^3*x4 + x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x4^2 + 2*x1*x2^3 + x1*x2^2*x3 + 3*x1*x2^2*x4 + 4*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 2*x1*x2*x4^2 + 4*x1*x3^3 + 3*x1*x3^2*x4 + x1*x3*x4^2 + 4*x1*x4^3 + 4*x2^3*x3 + 2*x2^3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 2*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3^2*x4^2 + 4*x3*x4^3 + 2*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_cpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_5 CPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()

    GC.gc()
    return
end

function bench_cpu_K3_7()
    GC.gc()
    n = 4
    p = 7
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 5*x1^4 + 6*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x4 +
    2*x1^2*x3^2 + 3*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 2*x1*x2*x3^2 +
    3*x1*x2*x3*x4 + 5*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 6*x2^4 +
    5*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 3*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 +
    3*x2*x4^3 + 5*x3^4 + 6*x3^2*x4^2 + 6*x3*x4^3 + 3*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_cpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:5
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_7 CPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    GC.gc()
    return
end

function bench_gpu_K3_3()
    GC.gc()
    n = 4
    p = 3
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 2*x1^4 + 2*x1^3*x2 + x1^3*x3 + 2*x1^2*x2^2 + x1^2*x2*x3 + 2*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1*x2^3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 2*x1*x2*x3*x4 + x1*x2*x4^2 + x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 2*x2^4 + x2^3*x3 + x2^2*x3^2 + x2*x3^2*x4 + 2*x3^4 + x3^3*x4 + x3^2*x4^2 + 2*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_3 GPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    GC.gc()
    return
end

function bench_gpu_K3_5()
    GC.gc()
    n = 4
    p = 5
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 3*x1^4 + x1^3*x4 + x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3^2 + x1^2*x4^2 + 2*x1*x2^3 + x1*x2^2*x3 + 3*x1*x2^2*x4 + 4*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 2*x1*x2*x4^2 + 4*x1*x3^3 + 3*x1*x3^2*x4 + x1*x3*x4^2 + 4*x1*x4^3 + 4*x2^3*x3 + 2*x2^3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 2*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3^2*x4^2 + 4*x3*x4^3 + 2*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_5 GPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    return
end

function bench_gpu_K3_7()
    GC.gc()
    n = 4
    p = 7
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 5*x1^4 + 6*x1^3*x2 + 2*x1^3*x3 + 3*x1^3*x4 + 4*x1^2*x2^2 + 3*x1^2*x2*x4 +
    2*x1^2*x3^2 + 3*x1^2*x3*x4 + 6*x1^2*x4^2 + 4*x1*x2^2*x3 + 6*x1*x2^2*x4 + 2*x1*x2*x3^2 +
    3*x1*x2*x3*x4 + 5*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + 5*x1*x3*x4^2 + 6*x2^4 +
    5*x2^3*x3 + 3*x2^2*x3^2 + 6*x2^2*x3*x4 + 3*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 +
    3*x2*x4^3 + 5*x3^4 + 6*x3^2*x4^2 + 6*x3*x4^3 + 3*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_7 GPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    return
end

function bench_gpu_K3_11()
    GC.gc()
    n = 4
    p = 11
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 10*x1^4 + x1^3*x2 + 6*x1^3*x3 + 3*x1^3*x4 + x1^2*x2^2 + 9*x1^2*x2*x3 + 6*x1^2*x2*x4 + 6*x1^2*x3^2 + 8*x1^2*x3*x4 + 4*x1^2*x4^2 + 3*x1*x2^3 + 7*x1*x2^2*x3 + 3*x1*x2^2*x4 + 7*x1*x2*x3^2 + 9*x1*x2*x3*x4 + 8*x1*x2*x4^2 + 7*x1*x3^3 + x1*x3*x4^2 + 7*x1*x4^3 + x2^4 + 3*x2^3*x3 + 7*x2^3*x4 + 5*x2^2*x3^2 + 7*x2^2*x3*x4 + 8*x2^2*x4^2 + 8*x2*x3^3 + 5*x2*x3^2*x4 + x2*x3*x4^2 + 9*x2*x4^3 + 7*x3^4 + 4*x3^3*x4 + 4*x3^2*x4^2 + 3*x3*x4^3

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_11 GPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    return
end

function bench_gpu_K3_13()
    GC.gc()
    n = 4
    p = 13
    R, vars = polynomial_ring(GF(p), n)

    (x1, x2, x3, x4) = vars

    h5 = 11*x1^4 + 4*x1^3*x2 + 12*x1^3*x3 + 4*x1^3*x4 + 6*x1^2*x2^2 + 10*x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 7*x1^2*x3*x4 + 4*x1^2*x4^2 + 6*x1*x2^3 + 11*x1*x2^2*x3 + 7*x1*x2^2*x4 + 8*x1*x2*x3^2 + 10*x1*x2*x4^2 + x1*x3^3 + 9*x1*x3^2*x4 + 8*x1*x3*x4^2 + 11*x1*x4^3 + 4*x2^4 + 8*x2^3*x3 + 5*x2^2*x3*x4 + 7*x2^2*x4^2 + 8*x2*x3^3 + 6*x2*x3^2*x4 + 5*x2*x4^3 + 2*x3^4 + 10*x3^3*x4 + 8*x3^2*x4^2 + 10*x3*x4^3 + 6*x4^4

    pregen = MMPSingularities.pregen_qfsheight(n, p)

    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    # prime jitter
    qfs_height_fn(h5)
    MMPSingularities.reset_times()
    for i in 1:10
        qfs_height_fn(h5)
    end

    hrule()
    println("K3_13 GPU Times: ")
    MMPSingularities.display_times()
    MMPSingularities.reset_times()
    return
end

run_benchmarks()