include("src/MMPSingularities.jl")

using Oscar
using CUDA
using Base.Threads
using Dates

function run_experiment(experimentThreads, n, p, localSamples)
    R, vars = polynomial_ring(GF(p), n)
    
    (x1, x2, x3, x4) = vars
    # restricted_mons = [x1^4, x1^3*x2, x1^3*x3, x1^3*x4, x1^2*x2^2, x1^2*x2*x3, x1^2*x2*x4, x1^2*x3^2, x1^2*x3*x4, x1^2*x4^2, x1*x2^3, x1*x2^2*x3, x1*x2^2*x4, x1*x2*x3^2, x1*x2*x3*x4, x1*x2*x4^2, x1*x3^3, x1*x3^2*x4, x1*x3*x4^2, x1*x4^3, x2^3*x3, x2^3*x4, x2^2*x3^2, x2^2*x3*x4, x2^2*x4^2, x2*x3^3, x2*x3^2*x4, x2*x3*x4^2, x2*x4^3, x3^3*x4, x3^2*x4^2, x3*x4^3]
    
    # randompoly() = p == 7 ? MMPSingularities.random_homog_poly_mod_restricted(p, vars, restricted_mons) : MMPSingularities.random_homog_poly_mod(p, vars, n)
    
    # pregen = p == 7 ? MMPSingularities.pregen_qfsheight(n, p, true) : MMPSingularities.pregen_qfsheight(n, p)

    randompoly() = MMPSingularities.random_homog_poly_mod(p, vars, n)
    pregen = MMPSingularities.pregen_qfsheight(n, p)

    Threads.@threads for i in 1:experimentThreads
        println("Thread $(Threads.threadid()) started...")
        while true
            f = randompoly()
            height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, f, 10, pregen)

            localSamples[Threads.threadid()] += 1
        end
    end
end

function periodic_writer(localSamples)
    # Write every x seconds
    x = 5
    lastSum = 0
    while true
        sleep(x)
        currentSum = sum(localSamples)
        diff = currentSum - lastSum
        SPS = diff / x
        println("Processed $SPS samples per second in past $x seconds, threads: $(Threads.nthreads())")
        lastSum = currentSum
    end
end

function run(n, p)
    println("Threads.nthreads(): $(Threads.nthreads())")
    experimentThreads = Threads.nthreads() - 1
    localSamples = zeros(Int, Threads.nthreads())
    
    # Start the periodic writer task
    @spawn periodic_writer(localSamples)
    
    # Run the experiment
    run_experiment(experimentThreads, n, p, localSamples)
end

run(4, 5)