include("../src/MMPSingularities.jl")

using Oscar
using CUDA
using Base.Threads
using Dates

include("supabase.jl")

client = include("key.jl")

function save_to_database(client, heights, heightstablename, all_dicts, polystablename)
    # println("Saving to database...")
    update_heights(client, heights, heightstablename)
    add_polys(client, all_dicts, polystablename)

    return nothing
end

function run_experiment(heights, thread_dicts, experimentThreads, n, p)
    R, vars = polynomial_ring(GF(p), n)
    
    (x1, x2, x3, x4) = vars
    restricted_mons = [x1^3*x2, x1^3*x3, x1^3*x4, x1^2*x2^2, x1^2*x2*x3, x1^2*x2*x4, x1^2*x3^2, x1^2*x3*x4, x1^2*x4^2, x1*x2^3, x1*x2^2*x3, x1*x2^2*x4, x1*x2*x3^2, x1*x2*x3*x4, x1*x2*x4^2, x1*x3^3, x1*x3^2*x4, x1*x3*x4^2, x1*x4^3, x2^3*x3, x2^3*x4, x2^2*x3^2, x2^2*x3*x4, x2^2*x4^2, x2*x3^3, x2*x3^2*x4, x2*x3*x4^2, x2*x4^3, x3^3*x4, x3^2*x4^2, x3*x4^3, x4^4]
    
    randompoly() = p >= 7 ? MMPSingularities.random_homog_poly_mod_restricted(p, vars, restricted_mons) : MMPSingularities.random_homog_poly_mod(p, vars, n)
    
    pregen = MMPSingularities.GPUDelta1.pregen_delta1(n, p)

    Threads.@threads for i in 1:experimentThreads
        println("Thread $(Threads.threadid()) started...")
        samples = 0
        localheights = zeros(Int, 11)
        while true
            f = randompoly()
            samples += 1
            isfsplit, fpminus1 = MMPSingularities.isFSplit2(p, f)
            if isfsplit
                localheights[1] += 1
            else
                height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu_second_step(p, fpminus1, 10, pregen)
                if height == 11 || height == 12
                    localheights[11] += 1
                else
                    localheights[height] += 1
                end
                if height >= 7 && height <= 10
                    push!(thread_dicts, Dict("height" => height, "polynomial" => string(f)))
                end
            end

            if samples % 1000 == 0
                heights[Threads.threadid(), :] .= localheights
                for i in eachindex(localheights)
                    localheights[i] = 0
                end
            end
        end
    end
end

function periodic_writer(heights, thread_dicts, p)
    sleep(30)
    # Write every x seconds
    x = 10
    while true
        sleep(x)
        reducedHeights = reduce(+, heights; dims = 1)
        heightstablename = "K3C$(p)Heights"
        polystablename = "K3C$(p)Polys"
        
        save_to_database(client, reducedHeights, heightstablename, thread_dicts, polystablename)

        fill!(heights, 0)
        empty!(thread_dicts)
    end
end

function run(n, p)
    println("Threads.nthreads(): $(Threads.nthreads())")
    experimentThreads = Threads.nthreads() - 1
    heights = zeros(Int, Threads.nthreads(), 11)
    
    thread_dicts = Vector{Dict{String, Any}}()
    
    # Start the periodic writer task
    @spawn periodic_writer(heights, thread_dicts, p)
    
    # Run the experiment
    run_experiment(heights, thread_dicts, experimentThreads, n, p)
end

run(4, 3)