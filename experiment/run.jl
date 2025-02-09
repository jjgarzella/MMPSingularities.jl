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

function write_to_runlog(string)
    open("runlog.txt", "a") do file
	write(file, string)
    end
    return nothing
end

function run_experiment(heights, thread_dicts, experimentThreads, n, p)
    R, vars = polynomial_ring(GF(p), n)
    
    (x1, x2, x3, x4) = vars
    # restricted_mons = [x1^4, x1^3*x2, x1^3*x3, x1^3*x4, x1^2*x2^2, x1^2*x2*x3, x1^2*x2*x4, x1^2*x3^2, x1^2*x3*x4, x1^2*x4^2, x1*x2^3, x1*x2^2*x3, x1*x2^2*x4, x1*x2*x3^2, x1*x2*x3*x4, x1*x2*x4^2, x1*x3^3, x1*x3^2*x4, x1*x3*x4^2, x1*x4^3, x2^3*x3, x2^3*x4, x2^2*x3^2, x2^2*x3*x4, x2^2*x4^2, x2*x3^3, x2*x3^2*x4, x2*x3*x4^2, x2*x4^3, x3^3*x4, x3^2*x4^2, x3*x4^3]
    
    # randompoly() = p == 7 ? MMPSingularities.random_homog_poly_mod_restricted(p, vars, restricted_mons) : MMPSingularities.random_homog_poly_mod(p, vars, n)
    
    # pregen = p == 7 ? MMPSingularities.pregen_qfsheight(n, p, true) : MMPSingularities.pregen_qfsheight(n, p)

    randompoly() = MMPSingularities.random_homog_poly_mod(p, vars, n)
    pregen = MMPSingularities.pregen_qfsheight(n, p)

    Threads.@threads for i in 1:experimentThreads
        println("Thread $(Threads.threadid()) started...")
        samples = 0
        while true
            f = randompoly()
            samples += 1
            height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p, f, 10, pregen)
            
            if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end
            if height >= 3
                push!(thread_dicts, Dict("height" => height, "polynomial" => string(f)))
            end
        end
    end
end

function periodic_writer(heights, thread_dicts, p)
    # Write every x seconds
    x = 150
    while true
        sleep(x)
        heightstablename = "K3C$(p)Heights"
        polystablename = "K3C$(p)Polys"
        
        save_to_database(client, heights, heightstablename, thread_dicts, polystablename)
        totalSamples = sum(heights)
        SPS = totalSamples / x
        write_to_runlog("Processed $SPS samples per second in past $x seconds, threads: $(Threads.nthreads()) \n")
        fill!(heights, 0)
        empty!(thread_dicts)
    end
end

function run(n, p)
    println("Threads.nthreads(): $(Threads.nthreads())")
    experimentThreads = Threads.nthreads() - 1
    heights = zeros(Int, 11)
    
    thread_dicts = Vector{Dict{String, Any}}()
    
    @spawn periodic_writer(heights, thread_dicts, p)
    
    run_experiment(heights, thread_dicts, experimentThreads, n, p)
end

run(4, 13)