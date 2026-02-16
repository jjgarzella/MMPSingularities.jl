using MMPSingularities

using Oscar
using CUDA
using Base.Threads
using Dates

include("supabase.jl")

client = include("key.jl")

function check_smoothness(f)
    p = characteristic(parent(f))
    nVars = length(gens(parent(f)))
    graded, _ = grade(parent(f))
    R, _ = quo(graded, ideal(graded, [graded(f)]))
    V = proj(R)
    is_smooth(V)
end

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
    cutoff = 0 # cutoff required to write to database
    if p == 3
	cutoff = 9
    elseif p == 5
	cutoff = 8
    elseif p == 7
	cutoff = 8
    elseif p == 11
	cutoff = 4
    elseif p == 13
	cutoff = 3
    end
    R, vars = polynomial_ring(GF(p), n)
    
    (x1, x2, x3, x4) = vars

    randompoly() = MMPSingularities.random_homog_poly_mod(p, vars, n)
    pregen = MMPSingularities.pregen_qfsheight(n, p)

    Threads.@threads for i in 1:experimentThreads
        println("Thread $(Threads.threadid()) started...")
        samples = 0
        while true
            f = randompoly()
            samples += 1
            height = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, f, 10, pregen)
            
            if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end
            if height >= cutoff
		if (check_smoothness(f))
		    push!(thread_dicts, Dict("height" => height, "polynomial" => string(f)))
		end
            end
        end
    end
end

function periodic_writer(heights, thread_dicts, p)
    # Write every x seconds
    x = 150
    while true
        sleep(x)
        heightstablename = "smooth_K3C$(p)Heights"
        polystablename = "smooth_K3C$(p)Polys"
        
        save_to_database(client, heights, heightstablename, thread_dicts, polystablename)
        totalSamples = sum(heights)
        SPS = totalSamples / x
        write_to_runlog("p: $(p). Processed $SPS samples per second in past $x seconds, threads: $(Threads.nthreads()) \n")
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
    wait(Condition())
end

println("hello")
