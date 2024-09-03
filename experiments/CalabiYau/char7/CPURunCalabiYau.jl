using Distributed

addprocs(16)

@everywhere begin
    using CSV
    using DataFrames
    using Oscar
    using Dates

    include("../../../src/MMPSingularities.jl")

    function run_experiment(numVars = 4, prime = 7, howHigh = 9, time = Hour(6))
	R, vars = polynomial_ring(GF(prime), numVars)

        dfheights = DataFrame(I = Int[], II = Int[], III = Int[], IV = Int[], V = Int[], VI = Int[], VII = Int[], VIII = Int[], IX = Int[], X = Int[], inf = Int[])
        df = DataFrame(QFSheight = Int[], polynomial = FqMPolyRingElem[])

        samples = 0
        println("Process $(getpid()) starting...")

        heights = zeros(Int, 11)
	startTime = now()
        while (now() - startTime) < time
            poly = MMPSingularities.random_homog_poly_mod(prime, vars, numVars)
	    sampletime = @timed begin            
		height = MMPSingularities.quasiFSplitHeight_CY_lift_sort(prime, poly, 10)
	    end 

	    if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end

            if height >= 8
                push!(df, [height, poly])
            end
            samples += 1

            if samples % 250 == 0
                println("Process $(getpid()): $samples Samples completed! Writing to file...")
                CSV.write("experiments/CalabiYau/char7/cpuheights.csv", df, append=true)
                CSV.write("experiments/CalabiYau/char7/cpuheightsbargraph.csv", dfheights, append=true)
                empty!(df)
                empty!(dfheights)
                heights = zeros(Int, length(heights))
            end
        end

        push!(dfheights, heights)
        CSV.write("experiments/CalabiYau/char7/cpuheights.csv", df, append=true)
        CSV.write("experiments/CalabiYau/char7/cpuheightsbargraph.csv", dfheights, append=true)
    end
end

@sync for pid in workers()
    @async remotecall_wait(run_experiment, pid)
end
