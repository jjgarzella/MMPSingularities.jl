using Distributed
addprocs(4)  # Add 4 worker processes

@everywhere begin
    using CSV
    using DataFrames
    using Oscar
    using Dates
    include("../../../src/MMPSingularities.jl")
    using .MMPSingularities

    function run_experiment(numVars = 4, prime = 7, howHigh = 9, time = Second(100))
        R, vars = polynomial_ring(GF(prime), numVars)

        startTime = now()
        dfheights = DataFrame(I = Int[], II = Int[], III = Int[], IV = Int[], V = Int[], VI = Int[], VII = Int[], VIII = Int[], IX = Int[], X = Int[], inf = Int[])
        df = DataFrame(QFSheight = Int[], polynomial = FqMPolyRingElem[])

        samples = 0
        println("Process $(getpid()) starting...")

        heights = zeros(Int, 11)

        while (now() - startTime) < time
            poly = MMPSingularities.random_homog_poly_mod(prime, vars, numVars)
            height = MMPSingularities.quasiFSplitHeight_CY_lift_matrix_combined(prime, poly, 10)
            if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end

            if height >= 8
                push!(df, [height, poly])
            end
            samples += 1

            println("Process $(getpid()): $samples Samples completed")

            if samples % 1000 == 0
                println("Process $(getpid()): $samples Samples completed")
                CSV.write("experiments/CalabiYau/char7/heights.csv", df, append=true)
                CSV.write("experiments/CalabiYau/char7/heightsbargraph.csv", dfheights, append=true)
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

@distributed for _ in 1:nworkers()
    run_experiment()
end

