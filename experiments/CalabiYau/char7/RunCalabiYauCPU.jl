using Distributed
addprocs()

@everywhere begin
    println("\nStarting imports...")
    include("../../../src/MMPSingularities.jl")
    using .MMPSingularities

    using CSV
    using DataFrames
    using Oscar
    using Dates

    println("Imports finished!")

    function run_experiment(process_id; numVars = 4, prime = 7, howHigh = 9, time = Second(100))
        R, vars = polynomial_ring(GF(prime), numVars)

        x1, x2, x3, x4 = vars

        startTime = now()
        dfheights = DataFrame(I = Int[], II = Int[], III = Int[], IV = Int[], V = Int[], VI = Int[], VII = Int[], VIII = Int[], IX = Int[], X = Int[], inf = Int[])
        df = DataFrame(QFSheight = Int[], polynomial = FqMPolyRingElem[])

        samples = 0
        println("Process $process_id starting...")

        heights = zeros(Int, 11)

        while (now() - startTime) < time
            poly = MMPSingularities.random_homog_poly_mod(prime, vars, numVars)
            height = MMPSingularities.quasiFSplitHeight_CY_lift_matrix_combined(prime, poly, 10)
            if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end

            if height >= 5
                push!(df, [height, poly])
            end
            samples += 1

            if samples % 1000 == 0
                println("Process $process_id: $samples Samples completed")
                lock(write_lock) do
                    CSV.write("experiments/CalabiYau/char7/heights.csv", df, append=true)
                    CSV.write("experiments/CalabiYau/char7/heightsbargraph.csv", dfheights, append=true)
                end
                empty!(df)
                empty!(dfheights)
                heights = zeros(Int, length(heights))
            end
        end

        push!(dfheights, heights)
        lock(write_lock) do
            CSV.write("experiments/CalabiYau/char7/cpuheights.csv", df, append=true)
            CSV.write("experiments/CalabiYau/char7/cpuheightsbargraph.csv", dfheights, append=true)
        end

        println("Process $process_id finished! $samples samples computed.")
    end
end

@distributed for process_id in 1:nworkers()
    run_experiment(process_id)
end