include("../../src/MMPSingularities.jl")
using .MMPSingularities

using CSV
using DataFrames
using Oscar
using Dates

function run_experiment(;numVars = 4, prime = 5, howHigh = 5, outputFileName = "heights.csv", time = Second(100))
    R, vars = polynomial_ring(GF(prime), numVars)

    startTime = now()
    df = DataFrame(QFSheight = Int[], polynomial = FqMPolyRingElem[])
    pregen = MMPSingularities.pregen_delta1(numVars, prime)
    samples = 0
    while (now() - startTime) < time
        numTerms = rand(1:binomial(2 * numVars - 1, numVars))
        poly = MMPSingularities.random_homog_poly_mod_k_coefs(prime, vars, numVars, numTerms)
        try
            height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(prime, poly, 10, pregen)
            if height >= howHigh && height <= 10
                push!(df, [height, poly])
            end
        catch (e)
            str = "$e $poly\n"
            open("experiments/CalabiYau/erroredpolynomials.txt", "a") do file
                write(file, str)
            end
        end
        samples += 1

        # So that memory doesn't explode
        if samples % 1000 == 0
            println("$samples samples completed")
            CSV.write("experiments/CalabiYau/$(outputFileName)", df, append=true)
            empty!(df)
        end
    end

    CSV.write("experiments/CalabiYau/$(outputFileName)", df, append=true)
    empty!(df)
    
    open("experiments/CalabiYau/experimentlog.log", "a") do file
        str = "Experiment ran for $time on $(now()), computed $samples samples\n"
        write(file, str)
    end
end