println("\nStarting imports...")
include("../../../src/MMPSingularities.jl")
using .MMPSingularities

using CSV
using DataFrames
using Oscar
using Dates

println("Imports finished!")
function run_experiment(;numVars = 4, prime = 7, howHigh = 9, time = Second(100))
    R, vars = polynomial_ring(GF(prime), numVars)

    x1, x2, x3, x4 = vars

    mons = [x1^3*x2, x1^3*x3, x1^3*x4, x1^2*x2^2, x1^2*x2*x3, x1^2*x2*x4, x1^2*x3^2, x1^2*x3*x4, x1^2*x4^2, x1*x2^3, x1*x2^2*x3, x1*x2^2*x4, x1*x2*x3^2, x1*x2*x3*x4, x1*x2*x4^2, x1*x3^3, x1*x3^2*x4, x1*x3*x4^2, x1*x4^3, x2^3*x3, x2^3*x4, x2^2*x3^2, x2^2*x3*x4, x2^2*x4^2, x2*x3^3, x2*x3^2*x4, x2*x3*x4^2, x2*x4^3, x3^3*x4, x3^2*x4^2, x3*x4^3, x4^4]

    startTime = now()
    # julia gets mad at me when I try to have plain numbers as column titles
    dfheights = DataFrame(I = Int[], II = Int[], III = Int[], IV = Int[], V = Int[], VI = Int[], VII = Int[], VIII = Int[], IX = Int[], X = Int[], inf = Int[])
    df = DataFrame(QFSheight = Int[], polynomial = FqMPolyRingElem[])
    println("Pregenerating...")
    pregen = MMPSingularities.pregen_delta1(numVars, prime)
    samples = 0
    println("Starting...")
    heights = zeros(Int, 11)

    while (now() - startTime) < time
        poly = MMPSingularities.random_homog_poly_mod_restricted(prime, vars, mons)
        height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(prime, poly, 10, pregen)
        if height == 11 || height == 12
            heights[11] += 1
        else
            heights[height] += 1
        end

        if height >= 8
            push!(df, [height, poly])
        end
        # try
        #     height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(prime, poly, 10, pregen)
        #     if height == 11 || height == 12
        #         heights[11] += 1
        #     else
        #         heights[height] += 1
        #     end
        # catch (e)
        #     str = "$e $poly\n"
        #     open("experiments/CalabiYau/char7/erroredpolynomials.txt", "a") do file
        #         write(file, str)
        #     end

        # end
        samples += 1
        # So that memory doesn't explode
        # if samples % 1000 == 0
        #     println("$samples Samples completed")
        # end
        if samples % 1000 == 0
            println("$samples Samples completed")
            push!(dfheights, heights) 
            CSV.write("experiments/CalabiYau/char7/heights.csv", df, append=true)
            CSV.write("experiments/CalabiYau/char7/heightsbargraph.csv", dfheights, append=true)
            empty!(df)
            empty!(dfheights)
            heights = zeros(Int, length(heights))
        end
    end
    push!(dfheights, heights) 
    CSV.write("experiments/CalabiYau/char7/heights.csv", df, append=true)
    CSV.write("experiments/CalabiYau/char7/heightsbargraph.csv", dfheights, append=true)

    println("Experiment finished! $samples samples computed.")
end
