println("\nStarting imports...")
include("../../src/MMPSingularities.jl")
using .MMPSingularities

using CSV
using DataFrames
using Oscar
using Dates

println("Imports finished!")
function run_experiment(;numVars = 4, prime = 5, howHigh = 9, time = Second(100))
    R, vars = polynomial_ring(GF(prime), numVars)

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
        poly = MMPSingularities.random_homog_poly(prime, vars, numVars)
        try
            height = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(prime, poly, 10, pregen)
            if height == 11 || height == 12
                heights[11] += 1
            else
                heights[height] += 1
            end


            # This can be commented out if quasiFSplitHeight_CY_lift_sort_gpu() gets fixed
            if height >= howHigh && height <= 10
	    	# Because of the bug, we also need to check with a guaranteed working version
		        realheight = MMPSingularities.quasiFSplitHeight_CY_gpu(prime, poly, 10, pregen)
                if (realheight >= height)
                    println("height $realheight found!")
                    push!(df, [realheight, poly])
                    if realheight > height
		        str = "Matrix method returned $height, safe method returned $height, polynomial: $poly"
		        println("Underestimated! $str")
                        open("experiments/CalabiYau/underestimated.txt", "a") do file
                            write(file, "Matrix method returned $height, safe method returned $height, polynomial: $poly\n")
                        end
                    end
                else
		    str = "Matrix method returned $height, safe method returned $height, polynomial: $poly"
                    println("Overestimated! $str")
		    open("experiments/CalabiYau/overestimated.txt", "a") do file
                        write(file, "Matrix method returned $height, safe method returned $height, polynomial: $poly\n")
                    end
                end
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
            println("$samples Samples completed")
        end
        if samples % 1000 == 0
            push!(dfheights, heights) 
            CSV.write("experiments/CalabiYau/heights.csv", df, append=true)
            CSV.write("experiments/CalabiYau/heightsbargraph.csv", dfheights, append=true)
            empty!(df)
            empty!(dfheights)
            heights = zeros(Int, length(heights))
        end
    end

    CSV.write("experiments/CalabiYau/heights.csv", df, append=true)
    CSV.write("experiments/CalabiYau/heightsbargraph.csv", dfheights, append=true)

    println("Experiment finished! $samples samples computed.")
end
