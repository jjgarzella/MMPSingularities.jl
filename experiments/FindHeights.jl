include("../src/MMPSingularities.jl")

using Oscar

# Gets 
function run_experiment(p, samples)
    n = 4
    R, vars = polynomial_ring(GF(p), n)
    (x1, x2, x3, x4) = vars

    pregen = MMPSingularities.pregen_qfsheight(n, p)
    qfs_height_fn(x) = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, x, 10, pregen)

    randompoly() = MMPSingularities.random_homog_poly_mod(p, vars, n)

    results = zeros(Int, 11)
    i = 0
    while i < samples
        f = randompoly()
        height = MMPSingularities.quasiFSplitHeight_CY_lift_wics_gpu(p, f, 10, pregen)
        
        if height == 11 || height == 12
            results[11] += 1
            # push!(results[11], string(f))
        else
            results[height] += 1
            # push!(results[height], string(f))
        end
        i += 1
    end

    println("results: ")
    println(results)

    return results
end

