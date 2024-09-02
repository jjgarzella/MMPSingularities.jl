include("MMPSingularities.jl")

using Distributed

addprocs(4)

using Oscar
using CUDA
using Test
using Profile
using Dates

function run_test()
    qfs_height_fn = MMPSingularities.quasiFSplitHeight_CY_lift_sort

    R, vars = polynomial_ring(GF(7), 4)
    f = MMPSingularities.random_homog_poly_mod(7, vars, 4)
    println("f: $f")

    height = qfs_height_fn(7, f, 10)

    println("height: $height")
end

run_test()