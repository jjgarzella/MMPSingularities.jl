# This file provides an easy interface to time gpu delta1 and 
# matmul on the cluster

include("src/MMPSingularities.jl")
include("src/RandomPolynomials.jl")
using CUDA
using Oscar

function time_delta1(n, p)
    R, vars = polynomial_ring(GF(p), n)

    for i in 1:10
        f = random_homog_poly_mod(p, vars, n)

        fpminus1 = f ^ (p - 1)

        pregen = MMPSingularities.pregen_delta1(n, p)

        fpminus1_homog = MMPSingularities.HomogeneousPolynomial(fpminus1)
        CUDA.@time Δ₁fpminus1 = MMPSingularities.memorysafe_delta1(fpminus1_homog, p; pregen = pregen).poly
    end

    return nothing
end

function p11delta1()
    n = 4
    p = 11
    R, vars = polynomial_ring(GF(p), n)

    f = random_homog_poly_mod(p, vars, n)
    println("f: $f")
    fpminus1 = f ^ (p - 1)
    
    pregen = MMPSingularities.pregen_delta1(n, p)

    fpminus1_homog = MMPSingularities.HomogeneousPolynomial(fpminus1)
    CUDA.@time Δ₁fpminus1 = MMPSingularities.memorysafe_delta1(fpminus1_homog, p; pregen = pregen).poly
end

function get_matrix(n, p)

end
