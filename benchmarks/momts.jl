include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")

using CUDA
using Oscar

function test_matrix()
    n = 4
    primes = [3, 5, 7, 11, 13]
    
    for p in primes
        R, vars = polynomial_ring(GF(p), n)
        (x1, x2, x3, x4) = vars

        pregen = MMPSingularities.pregen_delta1(n, p)

        f = random_homog_poly_mod(p, vars, n)
        fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
        Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly

        momtspregen = MMPSingularities.pregen_MOMTS(n, p)

        if p <= 7
            MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
        end
        MMPSingularities.matrix_of_multiply_then_split_gpu(Δ₁fpminus1, momtspregen)
        MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
        MMPSingularities.matrix_of_multiply_then_split_wics(Δ₁fpminus1)
        MMPSingularities.matrix_of_multiply_then_split_wics_gpu(Δ₁fpminus1, momtspregen)
        times = zeros(Float64, 5)
        for i in 1:10
            # @time mat0 = MMPSingularities.matrix_of_multiply_then_split_correct(Δ₁fpminus1)
            if p <= 7
                trivial = @timed MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
            else
                trivial = zero(Float64)
            end
            gputrivial = CUDA.@timed MMPSingularities.matrix_of_multiply_then_split_gpu(Δ₁fpminus1, momtspregen)
            merge = @timed MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
            wics = @timed MMPSingularities.matrix_of_multiply_then_split_wics(Δ₁fpminus1)
            gpuwics = CUDA.@timed MMPSingularities.matrix_of_multiply_then_split_wics_gpu(Δ₁fpminus1, momtspregen)

            times[1] += trivial.time
            times[2] += gputrivial.time
            times[3] += merge.time
            times[4] += wics.time
            times[5] += gpuwics.time
        end
        times ./= 10
        mystr = """
        p = $p:
        trivial: $(times[1]) s
        gputrivial: $(times[2]) s
        merge: $(times[3]) s
        wics: $(times[4]) s
        gpuwics: $(times[5]) s

        """
        open("matrixtimes.txt", "a") do file
            write(file, mystr) 
        end
    end
end

test_matrix()