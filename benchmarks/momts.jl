include("../src/MMPSingularities.jl")
include("../src/RandomPolynomials.jl")

using CUDA
using Oscar

"""
This is copied and pasted from 
DeRham.jl, so it may not be up to date.
"""
function gen_exp_vec(n, d, order=:lex)
    result = Vector{Vector{Int64}}(undef, binomial(n+d-1,d))
    for i in 1:binomial(n+d-1,d)
        result[i] = zeros(Int64,n)
    end
    if order == :lex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :neglex
        for i in 1:n
            dtemp = copy(d)
            k = 1
            while k <= length(result)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[k+j-i][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :invlex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][n-i+1] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][n-j+1] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    else
        throw(ArgumentError("Unsupported order '$order'"))
    end
    return result
end

function test_matrix()
    n = 4
    primes = [3, 5, 7]#, 11, 13]
    
    for p in primes
        R, vars = polynomial_ring(GF(p), n)
        (x1, x2, x3, x4) = vars

        #pregen = MMPSingularities.pregen_delta1(n, p)
        println("starting exp vecs")

        #f = random_homog_poly_mod(p, vars, n)
        #fpminus1 = MMPSingularities.HomogeneousPolynomial(f ^ (p - 1))
        #Δ₁fpminus1 = MMPSingularities.delta1(fpminus1, p; pregen = pregen).poly
        exp_vecs = gen_exp_vec(n,n*(p-1)*p)
        println("starting monomials")

        cxt = MPolyBuildCtx(R)
        for vec in exp_vecs
            push_term!(cxt,one(base_ring(R)),vec)
        end
        fake_delta1 = finish(cxt)
        #monomials = [prod(vars .^ vec) for vec in exp_vecs]#[1:div(end,2)]]
        #TODO the following is not efficient enough
        #just import DeRham.jl and use gen_exp_vec
        #monomials = allmonomialcombos(vars,n*(p-1)*p) 
        
        #fake_delta1 = sum(monomials)
        Δ₁fpminus1 = fake_delta1

        println("starting momts pregen")

        momtspregen = MMPSingularities.pregen_MOMTS(n, p)
        
        println("done with momts pregen")

        if p <= 7
            MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
        end
        println("trivial compiled")
        MMPSingularities.matrix_of_multiply_then_split_gpu(Δ₁fpminus1, momtspregen)
        println("gpu trivial compiled")
        MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
        println("merge compiled")
        MMPSingularities.matrix_of_multiply_then_split_wics(Δ₁fpminus1)
        println("wics compiled")
        MMPSingularities.matrix_of_multiply_then_split_wics_gpu(Δ₁fpminus1, momtspregen)
        println("gpu wics compiled")
        times = zeros(Float64, 5)
        for i in 1:10
            # @time mat0 = MMPSingularities.matrix_of_multiply_then_split_correct(Δ₁fpminus1)
            if p <= 7
                trivial = @timed MMPSingularities.matrix_of_multiply_then_split(Δ₁fpminus1)
            else
                trivial = @timed zero(Float64)
            end
            println("trivial done")
            gputrivial = CUDA.@timed MMPSingularities.matrix_of_multiply_then_split_gpu(Δ₁fpminus1, momtspregen)
            println("gpu trivial done")
            merge = @timed MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker(Δ₁fpminus1)
            println("merge done")
            wics = @timed MMPSingularities.matrix_of_multiply_then_split_wics(Δ₁fpminus1)
            println("wics done")
            gpuwics = CUDA.@timed MMPSingularities.matrix_of_multiply_then_split_wics_gpu(Δ₁fpminus1, momtspregen)
            println("gpuwics done")

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
