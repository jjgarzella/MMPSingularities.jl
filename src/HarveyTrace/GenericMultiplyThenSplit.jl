#using Oscar
#using Combinatorics
##using DeRham

## this code is originally due to Alex Pan


#function gen_exp_vec(n, d, order=:lex)
#    result = Vector{Vector{Int64}}(undef, binomial(n+d-1,d))
#    for i in 1:binomial(n+d-1,d)
#        result[i] = zeros(Int64,n)
#    end
#    if order == :lex
#        for i in 1:n
#            dtemp = copy(d)
#            k = 0
#            while k <= (length(result) - 1)
#                if i > 1
#                    dtemp = copy(d)
#                    for j in 1:n
#                        dtemp = dtemp - result[length(result)-k][j]
#                    end
#                end
#                if i == n && dtemp > 0
#                    result[length(result)-k][i] = dtemp
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 0
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 1 && i > 1
#                    for j in i:n
#                        result[length(result)-(k+j-i)][j] = 1
#                    end
#                    k = k + n - i + 1
#                    continue
#                end
#                dtemp2 = copy(d - dtemp)
#                if i == 1 || dtemp2 == 0
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
#                            result[length(result)-(j+k-1)][i] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
#                        dtemp = dtemp - 1
#                    end
#                else
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                            result[length(result)-(j+k-1)][i] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                        dtemp = dtemp - 1
#                    end
#                end
#            end
#        end
#    end
#    return result
#end


#"""
#Computes the monomials in n variables, of degree d, in the
#variable order order, as Oscar polynonmials.
#"""
## This version seems to be unused, and it breaks precompilation
##function compute_monomials(n,d,order=:lex)
##    S, vars = polynomial_ring(ZZ, ["x$i" for i in 1:n])
##
##	compute_monomials(n,d,PR,order=:lex)
##end

#function compute_monomials(n,d,PR,order=:lex)
#    if n < 0 || d < 0
#        return []
#    end
#    gen_mon(gen_exp_vec(n,d,order),base_ring(PR),PR)
#end

## Computes all monomials of degree `d` in `n` variables.
##function compute_monomials(n, d)
##    S, vars = polynomial_ring(ZZ, ["x$i" for i in 1:n])
##
##    result = []
##    
##    function backtrack(start, current)
##        if length(current) == d
##            push!(result, prod(S(var) for var in current))
##            return
##        end
##
##        for i in start:n
##            backtrack(i, [current..., vars[i]])
##        end
##    end
##
##    backtrack(1, [])
##
##    result
##end

#"""
#polynomial_to_vector(f, n, R, PR; order=:lex)

#Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.

#f - an oscar polynomial
#n - the number of variables (minus 1???)
#R - the base ring
#PR - the polynomial ring, i.e. it is parent(f)
#order - a symbol which denotes the term order
#"""
#function polynomial_to_vector(f, n, R, PR, order=:lex)
#    vars = gens(PR)

#    #TODO: if f turns out to be zero, we don't know what the degree should be.
#    #
#    #How best to fix this?
#    d = total_degree(f)

#    mon = compute_monomials(n, d,PR,order)
#    res = fill(R(0), length(mon))
#    for i in eachindex(mon)
#        res[i] = coeff(f, mon[i])
#    end

#    res
#end

function vector(f,d,order=:lex)
    R = parent(f)
    n = length(gens(R))
  
    F = coefficient_ring(R)
    f == zero(R) && return zeros(F,dim_of_homog_polys(n,d))
    @assert d == total_degree(f) "Expect d to be the degree of f"
    polynomial_to_vector(f, n, F, R,order)
end

function allmonomialcombos(vars,deg)
    repeated_vars = repeat(vars,deg)
    multiset_combinations(repeated_vars,deg)
  
end#function

function generic_homog_poly(R)
    nVars = length(gens(R))
  
    var_combos = collect(allmonomialcombos(gens(R),nVars))
    nMons = length(var_combos)
  
    # the line where it all happens
    coefs = gens(R.base_ring)
    @assert length(coefs) == length(var_combos)
    res = zero(R)
    for i in 1:nMons
      res = res + coefs[i] * prod(var_combos[i])
    end
  
    res
end#fucntion

function lift_poly(f::FqMPolyRingElem)
    return map_coefficients(x -> lift(ZZ, x), f)
end

function convert_to_FqMPolyRingElem(R::FqMPolyRing, poly::ZZMPolyRingElem)
    return map_coefficients(x -> R.base_ring(x), poly)
end

# function Δ₁l(p,poly)
#     R = parent(poly)
#     # display(R)
#     # display(fieldnames(typeof(R)))
#     # display(R.base_ring)
#     # display(typeof(R.base_ring))
#     # throw()
#     originallift = map_coefficients(x -> lift_poly(x),poly)
  
#     nocrossterms = sum(terms(originallift) .^p)
#     withcrossterms = originallift^p
  
#     crossterms = withcrossterms - nocrossterms
    
#     Δlift = map_coefficients(x -> divexact(x,p),crossterms)
#     Δ = map_coefficients(x -> convert_to_FqMPolyRingElem(R.base_ring, x), Δlift)

#     return Δ
# end#function

function encode_degs(degs, bits)
    result = zeros(UInt64, size(degs, 2))
    for i in eachindex(result)
        result[i] = base2kron(view(degs, :, i), bits)
    end

    return result
end

function base2kron(vec, bits)
    result = zero(UInt64)
    for i in eachindex(vec)
        result += vec[i] << (bits * (length(vec) - i))
    end
    return result
end

function base2divkron(num::T, m::T, numVars::Int, bits::Int) where T<:Unsigned
    result = zero(T)
    mask = (one(T) << bits) - one(T)
    for i in 0:(numVars - 1)
        element = (num >> (bits * i)) & mask
        divided = element ÷ m
        result += divided << (bits * i)
    end
    return result
end

function base2modkron(num::T, m::T, numVars::Int, bits::Int) where T<:Unsigned
    result = zero(T)
    mask = (one(T) << bits) - one(T)
    for i in 0:(numVars - 1)
        element = (num >> (bits * i)) & mask
        modded = element % m
        result += modded << (bits * i)
    end
    return result
end

function find_next_pminus1(num::T, nvars::Int, bits::Int, p::T) where T<:Number
    result = zero(T)
    mask = (one(T) << bits) - one(T)
    total = zero(T)
    added = zero(T)
    for i in 0:(nvars - 1)
        element = (num >> (bits * i)) & mask
        adjust = p - one(T) - (element % p)
        total += adjust
        added += adjust << (bits * i)
        result += (element + adjust) << (bits * i)
    end
    return result, added, total
end

function wics(n, k)
    x = fill(0, k)
    x[1] = n
    result = zeros(Int, k, binomial(n + k - 1, k - 1))
    idx = 1
    while true
        view(result, :, idx) .= x
        idx += 1
        v = x[end]
        if n == v
            break
        end
        x[end] = 0
        j = k - 1
        while x[j] == 0
            j -= 1
        end
        x[j] -= 1
        x[j + 1] = 1 + v
    end

    return result
end

function generic_matrix_of_multiply_then_split(f,p)
    #p = f.parent.base_ring.data.n
    p = UInt(p)
    coeffs = collect(coefficients(f))
    numVars = length(gens(f.parent))
    @assert numVars == 4
    bits = 16
    encodedDegs = encode_degs(f.exps, 16)
    @assert length(coeffs) == length(encodedDegs)
    
    d = Int(numVars * (p - 1))

    mons = gen_exp_vec(numVars, d)
    mons = reduce(hcat, mons)
    nMons = size(mons, 2)

    result = zeros(f.parent.base_ring, nMons, nMons)
    
    kron(vec) = base2kron(vec, bits)
    div_kron(n, m) = base2divkron(n, m, numVars, bits)
    mod_kron(n, m) = base2modkron(n, m, numVars, bits)

    reverseMons = Dict{UInt,Int}()
    encodedMons = encode_degs(mons, bits)
    for i in eachindex(encodedMons)
        reverseMons[encodedMons[i]] = i
    end

    weakintegercompositions = [encode_degs(wics(i, numVars) .* p, bits) for i in 0:fld(d, p)]

    relevant = kron(fill(p - 1, numVars))

    for term in eachindex(encodedDegs)
        initialDeg, initialMon, howmuchadded = find_next_pminus1(encodedDegs[term], numVars, bits, p)
        println("$howmuchadded")
        println("$(d - howmuchadded) / $p")
        weaks = divexact(d - howmuchadded, p)
        thingstoadd = weakintegercompositions[weaks + 1]
        for i in eachindex(thingstoadd)
            mon = initialMon + thingstoadd[i]
            new_exv = div_kron(initialDeg + thingstoadd[i] - relevant, p)
            result[reverseMons[new_exv], reverseMons[mon]] = coeffs[term]
        end
    end

    return result
end

function index_of_term_not_in_frobenius_power_CY(p,n,order=:lex)
    R, vars = polynomial_ring(GF(p),n)
    
    crit_term = prod(vars .^ (p-1))
  
    # perhaps assert this has only one element?
    findfirst(vector(crit_term,total_degree(crit_term)) .!= 0)
end

function run()
    n = 4
    k = n
    p = 2

    homogPolyDim = binomial(n + k - 1, k - 1)
    R, coeffVars = polynomial_ring(GF(p), homogPolyDim)
    R2, (x, y, z, w) = polynomial_ring(R, [:x, :y, :z, :w])

    f = generic_homog_poly(R2)

    g = f^(p - 1)
    @assert n == 4
    open("coefficients.txt", "a") do file
        write(file, "coefficients of critial term for h = 1: $(coeff(g, [p-1, p-1, p-1, p-1])) \n\n")
    end

    Δ = Δ₁l(p, g)
    # display(collect(coefficients(Δ)))
    mat = generic_matrix_of_multiply_then_split(Δ)
    start_vector = vector(g, (n * (p - 1)))
    critical_index = index_of_term_not_in_frobenius_power_CY(p, n)
    # display(start_vector)
    # println(count(x -> x == zero(R), mat))
    # println(size(mat, 1)^2)

    KTYideal_n_new_gen = mat * start_vector
    h = 2

    while h ≤ 10 
        open("coefficients.txt", "a") do file
            write(file, "coefficients of critial term for h = $h: $(KTYideal_n_new_gen[critical_index]) \n\n")
        end
        KTYideal_n_new_gen = mat * KTYideal_n_new_gen
        h = h + 1
    end
end

