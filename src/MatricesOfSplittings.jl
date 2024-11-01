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

function gen_mon(exp_vec, R, PR)
    result = []
    for i in axes(exp_vec,1)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), exp_vec[i])
        monomial = finish(B)
        push!(result,monomial)
    end
    result
end

function compute_monomials(n,d,PR,order=:lex)
    if n < 0 || d < 0
        return []
    end
    gen_mon(gen_exp_vec(n,d,order),base_ring(PR),PR)
end


"""
Wrapper for polynomial_to_vector

Converts the homogeneous polynomial poly
to a vector.

"""
function vector(f,d,order=:lex)
  R = parent(f)
  n = length(gens(R))

  F = coefficient_ring(R)
  f == zero(R) && return zeros(F,dim_of_homog_polys(n,d))
  @assert d == total_degree(f) "Expect d to be the degree of f"
  polynomial_to_vector(f, n, F, R,order)
end

function polynomial_to_vector(f, n, R, PR, order=:lex)
    vars = gens(PR)

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    d = total_degree(f)

    mon = compute_monomials(n, d,PR,order)
    res = fill(R(0), length(mon))
    for i in eachindex(mon)
        res[i] = coeff(f, mon[i])
    end

    res
end
function convert_to_gpu_representation(p)
    coeffs = coefficients(p)

    # julia by default doesn't realize that "ZZ" is not
    # an array, so insert it as a one-element tuple "(ZZ,)"
    # so that julia will know not to broadcast along it.
    coeffs_as_int64arr = UInt.(lift.((ZZ,),coeffs))

    exp_vecs = leading_exponent_vector.(terms(p))

    # shamelessly taken from 
    # https://discourse.julialang.org/t/how-to-convert-vector-of-vectors-to-matrix/72609/2 
    exponent_mat = reduce(hcat,exp_vecs)

    (coeffs_as_int64arr,exponent_mat)
end

"""
Takes a linear operator L on the space
of homogenous polynomials 
of degree d and computes 
the matrix representing it.

Currently uses lexographical order

L is a function, which is assumed to be a linear
endomoprhism on the vector space of homogeneous
polynomials.

d is the degree of the homogeneous polynomials.

R is the base ring.

"""
function matrix_of_lin_op(L,d,R,order=:lex)

  n = length(gens(R))
  monomials = compute_monomials(n,d,R,order)

  m = length(monomials) # will be an mxm matrix

  i = 0

  matrix = zeros(coefficient_ring(R),m,0)
  for monomial in monomials
    evaled = L(monomial)
    v = vector(evaled,d)
    matrix = [matrix v]

    if leading_exponent_vector(monomial) == [14,1,1,0]
      println(v)
    end

    i = i + 1
    if i % 50 == 0 
      println("50 rows completed")
    end
  end

  matrix
end

"""
Multiplies the two polynomials f and g together
and then applies `polynomial_frobenius_splitting`
to the result.

This algorithm only stores the relevant terms,
forgetting all intermediate ones. 
It *should* use less memory than the usual one.

"""
function multiply_then_split(p,f,g,indices)

  result = zero(f)

  vars = gens(parent(f))


  for i in 1:length(f)
    t = term(f,i)
    for j in 1:length(g)
      u = term(g,j)

      prodterm = t*u

      exps = exponent_vector(prodterm,1)

      if all((exps .% p) .== indices)

        coef = coeff(prodterm,1)

        new_exp_vec = divexact.(exps .- indices,p) # the difision should be exact by the if statement

        newterm = coef * prod(vars .^ new_exp_vec) 

        result = result + newterm
      end

    end

  end

  result
end#function

function multiply_then_split(p,f,g)
    nVars = length(gens(parent(f)))

    multiply_then_split(p,f,g,fill(p-1, nVars))
end#function

function matrix_of_multiply_then_split_correct(poly::FqMPolyRingElem)
    p = Int(poly.parent.data.n)
    n = poly.parent.data.nvars
    θFstar(a) = polynomial_frobenius_generator(p,poly*a)
    m = n * (p - 1)

    M = matrix_of_lin_op(θFstar,m,parent(poly))
    return lift_to_Int64(M)
end

function matrix_of_multiply_then_split(poly::FqMPolyRingElem)
    p = poly.parent.data.n
    n = poly.parent.data.nvars

    d = n * (p - 1)

    coeffs, degs = convert_to_gpu_representation(poly)
    return matrix_of_multiply_then_split(p, coeffs, degs, Int(d))
end

"""
Computes the matrix of the 
linear operator of multiplying
by the polynomnial f with coefficients
coefs and degrees degs and then applying 
polynomial_frobenius_generator
on the vector space of homogeneous polynomials
of degree d

This actually does a double for loop, thus 
it'll have slower time complexity than the
merge-based algorithms below which take
advantage of the order.

coefs - vector of coefficients
degs - 2d array of exponent vectors
"""
function matrix_of_multiply_then_split(p,coefs,degs,d)
    n = size(degs,1)
    mons = gen_exp_vec(n,d)
    reverseDict = Dict(mons[i] => i for i in eachindex(mons))
    mons = reduce(hcat, mons)
    nMons = size(mons, 2)
    result = zeros(eltype(coefs), nMons,nMons)
  
    for i in 1:nMons
      # compute column i
      for tInd in 1:size(degs,2)
  
        relevant = true
        for k in 1:n
          #prod_exp_vec_mod_p[k] = degs[tInd,k] + mons[i][k] % p
          if (degs[k, tInd] + mons[k, i]) % p != p-1
            relevant = false
          end
        end
  
        if relevant
        #if all((degs[tInd,:] .+ mons[i]) .% p .== fill(n,p-1))
          # this is a relevant term
          exv_in_prod = degs[:, tInd] .+ mons[:, i]
          new_exv = divexact.(exv_in_prod .- fill(p - 1, n),p)
          row = reverseDict[new_exv]
          result[row,i] += coefs[tInd]
  
        end
      end
    end
  
    result
end

"""
Lifts a matrix with entries in GF(p) to ZZ and converts the entries
to Julia integers
"""
lift_to_Int64(matrix) = Int64.(map(x -> lift(ZZ,x), matrix))

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

function matrix_of_multiply_then_split_sortmodp_kronecker(poly::FqMPolyRingElem)
    
    p = poly.parent.data.n
    n = poly.parent.data.nvars

    d = Int(n * (p - 1))

    coeffs = get_coeffs(poly)
    degs = get_exps(poly)

    return matrix_of_multiply_then_split_sortmodp_kronecker(p, coeffs, degs, d, n, poly.data.bits)
end

function matrix_of_multiply_then_split_sortmodp_kronecker(p::UInt, coefs::Vector{<:Unsigned}, encodedDegs::Vector{<:Unsigned}, d::Int, numVars::Int, bits::Int)
    mons = gen_exp_vec(numVars,d)
    mons = reduce(hcat,mons)

    nMons = size(mons,2)
    nTerms = length(encodedDegs)

    kron(vec) = base2kron(vec, bits)
    div_kron(n, m) = base2divkron(n, m, numVars, bits)
    mod_kron(n, m) = base2modkron(n, m, numVars, bits)

    reverseMons = Dict{UInt,Int}()
    encodedMons = encode_degs(mons, bits)
    for i in eachindex(encodedMons)
        reverseMons[encodedMons[i]] = i
    end
  
    encodedMonsModP = map(x -> mod_kron(x, p), encodedMons)
    mons_perm = sortperm(encodedMonsModP)
  
    left = true
  
    l = 1 # left index
    r = nMons # right index
  
    
    result = zeros(eltype(coefs),nMons,nMons)
    
    relevant = kron(fill(p - 1, numVars))

    reverseDegs = Dict{UInt,Int}()
    for i in eachindex(encodedDegs)
        reverseDegs[encodedDegs[i]] = i
    end
    encodedDegsModP = map(x -> mod_kron(x, p), encodedDegs)

    
    degs_perm = sortperm(encodedDegsModP)
    while l ≤ nTerms && 1 ≤ r
        monModP = encodedMonsModP[mons_perm[r]]
        termModP = encodedDegsModP[degs_perm[l]]
        if monModP + termModP == relevant
            nMatches = 1
            cmpTerm = encodedDegsModP[degs_perm[l + nMatches]]

            while l + nMatches ≤ nTerms && cmpTerm == termModP
                nMatches += 1
                if l + nMatches ≤ nTerms
                    cmpTerm = encodedDegsModP[degs_perm[l + nMatches]]
                end
            end
    

            cmpMon = encodedMonsModP[mons_perm[r]]
            # loop through all monomials and process each one
    
            while 1 <= r && cmpMon == monModP
            #mon = @view mons[mons_perm[r],:]
            mon = encodedMons[mons_perm[r]]
            for ll = l:(l + nMatches - 1)
                term = encodedDegs[degs_perm[ll]]
                newTerm = div_kron(mon + term - relevant, p)
                newcoefind = reverseDegs[term] 
                newcoef = coefs[newcoefind]
                row = reverseMons[newTerm]
                col = reverseMons[mon]
                result[row,col] += newcoef
            end
    
            r -= 1
            left = true
            
            #if 1 ≤ r
            #    cmpMon = encodedMonsModP[mons_perm[r]]
            #end
            # somehow this is erroring for me - Alex
            # NOTE: it should work if you put the assignment in parens - JJ
            (1 ≤ r) && (cmpMon = encodedMonsModP[mons_perm[r]])
        end
    
        l += nMatches - 1
        else
            if left
                l += 1
                left = false
            else
                r -= 1
                left = true
            end
        end
    end

    result
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

function matrix_of_multiply_then_split_alex(poly::FqMPolyRingElem)
    p = poly.parent.data.n
    n = poly.parent.data.nvars

    d = Int(n * (p - 1))

    coeffs = get_coeffs(poly)
    degs = get_exps(poly)

    return matrix_of_multiply_then_split_alex(p, coeffs, degs, d, n, poly.data.bits)
end

function matrix_of_multiply_then_split_alex(p::UInt, coeffs::Vector{<:Unsigned}, encodedDegs::Vector{<:Unsigned}, d::Int, numVars::Int, bits::Int)
    mons = gen_exp_vec(numVars,d)
    mons = reduce(hcat,mons)

    nMons = size(mons,2)

    result = zeros(eltype(coeffs), nMons, nMons)

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

function matrix_of_multiply_then_split_alex_gpu(poly::FqMPolyRingElem, pregen = nothing)
    
    p = poly.parent.data.n
    n = poly.parent.data.nvars

    if pregen === nothing
        pregen = generate_MOMTS(n, p)
    end

    d = Int(n * (p - 1))

    coeffs = CuArray(get_coeffs(poly))
    degs = CuArray(get_exps(poly))

    return matrix_of_multiply_then_split_alex_gpu(p, coeffs, degs, d, n, poly.data.bits, pregen)
end

function matrix_kernel(p::T, coeffs::CuDeviceVector{<:Integer}, encodedDegs::CuDeviceVector{T}, numVars::Int, weakintegercompositions::CuDeviceVector{T}, lengths::CuDeviceVector{Int}, startindices::CuDeviceVector{Int}, reverseMons::MyMap, bits::Int, d, div_kron, relevant::T, result) where T<:Unsigned
    term = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if term <= length(encodedDegs)
        initialDeg, initialMon, howmuchadded = find_next_pminus1(encodedDegs[term], numVars, bits, p)

        weaks = div(d - howmuchadded, p)
        numthingstoadd = lengths[weaks + 1]
        wicsstartidx = startindices[weaks + 1]
        for i in 0:numthingstoadd - 1
            mon = initialMon + weakintegercompositions[wicsstartidx + i]
            new_exv = div_kron(initialDeg + weakintegercompositions[wicsstartidx + i] - relevant, p)
            result[reverseMons[new_exv], reverseMons[mon]] = coeffs[term]
        end
    end

    return nothing
end

struct MOMTSPregen
    nMons::Int
    reverseMons::MyMap
    weakintegercompositions::CuVector{UInt}
    startindices::CuVector{Int}
    lengths::CuVector{Int}
end

function pregen_MOMTS(n, p)
    n = Int(n)
    p = Int(p)
    d = n * (p - 1)
    mons = gen_exp_vec(n, d)
    mons = reduce(hcat, mons)

    if n == 4
        bits = 16
    elseif n == 5
        bits = 12
    else
        throw("Pregeneration not implemented")
    end

    encodedMons = encode_degs(mons, bits)
    reverseMons = make_dict(encodedMons)
    encodedMons = CuArray(encodedMons)

    weakintegercompositions = [encode_degs(wics(i, n) .* p, bits) for i in 0:fld(d, p)]
    startindices = zeros(Int, length(weakintegercompositions))
    curridx = 1

    for i in eachindex(startindices)
        startindices[i] = curridx
        curridx += length(weakintegercompositions[i])
    end
    startindices = CuArray(startindices)
    lengths = CuArray([length(weakintegercompositions[i]) for i in eachindex(weakintegercompositions)])
    weakintegercompositions = reduce(vcat, weakintegercompositions)
    weakintegercompositions = CuArray(weakintegercompositions)

    return MOMTSPregen(length(encodedMons), reverseMons, weakintegercompositions, startindices, lengths)
end

function matrix_of_multiply_then_split_alex_gpu(p::UInt, coeffs::CuVector{<:Unsigned}, encodedDegs::CuVector{<:Unsigned}, d::Int, numVars::Int, bits::Int, pregen::MOMTSPregen)
    result = CUDA.zeros(eltype(coeffs), pregen.nMons, pregen.nMons)

    kron(vec) = base2kron(vec, bits)
    div_kron(n, m) = base2divkron(n, m, numVars, bits)
    mod_kron(n, m) = base2modkron(n, m, numVars, bits)

    relevant = kron(fill(p - 1, numVars))

    kernel = @cuda launch = false matrix_kernel(p, coeffs, encodedDegs, numVars, pregen.weakintegercompositions, pregen.lengths, pregen.startindices, pregen.reverseMons, bits, d, div_kron, relevant, result)
    config = launch_configuration(kernel.fun)
    threads = min(length(encodedDegs), config.threads)
    blocks = cld(length(encodedDegs), threads)

    kernel(p, coeffs, encodedDegs, numVars, pregen.weakintegercompositions, pregen.lengths, pregen.startindices, pregen.reverseMons, bits, d, div_kron, relevant, result; threads = threads, blocks = blocks)

    return result
end