lift_to_Int64(matrix) = Int64.(map(x -> lift(ZZ,x), matrix))

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

"""
    kronecker_opt(vec, n, kroneckerPregen::Vector{Int})

Version of kronecker() that utilizes pregeneration, see example in body of
matrix_of_multiply_then_split_sortmodp_kronecker2() to see how it works
"""
function kronecker_opt(vec, n, kroneckerPregen::Vector{Int})#,returntype=Int64)
    s = 0#zero(returntype)
    #println("$n")
    for i in eachindex(kroneckerPregen)
        s = s + vec[i] * kroneckerPregen[i]
    end
    s
end

"""
    mod_kronecker(num, m, numVars, kroneckerPregen::Vector{Int})

If a vector encodes to `num` with key generated from kroneckerPregen this method will compute
the vector .% m in the encoded space.
"""
function mod_kronecker(num, m, numVars, kroneckerPregen::Vector{Int})
    return num - div_kronecker(num, m, numVars, kroneckerPregen) * m
end

"""
    div_kronecker(num, m, numVars, kroneckerPregen::Vector{Int})

If a vector encodes to `num` with key generated from kroneckerPregen this method will compute
the vector .÷ m in the encoded space.
"""
function div_kronecker(num, m, numVars, kroneckerPregen::Vector{Int})
    result = 0
    for i in numVars:-1:1
        q, num = divrem(num, kroneckerPregen[i])
        result += (q ÷ m) * kroneckerPregen[i]
    end
    
    return result
end


function matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,d)
    numVars = size(degs,1)
    mons = gen_exp_vec(numVars,d)
    mons = reduce(hcat, mons)

    nMons = size(mons,2)
    nTerms = size(degs,2)

    # Everything needs to be encoded with the same key to make the is_relevant() check work
    maxdeg = p * (d + 1)

    # We want to avoid having to compute te same powers of the encoding key every time we call kronecker()
    kroneckerPregen = zeros(Int, numVars)
    for i in eachindex(kroneckerPregen)
        kroneckerPregen[i] = (maxdeg + 1) ^ (i - 1)
    end

    kron(v) = kronecker_opt(v, numVars, kroneckerPregen)
    div_kron(v, m) = div_kronecker(v, m, numVars, kroneckerPregen)
    mod_kron(v, m) = mod_kronecker(v, m, numVars, kroneckerPregen)

    reverseMons = Dict{Int,Int}()
    encodedMons = zeros(Int, nMons)
    tempmon = zeros(Int, numVars)
    for i in axes(mons, 2)
        for j in eachindex(tempmon)
            tempmon[j] = mons[j, i]
        end
        key = kron(tempmon)
        encodedMons[i] = key
        reverseMons[key] = i
    end

    reverseDegs = Dict{Int,Int}()
    encodedDegs = zeros(Int, nTerms)
    tempdeg = zeros(Int, numVars)
    for i in axes(degs, 2)
        for j in eachindex(tempdeg)
            tempdeg[j] = degs[j, i]
        end
        key = kron(tempdeg)
        encodedDegs[i] = key
        reverseDegs[key] = i
    end
  
    encodedMonsModP = map(x -> mod_kron(x, p), encodedMons)
    encodedDegsModP = map(x -> mod_kron(x, p), encodedDegs)
  
    mons_perm = sortperm(encodedMonsModP)
    degs_perm = sortperm(encodedDegsModP)
  
    # we need to traverse both arrays at once
    # we consider degs to be on the "left"
    left = true
  
    l = 1 # left index
    r = nMons # right index
  
  
    result = zeros(eltype(coefs),nMons,nMons)

    relevant = kron(fill(p - 1, numVars))
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
