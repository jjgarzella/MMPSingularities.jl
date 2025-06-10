
##include("../../src/FrobSplittingInfra.jl")

## BEGIN COPYPASTA
##
## TODO: remove copypasta

#exponent_vectors(poly) = leading_exponent_vector.(terms(poly))

#"""
#Evaluates the element of the Frobenius dual corresponding
#to 'indices' for a polynomial ring (this happens to also be a splitting)
#on the polynomial poly.

#An element of the Frobenius dual is an element of 
#Hom(F_*R, R), i.e. a map F_*R \\to R. 
#Here, R = k[x_1, ..., x_n] is a polynomial ring, where
#k is a field of characteristic p.

#Such maps are in 1-1 correspondence with n-tuples
#of integers mod p, where we have the maps
#\\phi_{indices[1], ..., indices[n]},
#as defined in the literature (see e.g. Ma-Polstra's notes).

#Now, these maps are defined by projecting to the direct sum component,
#which means we forget all but the terms which have exponent < indices[i] mod p
#for the variable x_i, then we subtract the exponent of x_i by indices[i] and divide them all by p.

#An example will be illustrive. Consider the element of the three-variable polynomial
#ring k[x,y,z]: 

#f = x^2y^3z + x^3y^2z^2 + x^8y^7z^7

#and say that p = 5, then 

#\\phi_{2,3,1}(f) = 1, and
#\\phi_{3,2,2}(f) = 1 + xyz

#Note: idk exactly what this is doing when you pass an element of a non-polynomial ring
#... it probably throws an error on account of the use of 'gens'

#Assumptions:

#length(indices) == number of variables in the ring

#indices \\in [0, ..., p-1]

#"""
#function polynomial_frobenius_splitting(p,poly,indices)

#  #coefs = collect(coefficients(poly)) # someday use iterator to make this more efficient?
#  #exp_vecs = exponent_vectors(poly)
#  vars = gens(parent(poly))

#  poly == zero(poly) && return 0
#  length(vars) != length(indices) && begin println("mismatched number of variables"); return end

#  result = zero(poly)
#  for i in 1:length(poly)
#    t = term(poly, i)
#    exp_vec = exponent_vector(t,1)

#    if all((exp_vec .% p) .== indices)


#      new_exp_vec = divexact.(exp_vec .- indices,p) # the difision should be exact by the if statement
      
#      c = coeff(t,1)
#      new_term = c * prod(vars .^ new_exp_vec) 
#      # uses that gens and leading_exponent_vector are using the same variable order

#      result = result + new_term
#    end
#  end

#  result

#end#function


#"""
#   gen_exp_vec(n, d, order)

#Returns all nonnegative integer lists of length n who entires sum to d

#These are the exponent vectors for all the homogeneous monomials of
#degree d, in n variables.

#TODO: give this function @memoize. perhaps for some big examples the
#storage required to store the result isn't worth it. However, for
#tests where we're running many similar examples of medium size,
#I think it'll be nicer.

#INPUTS:
#* "n" -- integer
#* "d" -- integer
#* "order" -- string, monomial ordering, defaulted to lexicographic ordering. Also supports neglex
#"""
##=
#function gen_exp_vec(n, d, order=:lex)
#    @assert (n >= 0) && (d >= 0) "n and d need to be non-negative"
#    result = Vector{Int64}[]
#    #=
#    if d == 0
#        return [1]
#    end
#    =#

#    if n == 1
#        return [[d]]
#    end

#    if order == :lex
#        if d == 1
#            for i in 1:n
#                s = zeros(Int64,n)
#                s[end-i+1] = 1
#                push!(result,s)
#            end
#            return result
#        end

#        for i in 0:d
#            y = gen_exp_vec(n-1,d-i,order)
#            for j in axes(y,1)
#                prepend!(y[j],i)
#            end
#            append!(result,y)
#        end

#    elseif order == :neglex
#        if d == 1
#            for i in 1:n
#                s = zeros(Int64,n)
#                s[i] = 1
#                push!(result,s)
#            end
#            return result
#        end

#        for i in 0:d
#            y = gen_exp_vec(n-1,d-i,order)
#            for j in axes(y,1)
#                prepend!(y[j],i)
#            end
#            prepend!(result,y)
#        end

#    elseif order == :invlex
#        vecs = gen_exp_vec(n,d,:lex)
#        result = reverse.(vecs)
#    else
#        throw(ArgumentError("Unsupported order '$order'"))
#    end

#    result
#end
# =#

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
#    elseif order == :neglex
#        for i in 1:n
#            dtemp = copy(d)
#            k = 1
#            while k <= length(result)
#                if i > 1
#                    dtemp = copy(d)
#                    for j in 1:n
#                        dtemp = dtemp - result[k][j]
#                    end
#                end
#                if i == n && dtemp > 0
#                    result[k][i] = dtemp
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 0
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 1 && i > 1
#                    for j in i:n
#                        result[k+j-i][j] = 1
#                    end
#                    k = k + n - i + 1
#                    continue
#                end
#                dtemp2 = copy(d - dtemp)
#                if i == 1 || dtemp2 == 0
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
#                            result[j+k-1][i] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
#                        dtemp = dtemp - 1
#                    end
#                else
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                            result[j+k-1][i] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                        dtemp = dtemp - 1
#                    end
#                end
#            end
#        end
#    elseif order == :invlex
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
#                    result[length(result)-k][n-i+1] = dtemp
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 0
#                    k = k + 1
#                    continue
#                end
#                if dtemp == 1 && i > 1
#                    for j in i:n
#                        result[length(result)-(k+j-i)][n-j+1] = 1
#                    end
#                    k = k + n - i + 1
#                    continue
#                end
#                dtemp2 = copy(d - dtemp)
#                if i == 1 || dtemp2 == 0
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
#                            result[length(result)-(j+k-1)][n-i+1] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
#                        dtemp = dtemp - 1
#                    end
#                else
#                    while dtemp >= 0
#                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                            result[length(result)-(j+k-1)][n-i+1] = dtemp
#                        end
#                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
#                        dtemp = dtemp - 1
#                    end
#                end
#            end
#        end
#    else
#        throw(ArgumentError("Unsupported order '$order'"))
#    end
#    return result
#end

#"""
#FIXME/DOCUMENTME: exp_vec seems to be a 2d array here?
#"""
#function gen_mon(exp_vec, R, PR)
#    result = []
#    for i in axes(exp_vec,1)
#        B = MPolyBuildCtx(PR)
#        push_term!(B, R(1), exp_vec[i])
#        monomial = finish(B)
#        push!(result,monomial)
#    end
#    result
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

#dim_of_homog_polys(n,d)  = binomial(n+d-1,n-1) # n+d-1 choose n-1

#"""
#Wrapper for polynomial_to_vector

#Converts the homogeneous polynomial poly
#to a vector.

#"""
#function vector(f,d,order=:lex)
#  R = parent(f)
#  n = length(gens(R))

#  F = coefficient_ring(R)
#  f == zero(R) && return zeros(F,dim_of_homog_polys(n,d))
#  @assert d == total_degree(f) "Expect d to be the degree of f"
#  polynomial_to_vector(f, n, F, R,order)
#end

#"""
#Takes a linear operator L on the space
#of homogenous polynomials 
#of degree d and computes 
#the matrix representing it.

#Currently uses lexographical order

#L is a function, which is assumed to be a linear
#endomoprhism on the vector space of homogeneous
#polynomials.

#d is the degree of the homogeneous polynomials.

#R is the base ring.

#"""
#function matrix_of_lin_op(L,d,R,order=:lex)

#  n = length(gens(R))
#  monomials = compute_monomials(n,d,R,order)

#  m = length(monomials) # will be an mxm matrix

#  i = 0

#  matrix = zeros(coefficient_ring(R),m,0)
#  println(monomials)
#  for monomial in monomials
#    evaled = L(monomial)
#    v = vector(evaled,d)
#    matrix = [matrix v]

#    if leading_exponent_vector(monomial) == [14,1,1,0]
#      println(v)
#    end

#    i = i + 1
#    if i % 50 == 0 
#      println("50 rows completed")
#    end
#  end

#  matrix
#end

## END COPYPASTA

# let's set p = q
# phi(a) = a^p

"""
For q = p, otherwise we need to apply
ϕ to the coefficients
"""
function ϕ(p,f)
    res = zero(f)
    vars = gens(parent(f))
    exps = exponent_vectors(f)
    ts = collect(terms(f))
    for i in eachindex(ts)
        mon = prod(vars .^ (p .* exps[i]))
        res += coeff(ts[i],1) * mon
    end
  
    res
end

"""
Again, if p ≠ q then this does not work
"""
function ψ(p,f)
    n = length(gens(parent(f)))
    polynomial_frobenius_splitting(p,f,zeros(Int,n))
end

function A(p,f)
    g -> ψ(p,g * f^(p-1))
end

function matrix_A(p,f)
    R = parent(f)
    d = total_degree(f)
    matrix_of_lin_op(A(p,f),d,R)
end

my_trace(A) = sum(A[i,i] for i in 1:size(A,1))

function first_harvey_trace(p,f)
    my_trace(matrix_A(p,f))
end


function αₛ(λ,s,τ)
    result = 1
    for t in 0:τ-1
        result += binomial(-λ,t) * binomial(λ,s-t)
    end
    result *= (-1)^s 
    result
end


#TODO:
#
#* implement/get multiply then split for A_f
#* implement α_s
#* implement the trace formula
#* implement zeta functions of k3 / cubic fourfold, plug into the formula
#
#* make it fast: implement lemma 3.2 with multiply then split
#
# λ is the precision, r is the order of the finite field extension,
# a is the exponent of p^a = q, and τ is the "fudge factor"



# MARK - zeta function utilities

function zeta_of_Pn(n,q,t)
  res = 1
  for i in 0:n
      res *= 1/(1 - (q^i)*t)
  end
  res
end

function L_poly_ell(fullzeta,q)
    t = gen(parent(fullzeta))
    L = fullzeta
    
    L *= (1-t)
    L *= (1-q*t)

    L
end

function L_poly_K3(fullzeta,q)
    t = gen(parent(fullzeta))
    L = fullzeta
    for i in 0:2
        L *= (1 - (q^i)*t)
    end

    1/L
end


function L_poly_cubicfourfold(fullzeta,q)
    t = gen(parent(fullzeta))
    L = fullzeta
    for i in 0:4
        L *= (1 - (q^i)*t)
    end

    1/L
end
