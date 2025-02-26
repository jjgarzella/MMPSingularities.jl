
include("../../src/FrobSplittingInfra.jl")
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
    g -> ψ(g * f^(p-1))
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
