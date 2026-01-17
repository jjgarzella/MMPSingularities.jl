using Combinatorics
using Oscar

"""
    h: homogeneous degree
    m: non-inclusive upper bound on coeffs
    n: # variables
    k: power
"""
function estimate_upper_bound(h, m, n, k)
    return ((BigInt(m) - 1) * binomial(BigInt(h) + n - 1, n - 1))^k
end

"""
    h: homogeneous degree
    m: non-inclusive upper bound
    n: # variables
    k: power
"""
function upper_bound(h, m, n, k)
    R, x = polynomial_ring(ZZ, ["x$i" for i in 1:n])
    c = ZZ(m) - 1
    # multiexponents(n, h) generates all exponents of length n summing to h
    P = sum(c * prod(x[i]^e[i] for i in 1:n) for e in multiexponents(n, h))
    return maximum(coefficients(P^k))
end

# log2(upper_bound(40, 11, 4, 4)) = 51.8 (g^4 for (11, 4))
# log2(upper_bound(40, 11, 4, 5)) = 68.3 (g^5 for (11, 4)) OVERFLOWS!!!

# log2(upper_bound(48, 13, 4, 4)) = 55.1 (g^4 for (13, 4))

"""
    h1: homogeneous degree of polynomial 1
    m1: non-inclusive upper bound on coefficients of polynomial 1
    h2: homogeneous degree of polynomial 2
    m2: non-inclusive upper bound on coefficients of polynomial 2
    n: # variables
"""
function upper_bound_mul(h1, m1, h2, m2, n)
    R, x = polynomial_ring(ZZ, ["x$i" for i in 1:n])
    c1 = ZZ(m1) - 1
    c2 = ZZ(m2) - 1
    P1 = sum(c1 * prod(x[i]^e[i] for i in 1:n) for e in multiexponents(n, h1))
    P2 = sum(c2 * prod(x[i]^e[i] for i in 1:n) for e in multiexponents(n, h2))
    return maximum(coefficients(P1 * P2))
end

# log2(upper_bound_mul())