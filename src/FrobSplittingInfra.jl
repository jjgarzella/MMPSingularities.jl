#module FrobSplittingInfra
#
#using Oscar
#using Memoize
#using Combinatorics
#
#
#include("griffiths-dwork-construction/Utils.jl")
#
#export polynomial_frobenius_splitting
#export polynomial_frobenius_generator
#export multiply_then_split
#export is_kernel_poly_frob_generator
#export multicomibnations
#export Δ₁, Δ₁l, Fstar_basis
#export dim_of_homog_polys
#export vector, matrix_of_lin_op
#export lift_to_Int64
#export index_of_term_not_in_frobenius_power_CY
#export inPowerOfVariableIdeal
#export isHomog, isFSplit
#export matrix_of_multiply_then_split
#export matrix_of_multiply_then_split_sortmodp



exponent_vectors(poly) = leading_exponent_vector.(terms(poly))

"""
Evaluates the element of the Frobenius dual corresponding
to 'indices' for a polynomial ring (this happens to also be a splitting)
on the polynomial poly.

An element of the Frobenius dual is an element of 
Hom(F_*R, R), i.e. a map F_*R \\to R. 
Here, R = k[x_1, ..., x_n] is a polynomial ring, where
k is a field of characteristic p.

Such maps are in 1-1 correspondence with n-tuples
of integers mod p, where we have the maps
\\phi_{indices[1], ..., indices[n]},
as defined in the literature (see e.g. Ma-Polstra's notes).

Now, these maps are defined by projecting to the direct sum component,
which means we forget all but the terms which have exponent < indices[i] mod p
for the variable x_i, then we subtract the exponent of x_i by indices[i] and divide them all by p.

An example will be illustrive. Consider the element of the three-variable polynomial
ring k[x,y,z]: 

f = x^2y^3z + x^3y^2z^2 + x^8y^7z^7

and say that p = 5, then 

\\phi_{2,3,1}(f) = 1, and
\\phi_{3,2,2}(f) = 1 + xyz

Note: idk exactly what this is doing when you pass an element of a non-polynomial ring
... it probably throws an error on account of the use of 'gens'

Assumptions:

length(indices) == number of variables in the ring

indices \\in [0, ..., p-1]

"""
function polynomial_frobenius_splitting(p,poly,indices)

  #coefs = collect(coefficients(poly)) # someday use iterator to make this more efficient?
  #exp_vecs = exponent_vectors(poly)
  vars = gens(parent(poly))

  poly == zero(poly) && return 0
  length(vars) != length(indices) && begin println("mismatched number of variables"); return end

  result = zero(poly)
  for i in 1:length(poly)
    t = term(poly, i)
    exp_vec = exponent_vector(t,1)

    if all((exp_vec .% p) .== indices)


      new_exp_vec = divexact.(exp_vec .- indices,p) # the difision should be exact by the if statement
      
      c = coeff(t,1)
      new_term = c * prod(vars .^ new_exp_vec) 
      # uses that gens and leading_exponent_vector are using the same variable order

      result = result + new_term
    end
  end

  result

end#function

"""
By a standard argument, e.g. see Ma-Polstra,
Hom(F_*R, R) is generated by a single element,
which may be denoted 
\\phi_{p-1, ..., p-1} 
in the previous notation.
In Kawakami-Takamatsu-Yoshikawa, they denote this by u.

This function evaluates this generator at poly.
"""
function polynomial_frobenius_generator(p,poly) 
  nVars = length(gens(parent(poly)))

  polynomial_frobenius_splitting(p,poly,fill(p-1, nVars))
end#function

"""
Returns true if the polynomial inputted is in the
kernel of the map u that gnerates the dual module
of F_*
"""
function in_kernel_poly_frob_generator(p,poly)
  nVars = length(gens(parent(poly)))
  indices = fill(p-1,nVars)

  for i in 1:length(poly)

    t = term(poly, i)
    exp_vec = exponent_vector(t,1)

    if all((exp_vec .% p) .== indices)

      # this term will not be zero!
      return false
    end
  end

  return true
end#function

# MARK - computing Δ_1

"""
This is not as space-efficient as it could be

Gives the ways of putting the numbers 1...sum(partition)
into length(unique(partition)) boxes where each box
is labeled by a part value and box labeled
by partition[i] has 
count(==(partition[i]),partition) elements

list is assumed to be a list of indices, i.e.
a list of integers with no repeated elements.

partition is assumed to be a list of positive
integers whose which is in increasing or decreasing
order 
(though I think one can get away with each value
being in a contiguous interval of indices, i.e.
[1;1;1;2;2;3] is ok but [1;1;2;1] is not)

The resulting output is an array which has the same
length as partition, whose i'th element is in the
bucket labeled by partition[i]
"""
@memoize function multicombinations(list,partition)
  #TODO: can this be more space-efficient as an iterator??
  # TODO: rewrite this function to allocate memory smarter

  if length(partition) == 0
    return [zeros(eltype(partition),0)]
  end
  
  first = partition[1]
  nFirst = count(==(first),partition)

  result = Vector{Vector{eltype(list)}}()
  for combo in Combinatorics.combinations(list,nFirst)
    withoutfirst = partition[partition .!= first]
    withoutcombo = list[list .∉ (combo,)] # ∉ is \notin 

    for multicombo in multicombinations(withoutcombo,withoutfirst)
      push!(result,[combo; multicombo])
    end
  end

  result
end#function

"""
finds the multicombinations of the set {1, ... n},
where n is the sum of the elements in partition
"""
function multicombinations(partition)
  multicombinations(collect(1:sum(partition)),partition)
end#function

"""
Computes Δ₁(poly) by iterating
through all possible multicombinations 
(cross terms).
"""
function Δ₁(p,poly)
  
  # I don't know how computationally intense these operations are
  #    do them outside any loops just in case.
  terms_iter = terms(poly)
  allterms = collect(terms_iter)
  nTerms = length(allterms)

  res = zero(poly)

  for partition in Combinatorics.partitions(p)
    # the algorithm will still work without this line 
    # because div(1,p) == 0 so coef == 0 below...
    #
    # ...but let's not do the extra work.
    length(partition) == 1 && continue

    distParts = unique(partition)
    nDistParts = length(distParts)
  
    # println("Partition", partition)

    for termlist in Combinatorics.combinations(allterms,length(partition))
      
      coef = div(multinomial(partition...),p)

      # println("Coef: ", coef)

      for multicombo in multicombinations(collect(1:length(partition)),partition)

        # println("Multicombo: ", multicombo)

        newterm = one(poly)
        for i in 1:length(partition)
          exponent = partition[i]
          termind = multicombo[i]
          newterm = newterm * (termlist[termind])^exponent
        end

        # println("Term added: ", newterm)
        res = res + coef * newterm
      end

    end

  end

  res
end#function


"""
Returns Δ₁(poly), computed by
lifting to characteristic zero,
computing the cross terms, and
reducing again
"""
function Δ₁l(p,poly)

  R = parent(poly)

  originallift = map_coefficients(x -> lift(ZZ,x),poly)

  #println("Original Lift: ", originallift)

  ZR = parent(originallift)



  nocrossterms = sum(terms(originallift) .^p)
  withcrossterms = originallift^p

  #println("No cross terms: ", nocrossterms)
  #println("With cross terms: ", withcrossterms)

  crossterms = withcrossterms - nocrossterms

  #println("Just cross terms: ", crossterms)
  #(x,y,z,w) = gens(ZR)
  #println("Interesting coefficient: $(coeff(withcrossterms,x^14*y^27*z^11*w^28))")

  Δlift = map_coefficients(x -> divexact(x,p),crossterms)

  change_coefficient_ring(coefficient_ring(R),Δlift,parent=R)

end#function

"""
If R is the polynomial ring in N variables,
returns a basis for R as an R^p - module,
equivalently returns a basis for
F_*R as an R-module.

"""
function Fstar_basis(p,poly)
  vars = gens(parent(poly))
  n = length(vars)

  numgens = p^n - 1

  generators = zeros(parent(poly),fill(p,n)...)
  for i in CartesianIndices(generators)
    exps = Tuple(i) .- 1
    generators[i] = prod(vars .^ exps)
  end

  vec(generators)
end#function

dim_of_homog_polys(n,d)  = binomial(n+d-1,n-1) # n+d-1 choose n-1

"""
Gives the index of the "critical term", i.e. the one
that determines whether or not something is quasi-F-split
for a CY manifold. This term is the degree (p-1)^n
monomial x_1^(p-1) ... x_n^(p-1)

"""
function index_of_term_not_in_frobenius_power_CY(p,n,order=:lex)
  R, vars = polynomial_ring(GF(p),n)
  
  crit_term = prod(vars .^ (p-1))

  # perhaps assert this has only one element?
  findfirst(vector(crit_term,total_degree(crit_term)) .!= 0)
end

"""
Returns true if the polynomial poly
is in the "frobenius power" \\frak{m}^[p],
where {m} is the ideal of variables of the ring.

"""
function inPowerOfVariableIdeal(p,m,poly)
  # don't need this because exponent_vectors will have 
  # no elements for the zero polynomial
  poly == zero(poly) && return true


  for i in 1:length(poly)
    ev = exponent_vector(poly,i)

    if all(ev .< m)
      #println("Found term not in the Frob power of the maximal ideal: " * string(exponent_vector))
      
      # We not in the power of the maximal ideal, we don't have any
      # powers that are big enough
      return false
    end
  end

  true
end#function

"""
Returns true if the polyonoimal is homogenous.

doesn't seem to be fully implemented in Oscar :(

try 'is_homogeneous', doesn't seem to be implemented 
for fpMPolyRingElem
"""
function isHomog(poly;ofdegree=-1)

  evs = exponent_vectors(poly)

  if ofdegree == -1 
    all(sum.(evs) .== sum(evs[1]))
  else
    all(sum.(evs) .== ofdegree)
  end
end#function


# MARK - calculations of quasi-F-split height in one form or another

"""
Calculates if the hypersuface defined by the 
polynomial poly is F-split

note that p must be prime for this to have mathematical meaning
"""
function isFSplit(p,poly)
  #maybe TODO: check that p is prime

  !inPowerOfVariableIdeal(p,p,poly^(p-1))

end#function


#end#module
