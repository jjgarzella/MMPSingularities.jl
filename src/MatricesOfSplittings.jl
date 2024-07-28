
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


"""
given two exponent vectors degs1 and degs2, 
will the term corresponding to them survive
the frobenius root? 

In other words, to these vectors sum to 
fill(p-1,n)?

degs1 and degs2 are assumed to be reduced mod p already

This function is meant to be used in tight inner loops,
so if there's any way to make it faster that would be good.

"""
function terms_are_relevant(p,degs1,degs2)
  n = length(degs1)
  
  #println("Testing $degs1 vs. $degs2 mod $p")
  relevant = true
  for k in 1:n

    #prod_exp_vec_mod_p[k] = degs[tInd,k] + mons[i][k] % p
    
    #println("Coord $k: Does $(degs1[k]) == $(degs2[k]) ???")
    if degs1[k] + degs2[k] != p-1
      #println("It does not! Test fail!")
      relevant = false
    end
  end
  #println("Tested to $relevant")

  relevant
end

"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This version uses Static vectors
"""
function matrix_of_multiply_then_split_sortmodp(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  #mons = reduce(vcat,transpose.(mons))

  nMons = length(mons)
  nTerms = size(degs,1)
    
  # Preprocessing
  mons = SVector{n}.(mons)
  degs = SVector{n}.(eachrow(degs))

  #println(typeof(degs))
  #println(typeof(mons))
  # convert to static array
  
  #degs = SMatrix{nTerms,n}(degs)
  #mons = SMatrix{nMons,n}(mons)

  # calculate going from exponent vector to index
  reverseMons = Dict(mons[i] => i for i in 1:size(mons,1))
  reverseDegs = Dict(degs[i] => i for i in 1:size(degs,1))

  
  degs_modp = map(v -> v .% p, degs)
  mons_modp = map(v -> v .% p, mons)

  degs_perm = sortperm(degs_modp)
  mons_perm = sortperm(mons_modp)

  # this uses a view to imptove performance later
  #degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms))
  #mons_perm = sortperm(view.(Ref(mons_modp),1:nMons))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))
  #
  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs[1]),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = 0
  #col = 0
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))

  while l ≤ nTerms && 1 ≤ r
    mon_modp = mons_modp[mons_perm[r]]
    term_modp = degs_modp[degs_perm[l]]
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!
      # println("Found match: ($l,$r)")
     

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      while l + nMatches ≤ nTerms && degs_modp[degs_perm[l+nMatches]] == term_modp
        nMatches = nMatches + 1
      end

      # loop through all monomials and process each one
      while 1 ≤ r && mons_modp[mons_perm[r]][:] == mon_modp

        mon = mons[mons_perm[r]]
        for ll = l:(l + nMatches - 1)
          term = degs[degs_perm[ll]]
          #print("found match ($ll,$r), ")
          #print("accessing term at ($(degs_perm[ll]),$(degs_perm[r])), ")
          #print("term $term multiplies with monomial $mon, ")
          
          # multiply then split the terms
          #newterm = div.(mon .+ term .- fill(p-1,n), p)
          for i in 1:n
              newterm[i] = div(mon[i] + term[i] - (p-1), p)
          end

          newcoefind = reverseDegs[term] 
          newcoef = coefs[newcoefind]
          row = reverseMons[newterm]
          col = reverseMons[mon]
          #println("setting matrix element ($row,$col)")
          result[row,col] += newcoef
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end

"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This version uses tuples
"""
function matrix_of_multiply_then_split_sortmodp_tuple(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  #mons = reduce(vcat,transpose.(mons))

  nMons = length(mons)
  nTerms = size(degs,1)
    
  # Preprocessing
  mons = Tuple.(mons)
  degs = Tuple.(eachrow(degs))

  #println(typeof(degs))
  #println(typeof(mons))
  # convert to static array
  
  #degs = SMatrix{nTerms,n}(degs)
  #mons = SMatrix{nMons,n}(mons)

  # calculate going from exponent vector to index
  reverseMons = Dict(mons[i] => i for i in 1:size(mons,1))
  reverseDegs = Dict(degs[i] => i for i in 1:size(degs,1))

  
  degs_modp = map(v -> v .% p, degs)
  mons_modp = map(v -> v .% p, mons)

  degs_perm = sortperm(degs_modp)
  mons_perm = sortperm(mons_modp)

  # this uses a view to imptove performance later
  #degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms))
  #mons_perm = sortperm(view.(Ref(mons_modp),1:nMons))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))
  #
  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs[1]),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = 0
  #col = 0
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))

  while l ≤ nTerms && 1 ≤ r
    mon_modp = mons_modp[mons_perm[r]]
    term_modp = degs_modp[degs_perm[l]]
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!
      # println("Found match: ($l,$r)")
     

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      while l + nMatches ≤ nTerms && degs_modp[degs_perm[l+nMatches]] == term_modp
        nMatches = nMatches + 1
      end

      # loop through all monomials and process each one
      while 1 ≤ r && mons_modp[mons_perm[r]][:] == mon_modp

        mon = mons[mons_perm[r]]
        for ll = l:(l + nMatches - 1)
          term = degs[degs_perm[ll]]
          #print("found match ($ll,$r), ")
          #print("accessing term at ($(degs_perm[ll]),$(degs_perm[r])), ")
          #print("term $term multiplies with monomial $mon, ")
          
          # multiply then split the terms
          #newterm = tuple((@. div(mon + term - (fill(p-1,n),)), p))
          for i in 1:n
              newterm[i] = div(mon[i] + term[i] - (p-1), p)
          end

          newcoefind = reverseDegs[term] 
          newcoef = coefs[newcoefind]
          row = reverseMons[tuple(newterm...)]
          col = reverseMons[mon]
          #println("setting matrix element ($row,$col)")
          result[row,col] += newcoef
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end



"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This version uses a dictionary of vectors, which allocates when it hashes.
That's causing the majority of the performance issues.

"""
function matrix_of_multiply_then_split_sortmodp_dict(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))

  nMons = size(mons,1)
  nTerms = size(degs,1)
    
  # Preprocessing
  reverseMons = Dict(mons[i,:] => i for i in 1:size(mons,1))
  reverseDegs = Dict(degs[i,:] => i for i in 1:size(degs,1))

#  println(reverseDegs[[24,24,24,8]])

  #println(reverseMons)

  degs_modp = degs .% p
  mons_modp = mons .% p

  # this uses a view to imptove performance later
  degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms,:))
  mons_perm = sortperm(view.(Ref(mons_modp),1:nMons,:))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))

  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = 0
  #col = 0
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))

  while l ≤ nTerms && 1 ≤ r
    mon_modp = @view mons_modp[mons_perm[r],:]
    term_modp = @view degs_modp[degs_perm[l],:]
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!
     

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      while l + nMatches ≤ nTerms && @view(degs_modp[degs_perm[l+nMatches],:]) == term_modp
        nMatches = nMatches + 1
      end

      # loop through all monomials and process each one
      while 1 ≤ r && @view(mons_modp[mons_perm[r],:]) == mon_modp

        mon = @view mons[mons_perm[r],:]
        for ll = l:(l + nMatches - 1)
          term = @view degs[degs_perm[ll],:]
          #print("accessing term at ($(degs_perm[ll]),$(degs_perm[r])), ")
          #print("term $term multiplies with monomial $mon, ")
          
          # multiply then split the terms
          #newterm = div.(mon .+ term .- fill(p-1,n), p)
          for i in 1:n
              newterm[i] = div(mon[i] + term[i] - (p-1), p)
          end

          newcoefind = reverseDegs[term] 
          newcoef = coefs[newcoefind]
          row = reverseMons[newterm]
          col = reverseMons[mon]
          #if row == 565
          #  print("found match ($ll,$r), ")
          #  print("coef ind $newcoefind, from term $term,")
          #  println("matrix elt ($row,$col)")
          #end
          result[row,col] += newcoef
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end


function kronecker(vec,d,n)#,returntype=Int64)
  s = 0#zero(returntype)
  #println("$n")
  for i = 1:n
    s = s + vec[i]*(d+1)^i
  end

  s
end

"""
we don't actually need this right now
"""
function undo_kronecker(p,d,vec)

end

"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This one uses Kronecker substitution
"""
function matrix_of_multiply_then_split_sortmodp_kronecker(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))
  #elementtype = eltype(degs)
  #TODO: perhaps if mons and degs are HybridArrays from HybridArrays.jl
  #  then we won't need to do all of the crazy tricks instead of slicing

  nMons = size(mons,1)
  nTerms = size(degs,1)
    
  monkron(v) = kronecker(v,d,n)
  degkron(v) = kronecker(v,p*d,n)

  # the following is processing it as an array of pairs instead of making the paris be the keys and values
  # Preprocessing
  reverseMons = Dict{Int,Int}()
  tempmon = zeros(eltype(mons),n)
  for i in 1:size(mons,1)
    for j = 1:n # do a for loop so we don't allocate
      tempmon[j] = mons[i,j]
    end
    key = monkron(tempmon)
    reverseMons[key] = i
  end

  reverseDegs = Dict{Int,Int}()
  tempdeg = zeros(eltype(degs),n)
  for i in 1:size(degs,1)
    for j = 1:n # do a for loop so we don't allocate
      tempdeg[j] = degs[i,j]
    end
    key = degkron(tempdeg)
    reverseDegs[key] = i
  end
  #reverseMons = Dict(monkron(mons[i,:]) => i for i in 1:size(mons,1))
  #reverseDegs = Dict(degkron(degs[i,:]) => i for i in 1:size(degs,1))
  # allocation from the above line will eventually be a bottleneck
  # so we replace with for-loop-based code that does not slice

  #jprintln(typeof(reverseMons))
  ## DID NOT REPRODUCE
  ##
  #iii = 1002320
# # mymon = mons[10,:]
  ##kmymon = monkron(mymon)
  #@time b = reverseMons[iii]
  #println(b)
  #error() 

  degs_modp = degs .% p
  mons_modp = mons .% p

  # this uses a view to imptove performance later
  degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms,:))
  mons_perm = sortperm(view.(Ref(mons_modp),1:nMons,:))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))

  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  # PRE-ALLOCATION

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = Ref{Int64}(0)
  #col = Ref{Int64}(0)
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))
  
  # this function is meant to be inlined
  function setslice!(target,source,ind)
    for i = 1:n
      target[i] = source[ind,i]
    end
  end

  mon_modp = zeros(eltype(degs),n)
  term_modp = zeros(eltype(degs),n)
  cmp_term = zeros(eltype(degs),n)
  cmp_mon = zeros(eltype(degs),n)
  term = zeros(eltype(degs),n)
  mon = zeros(eltype(degs),n)

  while l ≤ nTerms && 1 ≤ r
    #mon_modp = @view mons_modp[mons_perm[r],:]
    setslice!(mon_modp,mons_modp,mons_perm[r])
    #term_modp = @view degs_modp[degs_perm[l],:]
    setslice!(term_modp,degs_modp,degs_perm[l])
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      #cmp_term = @view(degs_modp[degs_perm[l+nMatches],:])
      setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      while l + nMatches ≤ nTerms && cmp_term == term_modp
        nMatches = nMatches + 1
        l + nMatches ≤ nTerms && setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      end

      # loop through all monomials and process each one

      #cmp_term = @view(mons_modp[mons_perm[r],:])
      setslice!(cmp_mon,mons_modp,mons_perm[r])

      while 1 ≤ r && cmp_mon == mon_modp

        #mon = @view mons[mons_perm[r],:]
        setslice!(mon,mons,mons_perm[r])
        for ll = l:(l + nMatches - 1)
            #term = @view degs[degs_perm[ll],:]
            setslice!(term,degs,degs_perm[ll])

            #print("found match ($ll,$r), ")
            #print("accessing at ($(degs_perm[ll]),$(degs_perm[r])), ")
            #println("term $term, monomial $mon, ")
            
            # multiply then split the terms
            #newterm = div.(mon .+ term .- fill(p-1,n), p)
            for i in 1:n
                newterm[i] = div(mon[i] + term[i] - (p-1), p)
            end

            #println("kterm")
            kterm = degkron(term)
            #println((kterm))
            #println("newcoefind")
            newcoefind = reverseDegs[kterm] 
            #println(typeof(newcoefind))
            newcoef = coefs[newcoefind]

            #println("knewterm")
            knewterm = monkron(newterm)
            #println((knewterm))
            #println("row")
            row = reverseMons[knewterm]
            #println(typeof(row))
            #println("kmon")
            kmon = monkron(mon)
            #println(typeof(kmon))
            #println("col")
            col = reverseMons[kmon]
            #println(typeof(col))
            #println("setting matrix element ($row,$col)")
            result[row,col] += newcoef
            #error()
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        1 ≤ r && setslice!(cmp_mon,mons_modp,mons_perm[r])
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end

using Dictionaries

"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This one uses Kronecker substitution
"""
function matrix_of_multiply_then_split_sortmodp_kronecker_dictionary(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))
  #elementtype = eltype(degs)

  nMons = size(mons,1)
  nTerms = size(degs,1)
    
  monkron(v) = kronecker(v,d,n)
  degkron(v) = kronecker(v,p*d,n)

  # the following is processing it as an array of pairs instead of making the paris be the keys and values

  # Preprocessing
  reverseMons = Dictionary{Int,Int}()
  tempmon = zeros(eltype(mons),n)
  for i in 1:size(mons,1)
    for j = 1:n # do a for loop so we don't allocate
      tempmon[j] = mons[i,j]
    end
    key = monkron(tempmon)
    insert!(reverseMons,key,i)
  end

  reverseDegs = Dictionary{Int,Int}()
  tempdeg = zeros(eltype(degs),n)
  for i in 1:size(degs,1)
    for j = 1:n # do a for loop so we don't allocate
      tempdeg[j] = degs[i,j]
    end
    key = degkron(tempdeg)
    insert!(reverseDegs,key,i)
  end
  #reverseMons = Dict(monkron(mons[i,:]) => i for i in 1:size(mons,1))
  #reverseDegs = Dict(degkron(degs[i,:]) => i for i in 1:size(degs,1))
  #TODO: allocation from the above line will eventually be a bottleneck
  #      please replace with for-loop-based code that does not slice

  #println(typeof(reverseMons))
  ## DID NOT REPRODUCE
  ##
  #iii = 1002320
# # mymon = mons[10,:]
  ##kmymon = monkron(mymon)
  #@time b = reverseMons[iii]
  #println(b)
  #error() 

  degs_modp = degs .% p
  mons_modp = mons .% p

  # this uses a view to imptove performance later
  degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms,:))
  mons_perm = sortperm(view.(Ref(mons_modp),1:nMons,:))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))

  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  # PRE-ALLOCATION

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = Ref{Int64}(0)
  #col = Ref{Int64}(0)
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))
  
  # this function is meant to be inlined
  function setslice!(target,source,ind)
    for i = 1:n
      target[i] = source[ind,i]
    end
  end

  mon_modp = zeros(eltype(degs),n)
  term_modp = zeros(eltype(degs),n)
  cmp_term = zeros(eltype(degs),n)
  cmp_mon = zeros(eltype(degs),n)
  term = zeros(eltype(degs),n)
  mon = zeros(eltype(degs),n)

  while l ≤ nTerms && 1 ≤ r
    #mon_modp = @view mons_modp[mons_perm[r],:]
    setslice!(mon_modp,mons_modp,mons_perm[r])
    #term_modp = @view degs_modp[degs_perm[l],:]
    setslice!(term_modp,degs_modp,degs_perm[l])
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      #cmp_term = @view(degs_modp[degs_perm[l+nMatches],:])
      setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      while l + nMatches ≤ nTerms && cmp_term == term_modp
        nMatches = nMatches + 1
        l + nMatches ≤ nTerms && setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      end

      # loop through all monomials and process each one

      #cmp_term = @view(mons_modp[mons_perm[r],:])
      setslice!(cmp_mon,mons_modp,mons_perm[r])

      while 1 ≤ r && cmp_mon == mon_modp

        #mon = @view mons[mons_perm[r],:]
        setslice!(mon,mons,mons_perm[r])
        for ll = l:(l + nMatches - 1)
            #term = @view degs[degs_perm[ll],:]
            setslice!(term,degs,degs_perm[ll])

            #print("found match ($ll,$r), ")
            #print("accessing at ($(degs_perm[ll]),$(degs_perm[r])), ")
            #println("term $term, monomial $mon, ")
            
            # multiply then split the terms
            #newterm = div.(mon .+ term .- fill(p-1,n), p)
            for i in 1:n
                newterm[i] = div(mon[i] + term[i] - (p-1), p)
            end

            #println("kterm")
            kterm = degkron(term)
            #println((kterm))
            #println("newcoefind")
            newcoefind = reverseDegs[kterm] 
            #println(typeof(newcoefind))
            newcoef = coefs[newcoefind]

            #println("knewterm")
            knewterm = monkron(newterm)
            #println((knewterm))
            #println("row")
            row = reverseMons[knewterm]
            #println(typeof(row))
            #println("kmon")
            kmon = monkron(mon)
            #println(typeof(kmon))
            #println("col")
            col = reverseMons[kmon]
            #println(typeof(col))
            #println("setting matrix element ($row,$col)")
            result[row,col] += newcoef
            #error()
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        1 ≤ r && setslice!(cmp_mon,mons_modp,mons_perm[r])
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end

struct FastHashInt{T<:Integer}; i::T; end

Base.:(==)(x::FastHashInt, y::FastHashInt) = x.i == y.i
Base.hash(x::FastHashInt, h::UInt) = xor(UInt(x.i), h)

"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This one uses Kronecker substitution and fast hash ints

"""
function matrix_of_multiply_then_split_sortmodp_kronecker_fasthashint(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))
  #elementtype = eltype(degs)

  nMons = size(mons,1)
  nTerms = size(degs,1)
    
  monkron(v) = kronecker(v,d,n)
  degkron(v) = kronecker(v,p*d,n)

  # Preprocessing
  #reverseMons = Dict{FastHashInt{Int64},Int64}()
  #for i in 1:size(mons,1)
  #  reverseMons[FastHashInt(monkron(

  reverseMons = Dict(FastHashInt(monkron(mons[i,:])) => i for i in 1:size(mons,1))
  reverseDegs = Dict(FastHashInt(degkron(degs[i,:])) => i for i in 1:size(degs,1))
  #TODO: allocation from the above line will eventually be a bottleneck
  #      please replace with for-loop-based code that does not slice

  # DID NOT REPRODUCE
  #
  #iii = 1002320
#  mymon = mons[10,:]
  #kmymon = monkron(mymon)
  #@time b = reverseMons[iii]
  #println(b)
  #error() 

  degs_modp = degs .% p
  mons_modp = mons .% p

  # this uses a view to imptove performance later
  degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms,:))
  mons_perm = sortperm(view.(Ref(mons_modp),1:nMons,:))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))

  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  # PRE-ALLOCATION

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = Ref{Int64}(0)
  #col = Ref{Int64}(0)
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))
  
  # this function is meant to be inlined
  function setslice!(target,source,ind)
    for i = 1:n
      target[i] = source[ind,i]
    end
  end

  mon_modp = zeros(eltype(degs),n)
  term_modp = zeros(eltype(degs),n)
  cmp_term = zeros(eltype(degs),n)
  cmp_mon = zeros(eltype(degs),n)
  term = zeros(eltype(degs),n)
  mon = zeros(eltype(degs),n)

  while l ≤ nTerms && 1 ≤ r
    #mon_modp = @view mons_modp[mons_perm[r],:]
    setslice!(mon_modp,mons_modp,mons_perm[r])
    #term_modp = @view degs_modp[degs_perm[l],:]
    setslice!(term_modp,degs_modp,degs_perm[l])
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      #cmp_term = @view(degs_modp[degs_perm[l+nMatches],:])
      setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      while l + nMatches ≤ nTerms && cmp_term == term_modp
        nMatches = nMatches + 1
        l + nMatches ≤ nTerms && setslice!(cmp_term,degs_modp,degs_perm[l+nMatches])
      end

      # loop through all monomials and process each one

      #cmp_term = @view(mons_modp[mons_perm[r],:])
      setslice!(cmp_mon,mons_modp,mons_perm[r])

      while 1 ≤ r && cmp_mon == mon_modp

        #mon = @view mons[mons_perm[r],:]
        setslice!(mon,mons,mons_perm[r])
        for ll = l:(l + nMatches - 1)
            #term = @view degs[degs_perm[ll],:]
            setslice!(term,degs,degs_perm[ll])

            #print("found match ($ll,$r), ")
            #print("accessing at ($(degs_perm[ll]),$(degs_perm[r])), ")
            #println("term $term, monomial $mon, ")
            
            # multiply then split the terms
            #newterm = div.(mon .+ term .- fill(p-1,n), p)
            for i in 1:n
                newterm[i] = div(mon[i] + term[i] - (p-1), p)
            end

            println("kterm")
            @time kterm = FastHashInt(degkron(term))
            #println(typeof(kterm))
            println("newcoefind")
            @time newcoefind = reverseDegs[kterm] 
            #println(typeof(newcoefind))
            newcoef = coefs[newcoefind]

            println("knewterm")
            @time knewterm = FastHashInt(monkron(newterm))
            #println((knewterm))
            println("row")
            @time row = reverseMons[knewterm]
            #println(typeof(row))
            println("kmon")
            @time kmon = FastHashInt(monkron(mon))
            #println(typeof(kmon))
            println("col")
            @time col = reverseMons[kmon]
            #println(typeof(col))
            #println("setting matrix element ($row,$col)")
            result[row,col] += newcoef
            error()
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        1 ≤ r && setslice!(cmp_mon,mons_modp,mons_perm[r])
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end


"""
Finds the matrix of multiplying by the polynomial with 
coefficients coefs and degrees degs

Performance TODO list:
 * try to ensure access inside the loops is regular (this will be important for gpu version)
 * In terms of allocations, the bottleneck seems to be allocating integers,
     and this seems to have to do with the way the dictionaries are used.
     Probably best to reengineer this in a more performant way.
 * Is it best to make the inner loop into it's own function? 
     Should try this out and see if Julia speeds up

For now this part isn't the bottleneck anymore, but I'll work on it later

This one uses Kronecker substitution

It doesn't give correct answers right now.
"""
function matrix_of_multiply_then_split_sortmodp_kronecker_noslice(p,coefs,degs,d)

  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))

  nMons = size(mons,1)
  nTerms = size(degs,1)

  monkron(v) = kronecker(v,d,n)
  degkron(v) = kronecker(v,p*d,n) # d is (p-1) * original d of the polynomail
    
  # Preprocessing
  reverseMons = Dict(monkron(mons[i,:]) => i for i in 1:size(mons,1))
  reverseDegs = Dict(degkron(degs[i,:]) => i for i in 1:size(degs,1))

  #println()
  #println(reverseDegs[793424])
  #println(kron([24,24,24,8]))
  #println(reverseDegs[kron([24,24,24,8])])
  #println()
  #println(reverseMons)

  degs_modp = degs .% p
  mons_modp = mons .% p

  # this uses a view to imptove performance later
  degs_perm = sortperm(view.(Ref(degs_modp),1:nTerms,:))
  mons_perm = sortperm(view.(Ref(mons_modp),1:nMons,:))
  #degs_perm = sortperm(collect(eachrow(degs_modp)))
  #mons_perm = sortperm(collect(eachrow(mons_modp)))

  #println("mons: $mons")
  #println("mons sorted: $(mons[mons_perm,:])")
  #println("mons sorted mod p: $(mons_modp[mons_perm,:])")
  #println("delta_1 sorted mod p: $(degs_modp[degs_perm,:])")

  # we need to traverse both arrays at once
  # we consider degs to be on the "left"
  left = true

  l = 1 # left index
  r = nMons # right index

  result = zeros(eltype(coefs),nMons,nMons)

  # used in the loop, but allocate here once
  newterm = zeros(eltype(degs),n)

  # Note: preallocating these things doesn't seem to change performance
  #row = 0
  #col = 0
  #newcoefind = 0
  #newcoef = zero(eltype(coefs))

  while l ≤ nTerms && 1 ≤ r
    mon_modp = @view mons_modp[mons_perm[r],:]
    term_modp = @view degs_modp[degs_perm[l],:]
    if terms_are_relevant(p,mon_modp,term_modp)
      # we have a match!
     

      # Short preprocessing step: how many terms in degs have this exponent vector?
      nMatches = 1
      while l + nMatches ≤ nTerms && @view(degs_modp[degs_perm[l+nMatches],:]) == term_modp
        nMatches = nMatches + 1
      end

      # loop through all monomials and process each one
      while 1 ≤ r && @view(mons_modp[mons_perm[r],:]) == mon_modp

        mon = @view mons[mons_perm[r],:]
        for ll = l:(l + nMatches - 1)
          term = @view degs[degs_perm[ll],:]
          #print("found match ($ll,$r), ")
          #print("accessing term at ($(degs_perm[ll]),$(degs_perm[r])), ")
          #print("term $term multiplies with monomial $mon, ")
          
          # multiply then split the terms
          #newterm = div.(mon .+ term .- fill(p-1,n), p)
          for i in 1:n
              newterm[i] = div(mon[i] + term[i] - (p-1), p)
          end

          kterm = degkron(term)
          newcoefind = reverseDegs[kterm] 
          newcoef = coefs[newcoefind]

          knewterm = monkron(newterm)
          row = reverseMons[knewterm]
          kmon = monkron(mon)
          col = reverseMons[kmon]

          #if row == 565
          #  print("found match ($ll,$r), ")
          #  print("coef ind $newcoefind, from term $term,")
          #  println("matrix elt ($row,$col)")
          #end
          #print("found match ($ll,$r), ")
          #print("coef ind $newcoefind, ")
          #println("matrix elt ($row,$col)")

          #println("setting matrix element ($row,$col)")
          result[row,col] += newcoef
        end

        r = r - 1
        left = true
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
      end

      l = l + nMatches - 1
    else
      if left
        l = l + 1
        left = false
      else
        r = r - 1
        #0 < r && println("$r, true monomial: $(mons_perm[r])")
        left = true
      end
    end
  end

  result
end
 


"""
Computes the matrix of the 
linear operator of multiplying
by the polynomnial f with coefficients
coefs and degrees degs and then applying 
polynomial_frobenius_generator
on the vector space of homogeneous polynomials
of degree d

coefs - vector of coefficients
degs - 2d array of exponent vectors
"""
function matrix_of_multiply_then_split(p,coefs,degs,d)
  n = size(degs,2)
  mons = gen_exp_vec(n,d)
  mons = reduce(vcat,transpose.(mons))
  nMons = size(mons,1)

  reverseDict = Dict(mons[i] => i for i in 1:length(mons))

  result = zeros(nMons,nMons)
  
  #p_minus_ones = fill(n,p-1)
  #prod_exp_vec_mod_p = zeros(n)

  for i in 1:nMons
    # compute column i
    for tInd in 1:size(degs,1)

      relevant = true
      for k in 1:n
        #prod_exp_vec_mod_p[k] = degs[tInd,k] + mons[i][k] % p
        if degs[tInd,k] + mons[i,k] % p != p-1
          relevant = false
        end
      end

      if relevant
      #if all((degs[tInd,:] .+ mons[i]) .% p .== fill(n,p-1))
        # this is a relevant term
        exv_in_prod = degs[tInd,:] .+ mons[i,:]
        new_exv = divexact.(exv_in_prod .- fill(n,p-1),p)
        row = reverseDict[new_exv]
        result[row,i] += coefs[tInd]

      end
    end
  end

  result
end

"""
Wrapper for polynomial_to_vector

Converts the homogeneous polynomial poly
to a vector.

"""
function vector(f,d,order=:lex)
  R = parent(f)
  n = length(gens(R))

  f == zero(R) && return zeros(R,dim_of_homog_polys(n,d))
  @assert d == total_degree(f) "Expect d to be the degree of f"
  F = coefficient_ring(R)
  polynomial_to_vector(f, n, F, R,order)
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

  matrix = zeros(R,m,0)
  for monomial in monomials
    evaled = L(monomial)
    v = vector(evaled,d)
    matrix = [matrix v]

    i = i + 1
    if i % 50 == 0 
      println("50 rows completed")
    end
  end
    
  matrix
end

"""
Lifts a matrix with entries in GF(p) to ZZ and converts the entries
to Julia integers
"""
lift_to_Int64(matrix) = Int64.(map(x -> lift(ZZ,x), matrix))

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
    mod_kronecker(num, m, numVars, kroneckerPregen::Vector{Int})

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
    numVars = size(degs,2)
    mons = gen_exp_vec(numVars,d)
    mons = reduce(vcat,transpose.(mons))

    nMons = size(mons,1)
    nTerms = size(degs,1)

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
    for i in axes(mons, 1)
        for j in eachindex(tempmon)
            tempmon[j] = mons[i, j]
        end
        key = kron(tempmon)
        encodedMons[i] = key
        reverseMons[key] = i
    end

    reverseDegs = Dict{Int,Int}()
    encodedDegs = zeros(Int, nTerms)
    tempdeg = zeros(Int, numVars)
    for i in axes(degs, 1)
        for j in eachindex(tempdeg)
            tempdeg[j] = degs[i, j]
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
            cmpTerm = encodedDegsModP[degs_perm[l + nMatches - 1]]

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
            
            if 1 ≤ r
                cmpMon = encodedMonsModP[mons_perm[r]]
            end
            # somehow this is erroring for me - Alex
            # (1 ≤ r) && cmpMon = encodedMonsModP[mons_perm[r]]
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