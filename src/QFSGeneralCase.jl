
"""
Calculates the quasi-F-split height
ht(f) of
f if ht(f) ≤ cutoff. 
Otherwise, 
return a value that is bigger than b.

"""
function quasiFSplitHeight(f,cutoff,pregen=nothing)
    n = length(gens(parent(f)))
    p = characteristic(parent(f))
    d = total_degree(f)

    #TODO: if the base field is not F_p, 
    #then throw an error

    if n == d
        quasiFSplitHeight_CY(p,f,cutoff,pregen)
    else
        #TODO: choose algorithm based
        #on how big the prime is
        quasiFSplitHeight_lift(p,f,cutoff)
    end
end

#TODO: retest example 7.8

"""
Computes the quasi-F-split height, if it is less
than or equal to cutoff, otherwise returns
a value which is strictly bigger than cutoff.

Uses lifting and the CPU to compute delta_1
"""
function quasiFSplitHeight_lift(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  fpminus1 = f^(p-1)
  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  Fstar_gens = Fstar_basis(p,f)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2

  KTY_ideal_generators = [fpminus1] # we don't actually use the initial value for anything but clarity
  KTY_pullback_generators = fpminus1 .* Fstar_gens

  while n ≤ cutoff

    # Step 2. Remove things not in the kernel of u

    #display(KTY_pullback_generators)
    println("Testing for height $n...")


    for i in eachindex(KTY_pullback_generators)
      g = KTY_pullback_generators[i]

      if !in_kernel_poly_frob_generator(p,g)
        # perhaps not the most efficient, but it will work
        KTY_pullback_generators[i] = zero(g)

        #println("Removed generator $i")

      end
    end

    # Step 3. Apply θFstar
    KTY_ideal_generators = θFstar.(KTY_pullback_generators)

    #TODO: Step 3.5. add fpminus1 and take minimal generating set

    # Step 4. Check whether the sequence of ideals terminatres here

    allzero = true

    for generator in KTY_ideal_generators

      if generator == zero(poly)
        continue
      else
        allzero = false
      end

      if !inPowerOfVariableIdeal(p,p,generator)
        # We are quasi-F split of height n! Yay!!
        return n
      end

    end

    allzero && return cutoff + 2 # the chain terminated early, provable infinity

    # We are not n-quasi-F-split, so check if we are n+1-quasi-F-split

    n = n + 1

    cutoff < n && continue # pretty much a break statement

    # Step 1. Compute the generators of the pullback.

    #TODO: remove once step 3.5 is implemented
    KTY_ideal_generators = [KTY_ideal_generators; fpminus1]

    KTY_pullback_generators = typeof(f)[]

    ##TODO: KTY_pullback_generators = zeros(parent(f),length(

    for gen in KTY_ideal_generators
      KTY_pullback_generators = [KTY_pullback_generators; gen .* Fstar_gens]
    end

  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear

end

"""
Takes the array fs and splits
it into chunks by the component of F_star in which it lies

"""
function Fstar_components(fs)

    p = characteristic(base_ring(fs[1]))
    
    vars = gens(parent(fs[1]))

    result = Dict()
    for f in fs

        f == 0 && continue

        component = get_common_variable_factors(f) .% p
        varterm = prod(vars .^ component)

        if haskey(result, component)
            push!(result[component],f / varterm)
        else
            result[component] = [f / varterm]
        end
    end

    result
end

"""
This one takes the minimal number of generators of each 
direct sum component of the module F_* R
"""
function quasiFSplitHeight_lift_mingens(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  R = parent(f)
  gR, _ = grade(R)

  f = gR(f)

  vars = gens(gR)

  fpminus1 = f^(p-1)
  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  Fstar_gens = Fstar_basis(p,f)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2

  KTY_ideal_generators = [fpminus1] # we don't actually use the initial value for anything but clarity
  KTY_pullback_generators = fpminus1 .* Fstar_gens


  while n ≤ cutoff

    # Step 2. Remove things not in the kernel of u
    println("Testing for height $n...")

    #display(KTY_pullback_generators)

    for i in eachindex(KTY_pullback_generators)
      g = KTY_pullback_generators[i]

      if !in_kernel_poly_frob_generator(p,g)
        # perhaps not the most efficient, but it will work
        KTY_pullback_generators[i] = zero(g)

        #println("Removed generator $i")

      end
    end

    

    # println("$(length(KTY_pullback_generators)) generators of the pullback")
    # display(KTY_pullback_generators)
    # I = ideal(KTY_pullback_generators[:])
    # println(I)
    # KTY_pullback_mingens = minimal_generating_set(I)
    # # this is the wrong thing
    # println(KTY_pullback_mingens)
    # println("$(length(KTY_pullback_mingens)) minimal generators")

    comps = Fstar_components(KTY_pullback_generators) 

    KTY_pullback_mingens = []

    for (component, hs) in comps

        I = ideal(hs)
        mingens = minimal_generating_set(I)
        varterm = prod(vars .^ component)
        newgens = varterm .* mingens

        i_gens = length(gens(I))
        min_i_gens = length(mingens)
        difference = i_gens - min_i_gens 
        if 1 ≤ difference
            println("Saved $difference generators")
        end

        for newgen in newgens
            push!(KTY_pullback_mingens, newgen)
        end
    end

    println("Starting to apply ΘFstar")
            
    # Step 3. Apply θFstar
    # KTY_ideal_generators = θFstar.(KTY_pullback_mingens)

    KTY_ideal_generators = similar(KTY_pullback_mingens)
    println(length(KTY_ideal_generators))
    for i in eachindex(KTY_ideal_generators)
        @time KTY_ideal_generators[i] = θFstar(KTY_pullback_mingens[i])
    end


    #TODO: Step 3.5. add fpminus1 and take minimal generating set

    # Step 4. Check whether the sequence of ideals terminatres here

    allzero = true

    for generator in KTY_ideal_generators

      if generator == zero(poly)
        continue
      else
        allzero = false
      end

      if !inPowerOfVariableIdeal(p,p,generator)
        # We are quasi-F split of height n! Yay!!
        return n
      end

    end

    allzero && return cutoff + 2 # the chain terminated early, provable infinity

    # We are not n-quasi-F-split, so check if we are n+1-quasi-F-split

    n = n + 1

    cutoff < n && continue # pretty much a break statement

    # Step 1. Compute the generators of the pullback.

    #TODO: remove once step 3.5 is implemented
    KTY_ideal_generators = [KTY_ideal_generators; fpminus1]

    KTY_pullback_generators = typeof(f)[]

    ##TODO: KTY_pullback_generators = zeros(parent(f),length(

    for gen in KTY_ideal_generators
      KTY_pullback_generators = [KTY_pullback_generators; gen .* Fstar_gens]
    end

  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear

end

"""
This one saves the ideals I_n  as we go along

"""
function quasiFSplitHeight_lift_mingens_idealsaves(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  R = parent(f)
  gR, _ = grade(R)

  f = gR(f)

  vars = gens(gR)

  fpminus1 = f^(p-1)
  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  Fstar_gens = Fstar_basis(p,f)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2

  KTY_ideal_generators = [fpminus1] # we don't actually use the initial value for anything but clarity
  KTY_pullback_generators = fpminus1 .* Fstar_gens

  idealsaves = []

  while n ≤ cutoff

    # Step 2. Remove things not in the kernel of u
    println("Testing for height $n...")

    #display(KTY_pullback_generators)

    for i in eachindex(KTY_pullback_generators)
      g = KTY_pullback_generators[i]

      if !in_kernel_poly_frob_generator(p,g)
        # perhaps not the most efficient, but it will work
        KTY_pullback_generators[i] = zero(g)

        #println("Removed generator $i")

      end
    end

    

    # println("$(length(KTY_pullback_generators)) generators of the pullback")
    # display(KTY_pullback_generators)
    # I = ideal(KTY_pullback_generators[:])
    # println(I)
    # KTY_pullback_mingens = minimal_generating_set(I)
    # # this is the wrong thing
    # println(KTY_pullback_mingens)
    # println("$(length(KTY_pullback_mingens)) minimal generators")

    comps = Fstar_components(KTY_pullback_generators) 

    KTY_pullback_mingens = []

    for (component, hs) in comps

        I = ideal(hs)
        mingens = minimal_generating_set(I)
        varterm = prod(vars .^ component)
        newgens = varterm .* mingens

        i_gens = length(gens(I))
        min_i_gens = length(mingens)
        difference = i_gens - min_i_gens 
        if 1 ≤ difference
            println("Saved $difference generators")
        end

        for newgen in newgens
            push!(KTY_pullback_mingens, newgen)
        end
    end

    # Step 3. Apply θFstar
    KTY_ideal_generators = θFstar.(KTY_pullback_mingens)

    push!(idealsaves,KTY_ideal_generators)

    if 2 < length(idealsaves)
        if idealsaves[end] == idealsaves[end-1]
            return cutoff + 2
        end
    end

    #TODO: Step 3.5. add fpminus1 and take minimal generating set

    # Step 4. Check whether the sequence of ideals terminatres here

    allzero = true

    for generator in KTY_ideal_generators

      if generator == zero(poly)
        continue
      else
        allzero = false
      end

      if !inPowerOfVariableIdeal(p,p,generator)
        # We are quasi-F split of height n! Yay!!
        # return n
        return (n, idealsaves)
      end

    end

    allzero && return cutoff + 2 # the chain terminated early, provable infinity

    # We are not n-quasi-F-split, so check if we are n+1-quasi-F-split

    n = n + 1

    cutoff < n && continue # pretty much a break statement

    # Step 1. Compute the generators of the pullback.

    #TODO: remove once step 3.5 is implemented
    KTY_ideal_generators = [KTY_ideal_generators; fpminus1]

    KTY_pullback_generators = typeof(f)[]

    ##TODO: KTY_pullback_generators = zeros(parent(f),length(

    for gen in KTY_ideal_generators
      KTY_pullback_generators = [KTY_pullback_generators; gen .* Fstar_gens]
    end

  end

  return (cutoff + 1, idealsaves)
  # return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear

end

function quasiFSplitHeight_lift_mingens_wics(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  R = parent(f)
  gR, _ = grade(R)

  f = gR(f)

  vars = gens(gR)

  fpminus1 = f^(p-1)
  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)
  Fstar_gens = Fstar_basis(p,f)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2

  KTY_ideal_generators = [fpminus1] # we don't actually use the initial value for anything but clarity
  KTY_pullback_generators = fpminus1 .* Fstar_gens


  while n ≤ cutoff

    # Step 2. Remove things not in the kernel of u
    println("Testing for height $n...")

    #display(KTY_pullback_generators)

    for i in eachindex(KTY_pullback_generators)
      g = KTY_pullback_generators[i]

      if !in_kernel_poly_frob_generator(p,g)
        # perhaps not the most efficient, but it will work
        KTY_pullback_generators[i] = zero(g)

        #println("Removed generator $i")

      end
    end

    

    # println("$(length(KTY_pullback_generators)) generators of the pullback")
    # display(KTY_pullback_generators)
    # I = ideal(KTY_pullback_generators[:])
    # println(I)
    # KTY_pullback_mingens = minimal_generating_set(I)
    # # this is the wrong thing
    # println(KTY_pullback_mingens)
    # println("$(length(KTY_pullback_mingens)) minimal generators")

    comps = Fstar_components(KTY_pullback_generators) 

    KTY_pullback_mingens = []

    for (component, hs) in comps

        I = ideal(hs)
        mingens = minimal_generating_set(I)
        varterm = prod(vars .^ component)
        newgens = varterm .* mingens

        i_gens = length(gens(I))
        min_i_gens = length(mingens)
        difference = i_gens - min_i_gens 
        if 1 ≤ difference
            println("Saved $difference generators")
        end

        for newgen in newgens
            push!(KTY_pullback_mingens, newgen)
        end
    end

    println("Starting to apply ΘFstar")
            
    # Step 3. Apply θFstar
    # KTY_ideal_generators = θFstar.(KTY_pullback_mingens)

    KTY_ideal_generators = similar(KTY_pullback_mingens)
    println(length(KTY_ideal_generators))
    for i in eachindex(KTY_ideal_generators)
        @time KTY_ideal_generators[i] = θFstar(KTY_pullback_mingens[i])
    end


        #TODO: Step 3.5. add fpminus1 and take minimal generating set

    # Step 4. Check whether the sequence of ideals terminatres here

    allzero = true

    for generator in KTY_ideal_generators

      if generator == zero(poly)
        continue
      else
        allzero = false
      end

      if !inPowerOfVariableIdeal(p,p,generator)
        # We are quasi-F split of height n! Yay!!
        return n
      end

    end

    allzero && return cutoff + 2 # the chain terminated early, provable infinity

    # We are not n-quasi-F-split, so check if we are n+1-quasi-F-split

    n = n + 1

    cutoff < n && continue # pretty much a break statement

    # Step 1. Compute the generators of the pullback.

    #TODO: remove once step 3.5 is implemented
    KTY_ideal_generators = [KTY_ideal_generators; fpminus1]

    KTY_pullback_generators = typeof(f)[]

    ##TODO: KTY_pullback_generators = zeros(parent(f),length(

    for gen in KTY_ideal_generators
      KTY_pullback_generators = [KTY_pullback_generators; gen .* Fstar_gens]
    end

  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear

end


# MARK - Quasi 1 f split (garbage)
#
#"""
#Calculates if the hypersuface defined by the
#polynomial poly is F-split
#
#note that p must be prime for this to have mathematical meaning
#"""
#function is1FSplit(p,poly)
#  #maybe TODO: check that p is prime
#
#  !inPowerOfVariableIdeal(p,p,poly^(2p-2))
#
#end#function
#
#
#"""
#Calculates the quasi-1-F-split height (I hope)
#in the case that deg(poly) = nvars(parent(poly))
#
#cutoff is inclusive, so it should be the highest possible height
#
#Uses the lift-based algorithm to calculate Δ₁
#
#This is a real guess here.
#"""
#function quasi1FSplitHeight_2CY_lift(p,poly,cutoff)
#  N = length(gens(parent(poly)))
#
#  #!isHomog(poly,ofdegree=N) && return -1 # type instability problem??
#
#  is1FSplit(p,poly) && return 1
#
#  f = poly
#
#  Δ₁fpminus1 = Δ₁l(p,f^(2p-2))
#  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)
#
#  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
#  # Honestly, just calling the ideals I_n could get confusing IMO
#
#  n = 2
#  # The newest generator in the KTY ideal I_2.
#  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
#  # be concatenating on new generator at each step until the chain terminates.
#  # See Theorem 5.8 in 2204.10076
#  KTYideal_n_new_gen = θFstar(f^(2p-2))
#
#  while n ≤ cutoff
#    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity
#
#    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
#      # We are quasi-F split of height n! Yay!!
#      return n
#    end
#
#    n = n + 1
#    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
#  end




#end
