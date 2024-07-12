#module QFSGeneralCase
#
#using Oscar
#using Memoize
#using Combinatorics
#
#include("FrobSplittingInfra.jl")
#using FrobSplittingInfra

#include("griffiths-dwork-construction/Utils.jl")

function quasiFSplitHeight(p,poly,cutoff)
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

  #FIXME FIXME this fails to do Example 7.8 correctly
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
