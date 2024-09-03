#module QFSCalabiYau
#
#using Oscar
#using Memoize
#using Combinatorics
#
#include("./FrobSplittingInfra.jl")
#using .FrobSplittingInfra
#
#include("../GPUPolynomials.jl/benchmarks/Benchmarks.jl")
#include("../GPUPolynomials.jl/src/Delta1.jl")
#include("griffiths-dwork-construction/Utils.jl")

using InteractiveUtils

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

"""
function quasiFSplitHeight_CY(p,poly,cutoff)
  println("Running CPU CY...")

  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly
  println("Computing delta1...")
  Δ₁fpminus1 = Δ₁(p,f^(p-1))
  println("delta1 finished!")
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076
  KTYideal_n_new_gen = θFstar(f^(p-1))
  
  println("Finding height...")
  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

"""
function quasiFSplitHeight_CY_lift(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  println("calculating Δ₁...")
  @time Δ₁fpminus1 = Δ₁l(p,f^(p-1))
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076
  println("trying height 2...")
  @time KTYideal_n_new_gen = θFstar(f^(p-1))

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    println("trying height $n...")
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    @time KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

Uses the matrix representaion of θFstar to compute the height.

"""
function quasiFSplitHeight_CY_lift_matrix(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
  start_vector = vector(fpminus1,m)

  M = matrix_of_lin_op(θFstar,m,parent(f))
  display(M)

  zzs = zeros(parent(start_vector[1]),m)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076

  @time KTYideal_n_new_gen = M * start_vector

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

    if KTYideal_n_new_gen[critical_ind] != 0
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    @time KTYideal_n_new_gen = M * KTYideal_n_new_gen
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

Uses the matrix representaion of θFstar to compute the height,
and uses the all-in-one-step method for getting this matrix,
rather than repeatedly evaluating..

Note that the method matrix_of_multiply_then_split is currently broken,
so this gives wrong results
"""
function quasiFSplitHeight_CY_lift_matrix_combined(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)

  M = matrix_of_multiply_then_split(p,coefs,degs,m)

  #@time M = matrix_of_lin_op(θFstar,m,parent(f))
  println("matrix finished:")
  display(M)

  nMonomials = length(start_vector)
  zzs = zeros(parent(start_vector[1]),nMonomials)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076

  println("trying height $n")
  @time KTYideal_n_new_gen = M * start_vector

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

    if KTYideal_n_new_gen[critical_ind] != 0
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    println("trying height $n")
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    @time KTYideal_n_new_gen = M * KTYideal_n_new_gen
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function


"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

Uses the matrix representaion of θFstar to compute the height,
and uses the all-in-one-step method for getting this matrix,
rather than repeatedly evaluating..

"""
function quasiFSplitHeight_CY_lift_sort(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  fpminus1 = f^(p-1)

  Δ₁fpminus1 = Δ₁l(p,fpminus1)
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  m = N*(p-1)
  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
  start_vector = lift_to_Int64(vector(fpminus1,m))
  (coefs,degs) = Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)

  M = matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,m)


  nMonomials = length(start_vector)
  zzs = zeros(parent(start_vector[1]),nMonomials)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076

  KTYideal_n_new_gen = (M * start_vector) .% p

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    #println("new generator: $KTYideal_n_new_gen")
    #println("zeros: $zzs")
    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

#    println(KTYideal_n_new_gen)


    if KTYideal_n_new_gen[critical_ind] != 0
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    @time KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function


"""
    isFSplit2(prime, poly)

Return tuple of whether poly is F-split or not and poly ^ (prime - 1)
This method exists to actually save f^(p-1) if f isn't F-split
"""
function isFSplit2(prime, poly)
    fpminus1 = poly ^ (prime - 1)

    return !inPowerOfVariableIdeal(prime, prime, fpminus1), fpminus1
end

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁, doing the
raising to powers on the gpu

Uses the matrix representaion of θFstar to compute the height,
and uses the all-in-one-step method for getting this matrix,
rather than repeatedly evaluating..

"""
function quasiFSplitHeight_CY_lift_sort_gpu(p,poly,cutoff,pregen=nothing)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isfsplit, fpminus1 = isFSplit2(p, poly)
  isfsplit && return 1

  fpminus1_homog = MMPSingularities.Polynomials.HomogeneousPolynomial(fpminus1)

  if pregen === nothing
    pregen = MMPSingularities.delta1.pregen_delta1(size(fpminus1_homog, 2),p)
  end
  MMPSingularities.Polynomials.sort_to_kronecker_order(fpminus1_homog, pregen.key1)
  
  Δ₁fpminus1 = MMPSingularities.delta1(fpminus1_homog,p; pregen = pregen)

  m = N*(p-1)
  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
  start_vector = lift_to_Int64(vector(fpminus1,m))
#   println("converting to gpu representation")
  #@time (coefs,degs) = Benchmarks.convert_to_gpu_representation(Δ₁fpminus1)
  coefs, degs = (Δ₁fpminus1.coeffs, Δ₁fpminus1.degrees)
#   println("Δ₁ has $(size(degs,1)) terms")

#   println("creating matrix...")
  #open("llvm_dump_1.ll","a") do io
  # println(io,@code_llvm matrix_of_multiply_then_split_sortmodp_kronecker(p,coefs,degs,m),a)
  #end
  #println()
  #error()
  M = matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,m)
  # @time M = matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,m)

  #@time M = matrix_of_lin_op(θFstar,m,parent(f))
#   println("matrix finished")
  #display(M)

  nMonomials = length(start_vector)
  zzs = zeros(parent(start_vector[1]),nMonomials)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076

  #println("matrix: $M")

#   println("critical index: $critical_ind")

#   println("trying height $n")
#   @time KTYideal_n_new_gen = (M * start_vector) .% p
  KTYideal_n_new_gen = (M * start_vector) .% p

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    #println("new generator: $KTYideal_n_new_gen")
    #println("zeros: $zzs")
    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

#    println(KTYideal_n_new_gen)

    # println("critical value: $(KTYideal_n_new_gen[critical_ind])")

    if KTYideal_n_new_gen[critical_ind] != 0
      # We are quasi-F split of height n! Yay!!
    #   println("height found!: $n")
      return n
    end

    n = n + 1
    # println("trying height $n")
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    # @time KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
    KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function


"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

This one uses multiply_then_split to only keep track of terms
that it needs.

"""
function quasiFSplitHeight_CY_lift_lazy(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  Δ₁fpminus1 = Δ₁l(p,f^(p-1))
  θFstar(a) = multiply_then_split(p,Δ₁fpminus1,a)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076
  KTYideal_n_new_gen = θFstar(f^(p-1))

  while n ≤ cutoff
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
  end

  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end#function


"""
Calculates the quasi-F-split height of the hypersurface
defined by poly.

Uses the naive formula in Theorem 5.8 (pg 42) of arXiv:2204.10076

This function seems to be broken right now, its results don't agree
with Table 2 in 2204.10076
FIXME: is this^^ still true?
"""
function quasiFSplitHeight_CY_naive_expansion(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  isFSplit(p,poly) && return 1

  f = poly

  n = 1
  fn = f^(p-1)

  while n < cutoff
    n = n + 1
    fn = fn * Δ₁(p,f^(p-1))^(p^(n-2))
    if !inPowerOfVariableIdeal(p,p^n,fn)
      # We are quasi-F-split of height n!! Yay!
      return n
    end
  end

  # We don't know what the quasi-F-split height is, it might
  # be infinity or it might just be greater than cutoff.
  return cutoff + 1
end#function

function quasiFSplitHeight_CY_gpu(p,poly,cutoff,pregen=nothing)
    # println("Doing check...")
    N = length(gens(parent(poly)))

    f = poly
    !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1

    fpminus1_gpu = Benchmarks.convert_to_gpu_representation(fpminus1)
    fpminus1_homog = HomogeneousPolynomial(fpminus1_gpu...)

    if pregen === nothing
        pregen = pregen_delta1(size(fpminus1_homog, 2),p)
    end
    sort_to_kronecker_order(fpminus1_homog, pregen.key1)
    
    Δ₁fpminus1_gpu = delta1(fpminus1_homog,p,pregen)

    R, (x, y, z, w) = polynomial_ring(poly.parent.base_ring, 4)
    Δ₁fpminus1 = zero(R)

    for (i, coeff) in enumerate(Δ₁fpminus1_gpu.coeffs)
        exp_row = Δ₁fpminus1_gpu.degrees[i, :]
        term = coeff * x^exp_row[1] * y^exp_row[2] * z^exp_row[3] * w^exp_row[4]
        Δ₁fpminus1 += term
    end

    θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

    n = 2
    KTYideal_n_new_gen = θFstar(f^(p-1))
  
    # println("Finding height...")
    while n ≤ cutoff
        #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
        KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
    # We are quasi-F split of height n! Yay!!
        return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
    end

    return cutoff + 1
end


#end
