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

  m = N*(p-1)
  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
  start_vector = lift_to_Int64(vector(fpminus1,m))
  hp = HomogeneousPolynomial(Δ₁fpminus1)

  coefs = hp.coeffs
  degs = hp.degrees

  M = matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,m)

  nMonomials = length(start_vector)
  zzs = zeros(parent(start_vector[1]),nMonomials)

  n = 2

  KTYideal_n_new_gen = (M * start_vector) .% p

  while n ≤ cutoff

    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

    if KTYideal_n_new_gen[critical_ind] != 0
      # We are quasi-F split of height n! Yay!!
      return n
    end

    n = n + 1
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
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

function quasiFSplitHeight_CY_lift_sort_gpu_second_step(p, fpminus1, cutoff, pregen)
    N = length(gens(parent(fpminus1)))
    
    fpminus1_homog = MMPSingularities.Polynomials.HomogeneousPolynomial(fpminus1)

    Δ₁fpminus1 = MMPSingularities.GPUDelta1.delta1(fpminus1_homog,p; pregen = pregen)
  
    m = N*(p-1)
    critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
    start_vector = lift_to_Int64(vector(fpminus1,m))

    coefs, degs = (Δ₁fpminus1.coeffs, Δ₁fpminus1.degrees)


    M = matrix_of_multiply_then_split_sortmodp_kronecker2(p,coefs,degs,m)
  
    nMonomials = length(start_vector)
    zzs = zeros(parent(start_vector[1]),nMonomials)
  
    n = 2

    KTYideal_n_new_gen = (M * start_vector) .% p
  
    while n ≤ cutoff
        KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity
        if KTYideal_n_new_gen[critical_ind] != 0
            return n
        end
    
        n = n + 1
        KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
    end
  
    return cutoff + 1 #
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
function quasiFSplitHeight_CY_lift_sort_gpu(p,poly,cutoff,pregen)
    N = length(gens(parent(poly)))
  
    !isHomog(poly,ofdegree=N) && return -1 # type instability problem??
  
    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1
  
    fpminus1_homog = MMPSingularities.Polynomials.HomogeneousPolynomial(fpminus1)

    Δ₁fpminus1 = MMPSingularities.GPUDelta1.delta1(fpminus1_homog,p; pregen = pregen)
  
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