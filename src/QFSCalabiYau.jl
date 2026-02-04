
"""
In the case that f is Calabi-Yau
(degree = number of variables)
Calculates the quasi-F-split height
ht(f) of
f if ht(f) ≤ cutoff. 
Otherwise, 
return a value that is bigger than b.

"""
function quasiFSplitHeight_CY(p,poly,cutoff,pregen=nothing)
    if p ≤ 3
        quasiFSplitHeight_CY_lift(p,poly,cutoff)
    else
        if pregen == nothing
            n = length(gens(parent(poly)))
            
            println("No pregen found. Creating one...")
            @time pregen = pregen_qfsheight(n, p)
        end
        h = quasiFSplitHeight_CY_lift_wics_gpu(p,poly,cutoff,pregen) 
        h
    end
end


"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁

"""
function quasiFSplitHeight_CY_lift(p,poly,cutoff)
  N = length(gens(parent(poly)))

  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

  push!(fpminus1_time, (@timed begin 
  isfsplit, fpminus1 = isFSplit2(p, poly)
  isfsplit && return 1
  end).time)
  
  push!(delta1_time, (@timed begin
  Δ₁fpminus1 = Δ₁lp²(fpminus1)
  end).time)
  
  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
  # Honestly, just calling the ideals I_n could get confusing IMO

  n = 2
  # The newest generator in the KTY ideal I_2.
  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
  # be concatenating on new generator at each step until the chain terminates.
  # See Theorem 5.8 in 2204.10076
  # println("trying height 2...")

  push!(stripe_mul_time, (@timed begin
    KTYideal_n_new_gen = θFstar(fpminus1)
  end).time)
  
  while n ≤ cutoff
    push!(if_time, (@timed begin
    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
      # We are quasi-F split of height n! Yay!!
      return n
    end
    end).time)
    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
    
    n = n + 1

    push!(stripe_mul_time, (@timed begin
      KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
    end).time)
    # println("trying height $n...")
    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
    
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

struct QFSHeightPregen
    Δ₁plan::Δ₁lp²Plan
    momtspregen::MOMTSPregen
end

function pregen_qfsheight(n, p)
    Δ₁plan = plan_Δ₁lp²(n, p)
    momtspregen = pregen_MOMTS(n, p)

    return QFSHeightPregen(Δ₁plan, momtspregen)
end

"""
Calculates the quasi-F-split height
in the case that deg(poly) = nvars(parent(poly))

cutoff is inclusive, so it should be the highest possible height

Uses the lift-based algorithm to calculate Δ₁, doing the
raising to powers on the gpu

Uses the matrix representaion of θFstar to compute the height,
and uses the all-in-one-step method for getting this matrix,
rather than repeatedly evaluating.
"""
function _pad_copy_matrix_float32(M::AbstractMatrix)
    m, n = size(M)
    dest = CUDA.zeros(Float32, m + 32, n + 32)
    src = CuArray(Float32.(M))
    copyto!(view(dest, 1:m, 1:n), src)
    return dest
end

function _pad_copy_vector_float32(v::AbstractVector)
    n = length(v)
    dest = CUDA.zeros(Float32, n + 32, 1 + 32)
    src = reshape(CuArray(Float32.(v)), n, 1)
    copyto!(view(dest, 1:n, 1), src)
    return dest
end

function quasiFSplitHeight_CY_lift_wics_gpu(p,poly,cutoff,pregen)
    N = length(gens(parent(poly)))
  
    !isHomog(poly,ofdegree=N) && return -1
  
    push!(fpminus1_time, (@timed begin 
    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1
    end).time)
  
    push!(delta1_time, (@timed begin
    fpminus1_gpu = CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = pregen.Δ₁plan
    Δ₁fpminus1 = Δ₁lp²(fpminus1_gpu)
    end).time)
  
    m = N*(p-1)
    critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
    start_vector = lift_to_Int64(vector(fpminus1,m))
  
    push!(move_matrix_time, (@timed begin
    M_data = matrix_of_multiply_then_split(Δ₁fpminus1; plan = pregen.momtspregen, alg = 6)
    m_rows, m_cols = size(M_data)
    m_rows -= 32
    m_cols -= 32
    M = CuModMatrix(M_data, p; new_size=(m_rows, m_cols))
    end).time)

    push!(move_vector_time, (@timed begin
    nMonomials = length(start_vector)
    
    start_vector_data = _pad_copy_vector_float32(start_vector)
    start_vector_gpu = CuModMatrix(start_vector_data, p; new_size=(nMonomials, 1))

    padded_out_rows = m_rows + GPUFiniteFieldMatrices.TILE_WIDTH
    padded_out_cols = 1 + GPUFiniteFieldMatrices.TILE_WIDTH

    KTYideal = CuModMatrix(CUDA.zeros(Float32, padded_out_rows, padded_out_cols), p; new_size=(m_rows, 1))
    KTYideal_next = CuModMatrix(CUDA.zeros(Float32, padded_out_rows, padded_out_cols), p; new_size=(m_rows, 1))
    end).time)
  
    n = 2
  
    push!(stripe_mul_time, (@timed begin
    GPUFiniteFieldMatrices.stripe_mul!(KTYideal, M, start_vector_gpu)
    end).time)
  
    while n ≤ cutoff
        push!(if_time, (@timed begin
      if CUDA.all(iszero, @view KTYideal.data[1:m_rows, 1])
        return cutoff + 2
      end
  
      if CUDA.@allowscalar KTYideal.data[critical_ind, 1] != 0
        return n
      end
        end).time)
  
      n = n + 1
  
      push!(stripe_mul_time, (@timed begin
      GPUFiniteFieldMatrices.stripe_mul!(KTYideal_next, M, KTYideal)
      end).time)
      KTYideal, KTYideal_next = KTYideal_next, KTYideal
    end
    return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end

function quasiFSplitHeight_CY_lift_wics_gpu_cpu_check(p, poly, cutoff, pregen)
    N = length(gens(parent(poly)))

    !isHomog(poly, ofdegree = N) && return -1

    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1

    fpminus1_gpu = CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = pregen.Δ₁plan
    Δ₁fpminus1 = Δ₁lp²(fpminus1_gpu)

    m = N * (p - 1)
    critical_ind = index_of_term_not_in_frobenius_power_CY(p, N)
    start_vector_cpu = lift_to_Int64(vector(fpminus1, m))

    M_cpu = Array(matrix_of_multiply_then_split(Δ₁fpminus1; plan = pregen.momtspregen, alg = 5))
    nMonomials = length(start_vector_cpu)
    zzs = zeros(parent(start_vector_cpu[1]), nMonomials)
    KTY_cpu = (M_cpu * start_vector_cpu) .% p

    M_original = matrix_of_multiply_then_split(Δ₁fpminus1; plan = pregen.momtspregen, alg = 5)
    m_rows, m_cols = size(M_original)
    M_data = _pad_copy_matrix_float32(M_original)
    M = CuModMatrix(M_data, p; new_size = (m_rows, m_cols))
    mod_elements!(M)

    start_vector_data = _pad_copy_vector_float32(start_vector_cpu)
    start_vector_gpu = CuModMatrix(start_vector_data, p; new_size = (nMonomials, 1))
    mod_elements!(start_vector_gpu)

    padded_out_rows = m_rows + GPUFiniteFieldMatrices.TILE_WIDTH
    padded_out_cols = 1 + GPUFiniteFieldMatrices.TILE_WIDTH
    KTY_gpu = CuModMatrix(CUDA.zeros(Float32, padded_out_rows, padded_out_cols), p; new_size = (m_rows, 1))
    KTY_gpu_next = CuModMatrix(CUDA.zeros(Float32, padded_out_rows, padded_out_cols), p; new_size = (m_rows, 1))

    M_gpu_cpu = Array(@view M.data[1:m_rows, 1:m_cols])
    if M_cpu != M_gpu_cpu
        mismatch_idx = findfirst(i -> M_cpu[i] != M_gpu_cpu[i], 1:length(M_cpu))
        cpu_val = M_cpu[mismatch_idx]
        gpu_val = M_gpu_cpu[mismatch_idx]
        error("GPU/CPU mismatch at M index=$mismatch_idx cpu=$cpu_val gpu=$gpu_val")
    end

    start_vec_gpu_cpu = vec(Array(@view start_vector_gpu.data[1:nMonomials, 1]))
    if start_vector_cpu != start_vec_gpu_cpu
        mismatch_idx = findfirst(i -> start_vector_cpu[i] != start_vec_gpu_cpu[i], 1:length(start_vector_cpu))
        cpu_val = start_vector_cpu[mismatch_idx]
        gpu_val = start_vec_gpu_cpu[mismatch_idx]
        error("GPU/CPU mismatch at start_vector index=$mismatch_idx cpu=$cpu_val gpu=$gpu_val")
    end

    GPUFiniteFieldMatrices.stripe_mul!(KTY_gpu, M, start_vector_gpu)

    KTY_gpu_cpu = vec(Array(@view KTY_gpu.data[1:m_rows, 1]))
    if KTY_cpu != KTY_gpu_cpu
        mismatch_idx = findfirst(i -> KTY_cpu[i] != KTY_gpu_cpu[i], 1:length(KTY_cpu))
        cpu_val = KTY_cpu[mismatch_idx]
        gpu_val = KTY_gpu_cpu[mismatch_idx]
        println("KTY_cpu: ", Int.(KTY_cpu))
        println("KTY_gpu_cpu: ", KTY_gpu_cpu)
        error("GPU/CPU mismatch at n=1 index=$mismatch_idx cpu=$cpu_val gpu=$gpu_val")
    end

    n = 2
    while n ≤ cutoff

        KTY_gpu_cpu = vec(Array(@view KTY_gpu.data[1:m_rows, 1]))
        if KTY_cpu != KTY_gpu_cpu
            mismatch_idx = findfirst(i -> KTY_cpu[i] != KTY_gpu_cpu[i], 1:length(KTY_cpu))
            cpu_val = KTY_cpu[mismatch_idx]
            gpu_val = KTY_gpu_cpu[mismatch_idx]
            println("KTY_cpu: ", Int.(KTY_cpu))
            println("KTY_gpu_cpu: ", KTY_gpu_cpu)
            error("GPU/CPU mismatch at n=$n index=$mismatch_idx cpu=$cpu_val gpu=$gpu_val")
            
        end

        KTY_cpu == zzs && return cutoff + 2

        if KTY_cpu[critical_ind] != 0
            return n
        end

        n = n + 1
        
        KTY_cpu = (M_cpu * KTY_cpu) .% p
        GPUFiniteFieldMatrices.stripe_mul!(KTY_gpu_next, M, KTY_gpu)
        KTY_gpu, KTY_gpu_next = KTY_gpu_next, KTY_gpu
    end
    return cutoff + 1
end

function quasiFSplitHeight_CY_lift_wics_cpu(p,poly,cutoff,pregen)
    N = length(gens(parent(poly)))
  
    !isHomog(poly,ofdegree=N) && return -1
  
    push!(fpminus1_time, (@timed begin 
    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1
    end).time)
  
    push!(delta1_time, (@timed begin
    Δ₁fpminus1 = Δ₁lp²(fpminus1)
    end).time)
  
    m = N*(p-1)
    critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
    start_vector = lift_to_Int64(vector(fpminus1,m))
  
    push!(move_matrix_time, (@timed begin
    M = matrix_of_multiply_then_split(Δ₁fpminus1; plan = pregen.momtspregen, alg = 4)
    end).time)
    nMonomials = length(start_vector)
    zzs = zeros(parent(start_vector[1]),nMonomials)
  
    n = 2
    
    push!(stripe_mul_time, (@timed begin
    KTYideal_n_new_gen = (M * start_vector) .% p
    end).time)
  
    while n ≤ cutoff
      push!(if_time, (@timed begin
      KTYideal_n_new_gen == zzs && return cutoff + 2
  
      if KTYideal_n_new_gen[critical_ind] != 0
        return n
      end
      end).time)
  
      n = n + 1
      
      push!(stripe_mul_time, (@timed begin
      KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
      end).time)
    end
    return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end

function quasiFSplitHeight_CY_lift_wics_cpu_momts(p,poly,cutoff,pregen)
    N = length(gens(parent(poly)))
  
    !isHomog(poly,ofdegree=N) && return -1
  
    isfsplit, fpminus1 = isFSplit2(p, poly)
    isfsplit && return 1
  
    fpminus1_gpu = CufpMPolyRingElem(fpminus1.data, UInt64)
    fpminus1_gpu.opPlan = pregen.Δ₁plan
    Δ₁fpminus1 = Δ₁lp²(fpminus1_gpu)
  
    m = N*(p-1)
    critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
    start_vector = lift_to_Int64(vector(fpminus1,m))
  
  
    M = Array(matrix_of_multiply_then_split(Δ₁fpminus1; plan = pregen.momtspregen, alg = 5))
    nMonomials = length(start_vector)
    zzs = zeros(parent(start_vector[1]),nMonomials)
  
    n = 2
  
    KTYideal_n_new_gen = (M * start_vector) .% p
  
    while n ≤ cutoff
      KTYideal_n_new_gen == zzs && return cutoff + 2
  
      if KTYideal_n_new_gen[critical_ind] != 0
        return n
      end
  
      n = n + 1
  
      KTYideal_n_new_gen = (M * KTYideal_n_new_gen) .% p
    end
    return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
end

# """
# Uses the gpu to calculate Δ_1 and then 
# finds the quasi-F-split height using the classical
# polynomial multiplication algorithm (in OSCAR)
# """
# function quasiFSplitHeight_CY_gpu(p,poly,cutoff,pregen=nothing)
#     # println("Doing check...")
#     N = length(gens(parent(poly)))

#     f = poly
#     !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

#     isfsplit, fpminus1 = isFSplit2(p, poly)
#     isfsplit && return 1

#     fpminus1_gpu = convert_to_gpu_representation(fpminus1)
#     fpminus1_homog = GPUPolynomials.HomogeneousPolynomial(fpminus1_gpu...)

#     if pregen === nothing
#         pregen = pregen_delta1(size(fpminus1_homog, 2),p)
#     end
#     GPUPolynomials.sort_to_kronecker_order(fpminus1_homog, pregen.key1)
    
#     Δ₁fpminus1_gpu = delta1(fpminus1_homog,p;pregen)

#     R, (x, y, z, w) = polynomial_ring(poly.parent.base_ring, 4)
#     Δ₁fpminus1 = zero(R)

#     for (i, coeff) in enumerate(Δ₁fpminus1_gpu.coeffs)
#         exp_row = Δ₁fpminus1_gpu.degrees[i, :]
#         term = coeff * x^exp_row[1] * y^exp_row[2] * z^exp_row[3] * w^exp_row[4]
#         Δ₁fpminus1 += term
#     end

#     θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

#     n = 2
#     KTYideal_n_new_gen = θFstar(f^(p-1))
  
#     # println("Finding height...")
#     while n ≤ cutoff
#         #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#         KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

#         if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
#         # We are quasi-F split of height n! Yay!!
#             return n
#         end

#         n = n + 1
#         #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#         KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
#     end

#     return cutoff + 1
# end

# MARK - other methods
#
# The following are also algorithms that can be used to calculate
# the quasi-F-split height. 
# They aren't as efficient as the methods above, but they can 
# be used to double check results.

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

#"""
#Calculates the quasi-F-split height
#in the case that deg(poly) = nvars(parent(poly))

#cutoff is inclusive, so it should be the highest possible height

#Uses the lift-based algorithm to calculate Δ₁

#This one uses multiply_then_split to only keep track of terms
#that it needs.

#"""
#function quasiFSplitHeight_CY_lift_lazy(p,poly,cutoff)
#  N = length(gens(parent(poly)))

#  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

#  isFSplit(p,poly) && return 1

#  f = poly

#  Δ₁fpminus1 = Δ₁l(p,f^(p-1))
#  θFstar(a) = multiply_then_split(p,Δ₁fpminus1,a)

#  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
#  # Honestly, just calling the ideals I_n could get confusing IMO

#  n = 2
#  # The newest generator in the KTY ideal I_2.
#  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
#  # be concatenating on new generator at each step until the chain terminates.
#  # See Theorem 5.8 in 2204.10076
#  KTYideal_n_new_gen = θFstar(f^(p-1))

#  while n ≤ cutoff
#    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

#    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
#      # We are quasi-F split of height n! Yay!!
#      return n
#    end

#    n = n + 1
#    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
#  end

#  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
#end#function

# """
# Calculates the quasi-F-split height
# in the case that deg(poly) = nvars(parent(poly))

# cutoff is inclusive, so it should be the highest possible height

# Uses the lift-based algorithm to calculate Δ₁

# Uses the matrix representaion of θFstar to compute the height,
# and uses the all-in-one-step method for getting this matrix,
# rather than repeatedly evaluating..

# This uses the method `matrix_of_multiply_then_split` to get that
# matrix.

# Note that the method matrix_of_multiply_then_split is currently broken,
# so this gives wrong results
# """
# function quasiFSplitHeight_CY_lift_matrix_combined(p,poly,cutoff)
#   N = length(gens(parent(poly)))

#   !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

#   isFSplit(p,poly) && return 1

#   f = poly

#   fpminus1 = f^(p-1)

#   Δ₁fpminus1 = Δ₁l(p,fpminus1)
#   θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

#   m = N*(p-1)
#   critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
#   start_vector = lift_to_Int64(vector(fpminus1,m))
#   @time (coefs,degs) = convert_to_gpu_representation(Δ₁fpminus1)
#   println("Δ₁ has $(size(degs,1)) terms")

#   println("creating matrix...")
#   @time M = matrix_of_multiply_then_split(p,coefs,degs,m)

#   #@time M = matrix_of_lin_op(θFstar,m,parent(f))
#   println("matrix finished:")
#   display(M)

#   nMonomials = length(start_vector)
#   zzs = zeros(parent(start_vector[1]),nMonomials)

#   # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
#   # Honestly, just calling the ideals I_n could get confusing IMO

#   n = 2
#   # The newest generator in the KTY ideal I_2.
#   # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
#   # be concatenating on new generator at each step until the chain terminates.
#   # See Theorem 5.8 in 2204.10076

#   println("trying height $n")
#   @time KTYideal_n_new_gen = M * start_vector

#   while n ≤ cutoff
#     #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#     KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

#     if KTYideal_n_new_gen[critical_ind] != 0
#       # We are quasi-F split of height n! Yay!!
#       return n
#     end

#     n = n + 1
#     println("trying height $n")
#     #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#     @time KTYideal_n_new_gen = M * KTYideal_n_new_gen
#   end

#   return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
# end#function


#"""
#Calculates the quasi-F-split height
#in the case that deg(poly) = nvars(parent(poly))

#cutoff is inclusive, so it should be the highest possible height

#Uses the lift-based algorithm to calculate Δ₁

#Uses the matrix representaion of θFstar to compute the height.

#this uses the method `matrix_of_lin_op` to calculate
#the matrix of first multiplying and then applying the splitting.
#"""
#function quasiFSplitHeight_CY_lift_matrix(p,poly,cutoff)
#  N = length(gens(parent(poly)))

#  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

#  isFSplit(p,poly) && return 1

#  f = poly

#  fpminus1 = f^(p-1)

#  Δ₁fpminus1 = Δ₁l(p,fpminus1)
#  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

#  m = N*(p-1)
#  critical_ind = index_of_term_not_in_frobenius_power_CY(p,N) # lex order (i.e. the default)
#  start_vector = vector(fpminus1,m)
#  println("creating matrix...")
#  @time M = matrix_of_lin_op(θFstar,m,parent(f))
#  println("matrix finished:")
#  display(M)

#  zzs = zeros(parent(start_vector[1]),m)

#  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
#  # Honestly, just calling the ideals I_n could get confusing IMO

#  n = 2
#  # The newest generator in the KTY ideal I_2.
#  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
#  # be concatenating on new generator at each step until the chain terminates.
#  # See Theorem 5.8 in 2204.10076

#  println("trying height $n")
#  @time KTYideal_n_new_gen = M * start_vector

#  while n ≤ cutoff
#    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#    KTYideal_n_new_gen == zzs && return cutoff + 2 # the chain terminated early, provable infinity

#    if KTYideal_n_new_gen[critical_ind] != 0
#      # We are quasi-F split of height n! Yay!!
#      return n
#    end

#    n = n + 1
#    println("trying height $n")
#    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#    @time KTYideal_n_new_gen = M * KTYideal_n_new_gen
#  end

#  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
#end#function

#"""
#Calculates the quasi-F-split height
#in the case that deg(poly) = nvars(parent(poly))

#cutoff is inclusive, so it should be the highest possible height

#This uses calculates Δ_1 by using a formula for the coefficients,
#i.e. by using multinomial coefficients.
#"""
#function quasiFSplitHeight_CY_formula(p,poly,cutoff)
#  N = length(gens(parent(poly)))

#  !isHomog(poly,ofdegree=N) && return -1 # type instability problem??

#  isFSplit(p,poly) && return 1

#  f = poly

#  Δ₁fpminus1 = Δ₁(p,f^(p-1))
#  θFstar(a) = polynomial_frobenius_generator(p,Δ₁fpminus1*a)

#  # KTY is for Kawakami, Takamatsu, and Yoshikawa, the authors of 2204.10076
#  # Honestly, just calling the ideals I_n could get confusing IMO

#  n = 2
#  # The newest generator in the KTY ideal I_2.
#  # For Calabi-Yau varieties, one has that the sequence I_n can be seen to
#  # be concatenating on new generator at each step until the chain terminates.
#  # See Theorem 5.8 in 2204.10076
#  KTYideal_n_new_gen = θFstar(f^(p-1))

#  while n ≤ cutoff
#    #println("New Generator of KTY ideal I_n: ", KTYideal_n_new_gen)
#    KTYideal_n_new_gen == zero(poly) && return cutoff + 2 # the chain terminated early, provable infinity

#    if !inPowerOfVariableIdeal(p,p,KTYideal_n_new_gen)
#      # We are quasi-F split of height n! Yay!!
#      return n
#    end

#    n = n + 1
#    #println("next one should be: ", θFstar(KTYideal_n_new_gen))
#    KTYideal_n_new_gen = θFstar(KTYideal_n_new_gen)
#  end

#  return cutoff + 1 # we didn't see the chain terminate, conclusion is unclear
#end#function
