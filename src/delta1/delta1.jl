import .GPUPolynomials.OperationPlan

include("int128stuff.jl")

nvars(x::CufpMPolyRingElem) = x.parent.nvars

struct Δ₁Plan <: OperationPlan
    numVars::Int
    prime::Int
    key::Int
    fftLen::Int
    totalDegree::Int
    primeArray::Vector
    nttType::DataType
    nttPowPlans::Vector{GPUPolynomials.NTTPowPlan}
    crtPlan::CuMatrix
    memoryefficient::Bool
end

function plan_Δ₁(numVars, prime)::Δ₁Plan
    memoryefficient = false
    if (numVars, prime) == (4, 2)
        primeArray = UInt.([12289])
    elseif (numVars, prime) == (4, 3)
        # primeArray = UInt64.([114689])
        primeArray = UInt.([0x3ffffff960000001])
    elseif (numVars, prime) == (4, 5)
        # primeArray = UInt32.([13631489, 23068673])
        primeArray = UInt.([0x3ffffff960000001])
    elseif (numVars, prime) == (4, 7)
        primeArray = UInt.([0x3ffffff960000001, 0x3ffffff760000001])
    elseif (numVars, prime) == (4, 11)
        primeArray = UInt.([0x3ffffff960000001, 0x3ffffff760000001, 0x3fffffeec0000001,  0x3fffffee60000001])
        memoryefficient = true
    elseif (numVars, prime) == (4, 13)
        primeArray = UInt.([0x3ffffff960000001, 0x3ffffff760000001, 0x3fffffeec0000001,  0x3fffffee60000001, 0x3fffffee00000001])
        memoryefficient = true
    else
        throw(ArgumentError("I haven't figured out bounds for this yet!"))
    end

    resultTotalDegree = numVars * (prime - 1) * prime
    key = resultTotalDegree + 1
    fftLen = Base._nextpow2(resultTotalDegree * key^(numVars - 2) + 1)
    
    nttPowPlans = GPUPolynomials.NTTPowPlan[]
    # @assert all(isprime.(primeArray)) # yeah idk
    for p in primeArray
        nttPowPlan = GPUPolynomials.NTTPowPlan(fftLen, prime, p; memoryefficient = memoryefficient)
        push!(nttPowPlans, nttPowPlan)
    end
    resultDataType = GPUPolynomials.get_uint_type(max(Base._nextpow2(Int(ceil(log2(prod(BigInt.(primeArray)))))), 64))
    crtPlan = GPUPolynomials.plan_crt(resultDataType.(primeArray))

    return Δ₁Plan(numVars, prime, key, fftLen, resultTotalDegree, primeArray, eltype(primeArray), nttPowPlans, crtPlan, memoryefficient)
end

function Δ₁(g::CufpMPolyRingElem)
    if !(g.opPlan isa Δ₁Plan)
        throw(ArgumentError("Input polynomial needs an OperationPlan!"))
    end

    if g.opPlan.memoryefficient
        memoryefficient_Δ₁(g)
    else
        fast_Δ₁(g)
    end
end

function fast_Δ₁(g::CufpMPolyRingElem)
    numVars = nvars(g)

    vecs = GPUPolynomials.get_dense_representation(g, g.opPlan.fftLen, g.bits, g.opPlan.nttType, g.opPlan.key, length(g.opPlan.nttPowPlans))

    currPtr = pointer(vecs)
    for planNum in eachindex(g.opPlan.nttPowPlans)
        vect = CUDA.unsafe_wrap(CuVector{g.opPlan.nttType}, currPtr, g.opPlan.fftLen)
        GPUPolynomials.ntt_pow(vect, g.opPlan.nttPowPlans[planNum])
        currPtr += sizeof(g.opPlan.nttType) * g.opPlan.fftLen
    end

    remove_pth_power_terms(g, g.opPlan.key, vecs, g.opPlan.prime, g.opPlan.primeArray)

    multimodResultCoeffs, encodedDegs = GPUPolynomials.sparsify(vecs)

    resultCoeffs = GPUPolynomials.build_result(multimodResultCoeffs, g.opPlan.crtPlan)
    # @assert all(x -> x % eltype(resultCoeffs)(g.opPlan.prime) == zero(eltype(resultCoeffs)), Array(resultCoeffs))
    divide_and_mod!(resultCoeffs, g.opPlan.prime)
    resultCoeffs = UInt64.(resultCoeffs)

    resultDegs = GPUPolynomials.kronecker_to_bitpacked(encodedDegs, g.opPlan.key, numVars, g.opPlan.totalDegree, g.bits, UInt)

    return CufpMPolyRingElem(resultCoeffs, resultDegs, g.bits, true, g.opPlan.totalDegree, g.parent, GPUPolynomials.EmptyPlan())
end


function memoryefficient_Δ₁(g::CufpMPolyRingElem)
    numVars = nvars(g)

    vecs = GPUPolynomials.cpu_get_dense_representation(g, g.opPlan.fftLen, g.bits, g.opPlan.nttType, g.opPlan.key, length(g.opPlan.nttPowPlans))

    currPtr = pointer(vecs)
    gpualloc = CUDA.zeros(g.opPlan.nttType, g.opPlan.fftLen)
    for planNum in eachindex(g.opPlan.nttPowPlans)
        cpuvec = unsafe_wrap(Vector{g.opPlan.nttType}, currPtr, g.opPlan.fftLen)
        copyto!(gpualloc, cpuvec)
        GPUPolynomials.ntt_pow(gpualloc, g.opPlan.nttPowPlans[planNum])
        copyto!(cpuvec, gpualloc)
        currPtr += sizeof(g.opPlan.nttType) * g.opPlan.fftLen
    end
    
    remove_pth_power_terms(g, g.opPlan.key, vecs, g.opPlan.prime, g.opPlan.primeArray)

    multimodResultCoeffs, encodedDegs = GPUPolynomials.sparsify(Array(vecs))
    encodedDegs = CuArray(encodedDegs)

    # resultCoeffs = GPUPolynomials.cpu_build_result(multimodResultCoeffs, Array(g.opPlan.crtPlan))
    crtPlan = Array(g.opPlan.crtPlan)
    resultCoeffs = zeros(eltype(crtPlan), size(multimodResultCoeffs, 1))

    for i in axes(multimodResultCoeffs, 1)
        subarr = view(multimodResultCoeffs, i, :)
        x = eltype(crtPlan)(subarr[1])
        for j in axes(crtPlan, 2)
            a = mul_mod(x, crtPlan[2, j], crtPlan[3, j])
            b = mul_mod(eltype(crtPlan)(subarr[j + 1]), crtPlan[1, j], crtPlan[3, j])
            x = add_mod(a, b, crtPlan[3, j])
        end

        resultCoeffs[i] = x
    end

    @assert all(x -> x % eltype(resultCoeffs)(g.opPlan.prime) == zero(eltype(resultCoeffs)), resultCoeffs)

    p = eltype(resultCoeffs)(g.opPlan.prime)
    cpu_resultCoeffs = Array(resultCoeffs)
    @assert all(x -> x % p == 0, cpu_resultCoeffs)
    cpu_resultCoeffs .÷= p
    cpu_resultCoeffs .%= p
    cpu_resultCoeffs = UInt64.(cpu_resultCoeffs)
    resultCoeffs = CuArray(cpu_resultCoeffs)

    resultDegs = GPUPolynomials.kronecker_to_bitpacked(encodedDegs, g.opPlan.key, numVars, g.opPlan.totalDegree, g.bits, UInt)

    return CufpMPolyRingElem(resultCoeffs, resultDegs, g.bits, true, g.opPlan.totalDegree, g.parent, GPUPolynomials.EmptyPlan())
end

function divide_and_mod!(coeffs::CuVector{T}, prime::Integer) where T<:Unsigned
    p = T(prime)

    kernel = @cuda launch=false divide_and_mod_kernel!(coeffs, p)
    config = launch_configuration(kernel.fun)
    threads = min(length(coeffs), config.threads)
    blocks = cld(length(coeffs), threads)

    CUDA.@sync kernel(coeffs, p; threads = threads, blocks = blocks)
end

function divide_and_mod_kernel!(coeffs::CuDeviceVector{T}, prime::T) where T<:Unsigned
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    
    if idx <= length(coeffs)
        @inbounds begin
            coeffs[idx] = unchecked_div(coeffs[idx], prime)
            coeffs[idx] = unchecked_mod(coeffs[idx], prime)
        end
    end

    return nothing
end

function generate_remove_indices(intermediate::CufpMPolyRingElem, key::Int, pow::Int)
    coeffs = Array(intermediate.coeffs)
    degrees = Array(intermediate.exps)
    bits = intermediate.bits
    mask = (one(eltype(coeffs)) << bits) - 1

    keyPowers = [key ^ i for i in 0:nvars(intermediate) - 2]

    indices = zeros(Int, length(coeffs))
    subtract = zeros(eltype(coeffs), length(coeffs))
    degrees .*= pow
    for i in eachindex(degrees)
        resultidx = 1
        deg = degrees[i]
        for i in 1:nvars(intermediate) - 1
            resultidx += (deg & mask) * keyPowers[i]
            deg >>= bits
        end
        indices[i] = resultidx
        subtract[i] = coeffs[i] ^ pow
    end

    return indices, subtract
end

function remove_pth_power_terms(intermediate::CufpMPolyRingElem, key::Int, vec, p::Int, primeArray)
    removeindices, subtract = generate_remove_indices(intermediate, key, p)
    if vec isa CuArray
        gpu_remove_pth_power_terms(CuArray(removeindices), CuArray(subtract), vec, CuArray(primeArray))
    else
        cpu_remove_pth_power_terms(removeindices, subtract, vec, primeArray)
    end

    return nothing
end

function cpu_remove_pth_power_terms(removeindices, subtract, arr, primeArray)
    for i in eachindex(removeindices)
        removeidx = removeindices[i]
        for p in axes(arr, 2)
            arr[removeidx, p] = sub_mod(arr[removeidx, p], subtract[i], primeArray[p])
        end
    end

    return nothing
end

function gpu_remove_pth_power_terms(removeindices::CuArray, subtract::CuArray, arr::CuArray, primeArray::CuArray)
    kernel = @cuda launch=false gpu_remove_pth_power_terms_kernel!(removeindices, subtract, arr, primeArray)
    config = launch_configuration(kernel.fun)
    threads = min(length(removeindices), config.threads)
    blocks = cld(length(removeindices), threads)

    kernel(removeindices, subtract, arr, primeArray; threads = threads, blocks = blocks)
end

function gpu_remove_pth_power_terms_kernel!(removeindices, subtract, vec, primeArray)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if idx <= length(removeindices)
        removeidx = removeindices[idx]
        for p in axes(vec, 2)
            vec[removeidx, p] = sub_mod(vec[removeidx, p], subtract[idx], primeArray[p])
        end
    end

    return nothing
end