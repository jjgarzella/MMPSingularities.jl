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
    memorySafe::Bool
end

function plan_Δ₁(numVars, prime)::Δ₁Plan
    memorySafe = false
    if (numVars, prime) == (4, 2)
        primeArray = UInt32.([12289])
    elseif (numVars, prime) == (4, 3)
        primeArray = UInt32.([114689])
    elseif (numVars, prime) == (4, 5)
        primeArray = UInt32.([13631489, 23068673])
    elseif (numVars, prime) == (4, 7)
        primeArray = UInt32.([167772161, 377487361, 469762049])
    elseif (numVars, prime) == (4, 11)
        primeArray = UInt.([2033101286932481, 2033107326730241, 2033107863601153, 2033108266254337])
        memorySafe = true
    elseif (numVars, prime) == (4, 13)
        primeArray = UInt.([4089429488566273, 4089440225984513, 4089445594693633, 4089451231838209, 4089452037144577])
        memorySafe = true
    else
        throw(ArgumentError("I haven't figured out bounds for this yet!"))
    end

    resultTotalDegree = numVars * (prime - 1) * prime
    key = resultTotalDegree + 1
    fftLen = Base._nextpow2(resultTotalDegree * key^(numVars - 2) + 1)
    
    nttPowPlans = GPUPolynomials.NTTPowPlan[]
    for p in primeArray
        nttPowPlan = GPUPolynomials.NTTPowPlan(fftLen, prime, p)
        push!(nttPowPlans, nttPowPlan)
    end
    resultDataType = GPUPolynomials.get_uint_type(max(Base._nextpow2(Int(ceil(log2(prod(BigInt.(primeArray)))))), 32))
    crtPlan = GPUPolynomials.plan_crt(resultDataType.(primeArray))

    return Δ₁Plan(numVars, prime, key, fftLen, resultTotalDegree, primeArray, eltype(primeArray), nttPowPlans, crtPlan, memorySafe)
end

function Δ₁(g::CufpMPolyRingElem)
    if !(g.opPlan isa Δ₁Plan)
        throw(ArgumentError("Input polynomial needs an OperationPlan!"))
    end

    if g.opPlan.memorySafe
        throw("")
    else
        memory_unsafe_Δ₁(g)
    end
end

function memory_unsafe_Δ₁(g::CufpMPolyRingElem)
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
    resultCoeffs .÷= eltype(resultCoeffs)(g.opPlan.prime)
    resultCoeffs .%= eltype(resultCoeffs)(g.opPlan.prime)
    resultCoeffs = UInt32.(resultCoeffs)

    resultDegs = GPUPolynomials.kronecker_to_bitpacked(encodedDegs, g.opPlan.key, numVars, g.opPlan.totalDegree, g.bits, UInt)

    return CufpMPolyRingElem(resultCoeffs, resultDegs, g.bits, true, g.opPlan.totalDegree, g.parent, GPUPolynomials.EmptyPlan())
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