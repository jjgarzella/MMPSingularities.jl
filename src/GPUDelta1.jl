module GPUDelta1

export Delta1Pregen, pregen_delta1, pregen_delta1_unrestricted, delta1, sort_to_kronecker_order

using CUDA
using GPUPolynomials

struct Delta1Pregen
    numVars::Int
    prime::Int
    key::Int
    encodedLen::Int
    totalDegree::Int
    gpupregen::GPUPowPregen
end

function pregen_delta1(numVars, prime)
    if (numVars, prime) == (4, 2)
        return pregen_delta1_unrestricted(numVars, prime)
    elseif (numVars, prime) == (4, 3)
        return pregen_delta1_unrestricted(numVars, prime)
    elseif (numVars, prime) == (4, 5)
        return pregen_delta1_unrestricted(numVars, prime)
    elseif (numVars, prime) == (4, 7)
        return pregen_delta1_restricted(numVars, prime)
    elseif (numVars, prime) == (4, 11)
        return pregen_delta1_restricted(numVars, prime)
    else
        throw(ArgumentError("I haven't figured out bounds for this yet!"))
    end
end

function pregen_delta1_unrestricted(numVars, prime)
    if (numVars, prime) == (4, 2)
        primeArray = UInt.([12289])
    elseif (numVars, prime) == (4, 3)
        primeArray = UInt.([114689])
    elseif (numVars, prime) == (4, 5)
        primeArray = UInt.([13631489, 23068673])
    elseif (numVars, prime) == (4, 7)
        primeArray = UInt.([167772161, 377487361, 469762049])
    elseif (numVars, prime) == (4, 11)
        primeArray = UInt.([2281701377, 3221225473, 3489660929, 3892314113, 7918845953, 8858370049])
    else
        throw(ArgumentError("I haven't figured out bounds for this yet!"))
    end

    pregentime = CUDA.@timed begin
        resultTotalDegree = numVars * (prime - 1) * prime
        key = resultTotalDegree + 1
        encodedLen = numVars * (prime - 1) * key^(numVars - 2) + 1
        gpupregen = pregen_gpu_pow(primeArray, get_fft_size(encodedLen, prime))
    end
    # println("Delta1Pregen took $(pregentime.time) s to pregenerate")

    return Delta1Pregen(numVars, prime, key, encodedLen, resultTotalDegree, gpupregen)
end

function pregen_delta1_restricted(numVars, prime)
    if (numVars, prime) == (4, 2)
        primeArray = UInt.([12289])
    elseif (numVars, prime) == (4, 3)
        primeArray = UInt.([114689])
    elseif (numVars, prime) == (4, 5)
        primeArray = UInt.([13631489, 23068673])
    elseif (numVars, prime) == (4, 7)
        primeArray = UInt.([167772161, 377487361, 469762049])
    elseif (numVars, prime) == (4, 11)
        primeArray = UInt.([2281701377, 3221225473, 3489660929, 3892314113, 7918845953, 8858370049])
    else
        throw(ArgumentError("I haven't figured out bounds for this yet!"))
    end

    resultTotalDegree = numVars * (prime - 1) * prime
    key = (numVars - 1) * prime * (prime - 1) + 1
    encodedLen = (numVars - 1) * (prime - 1) * key^(numVars - 2) + (prime - 1) * key^(numVars - 3) + 1
    gpupregen = pregen_gpu_pow(primeArray, get_fft_size(encodedLen, prime))

    return Delta1Pregen(numVars, prime, key, encodedLen, resultTotalDegree, gpupregen)
end

function get_fft_size(veclength::Int, pow)
    finalLength = (veclength - 1) * pow + 1
    return Base._nextpow2(finalLength)
end

function generate_remove_indices(intermediate, key, p)
    indices = zeros(Int, length(intermediate.coeffs))
    subtract = zeros(eltype(intermediate.coeffs), length(intermediate.coeffs))
    intermediate.degrees .*= p
    for i in eachindex(intermediate.coeffs)
        indices[i] = homogkron(view(intermediate.degrees, :, i), key)
        subtract[i] = intermediate.coeffs[i] ^ p
    end
    return CuArray(indices), CuArray(subtract)
end

function gpu_remove_pth_power_terms(intermediate, key, vec, p)
    removeindices, subtract = generate_remove_indices(intermediate, key, p)

    kernel = @cuda launch=false gpu_remove_pth_power_terms_kernel!(removeindices, subtract, vec)
    config = launch_configuration(kernel.fun)
    threads = min(length(removeindices), config.threads)
    blocks = cld(length(removeindices), threads)

    kernel(removeindices, subtract, vec; threads = threads, blocks = blocks)

    return nothing
end

function gpu_remove_pth_power_terms_kernel!(removeindices, subtract, vec)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if idx <= length(removeindices)
        vec[removeindices[idx]] -= subtract[idx]
    end

    return nothing
end

function delta1(intermediate::HomogeneousPolynomial, prime::Int; pregen = nothing)
    numVars = size(intermediate.degrees, 1)
    if pregen === nothing
        pregen = pregen_delta1(numVars, prime)
    end
    @assert (pregen.numVars, pregen.prime) == (numVars, prime)

    vec = kronecker_substitution(intermediate, pregen.key, pregen.encodedLen, eltype(pregen.gpupregen.primeArray))
    # cpuvec = cpu_kronecker_substitution(intermediate, pregen.key, pregen.encodedLen, eltype(pregen.gpupregen.primeArray))
    # @assert Array(vec) == cpuvec
    result = gpu_ntt_pow(vec, prime; pregen = pregen.gpupregen)

    gpu_remove_pth_power_terms(intermediate, pregen.key, result, prime)

    hp_result = decode_kronecker_substitution(result, pregen.key, pregen.numVars, pregen.totalDegree)    

    for i in eachindex(hp_result.coeffs)
        if hp_result.coeffs[i] % prime != 0
            println("$(hp_result.coeffs[i]), $(hp_result.degrees[:, i]) isn't equal to 0 mod $prime !")
        end
    end
    hp_result.coeffs .= map(x -> x ÷ prime, hp_result.coeffs)
    hp_result.coeffs .= map(x -> mod(x, prime), hp_result.coeffs)
    hp_result.coeffs = eltype(intermediate.coeffs).(hp_result.coeffs)

    return hp_result
end

function cpu_decode_kronecker_substitution(vec, key, numVars, totalDegree)
    flags = map(x -> x != 0 ? 1 : 0, vec)
    indices = accumulate(+, flags)

    resultLen = indices[end]

    resultCoeffs = zeros(eltype(vec), resultLen)
    resultDegrees = zeros(Int, resultLen, numVars)

    for i in eachindex(vec)
        if i != 0
            return false
        end
    end
end

function decode_kronecker_substitution(vec, key, numVars, totalDegree)
    flags = map(x -> x != 0 ? 1 : 0, vec)
    indices = accumulate(+, flags)

    CUDA.@allowscalar resultLen = indices[end]

    resultCoeffs = CUDA.zeros(eltype(vec), resultLen)
    resultDegrees = CUDA.zeros(Int, numVars, resultLen)

    kernel = @cuda launch = false decode_kronecker_kernel!(resultCoeffs, resultDegrees, vec, flags, indices, key, numVars, totalDegree)
    config = launch_configuration(kernel.fun)
    threads = min(length(vec), config.threads)
    blocks = cld(length(vec), threads)

    kernel(resultCoeffs, resultDegrees, vec, flags, indices, key, numVars, totalDegree; threads = threads, blocks = blocks)

    return HomogeneousPolynomial(Array(resultCoeffs), Array(resultDegrees))
end

function decode_kronecker_kernel!(resultCoeffs, resultDegrees, vec, flags, indices, key, numVars, totalDegree)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if idx <= length(vec)
        if flags[idx] != 0
            resultidx = indices[idx]
            num = idx - 1
            invhomogkron(num, key, numVars, view(resultDegrees, :, resultidx), totalDegree)
            resultCoeffs[resultidx] = vec[idx]
        end
    end

    return nothing
end

function cpu_kronecker_substitution(hp, key, len, nttType)
    result = zeros(nttType, len)
    for term in axes(hp.degrees, 2)
        resultidx = homogkron(view(hp.degrees, :, term), key)
        result[resultidx] = nttType(hp.coeffs[term])
    end

    return result
end

function kronecker_substitution(hp, key, len, nttType)
    result = CUDA.zeros(nttType, len)
    coeffs = CuArray(hp.coeffs)
    degrees = CuArray(hp.degrees)

    kernel = CUDA.@sync @cuda launch=false kronecker_substitution_kernel!(coeffs, degrees, key, result, nttType)
    config = launch_configuration(kernel.fun)
    threads = min(length(coeffs), config.threads)
    blocks = cld(length(coeffs), threads)

    kernel(coeffs, degrees, key, result, nttType; threads = threads, blocks = blocks)
    return result
end

function kronecker_substitution_kernel!(coeffs, degrees, key::Int, result::CuDeviceVector, nttType::DataType)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if idx <= length(coeffs)
        resultidx = homogkron(view(degrees, :, idx), key)
        result[resultidx] = nttType(coeffs[idx])
    end

    return nothing
end

function sort_to_kronecker_order(hp, key)
    encoded = encode_degrees(hp.degrees, key, true)

    perm = sortperm(encoded)
    hp.coeffs .= hp.coeffs[perm]
    hp.degrees .= hp.degrees[:, perm]
end

function cpu_remove_pth_power_terms!(big,small,p)
    i = 1
    k = 1

    n = size(small.degrees,1)

    smalldegs = zeros(eltype(small.degrees),n)
    bigdegs = zeros(eltype(big.degrees),n)

    function setslice_noalloc!(target,source,k)
        for j = 1:n
            target[j] = source[j,k]
        end
    end

    while i ≤ length(small.coeffs)
        # for a standard addition, remove the p
        
        #smalldegs = p .* small.degrees[i,:]
        for j = 1:n
            smalldegs[j] = p * small.degrees[j,k]
        end

        #bigdegs = big.degrees[k,:]
        setslice_noalloc!(bigdegs,big.degrees,k)
        while smalldegs != bigdegs
            k = k + 1
            setslice_noalloc!(bigdegs,big.degrees,k)
        end
        # now we know that the term in row k of big is a pth power of 
        # the term in row i of small
        # this is a subtraction
        big.coeffs[k] -= (small.coeffs[i])^p
        # big.coeffs[k] = 0
        i = i + 1
    end

    nothing
end

end