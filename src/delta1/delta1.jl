nvars(x) = GPUPolynomials.nvars(x)

include("int128stuff.jl")

struct Delta1Pregen
    numVars::Int
    prime::Int
    key::Int
    encodedLen::Int
    totalDegree::Int
    gpupregen::GPUPowPregen
end

function serialize_pregen(n, p, restricted = false)
    if !restricted
        pregen = pregen_delta1_unrestricted(n, p)
    else
        pregen = pregen_delta1_restricted(n, p)
    end
    temp = CUDA.zeros(Int, 1)
    pregen.gpupregen.nttpregen.butterfly = temp
    pregen.gpupregen.inttpregen.nttpregen.butterfly = temp
    
    str = restricted ? "restricted_pregen/" : "unrestricted_pregen/"
    open("src/delta1/$str$(n)_$(p).jls", "w") do io
        serialize(io, pregen)
    end
end

function deserialize_pregen(n, p, restricted = false)
    str = restricted ? "restricted_pregen/" : "unrestricted_pregen/"
    pregen = open("src/delta1/$str$(n)_$(p).jls", "r") do io
        deserialize(io)
    end
    butterfly = generate_butterfly_permutations(pregen.gpupregen.nttpregen.len)
    pregen.gpupregen.nttpregen.butterfly = butterfly
    pregen.gpupregen.inttpregen.nttpregen.butterfly = butterfly

    return pregen
end

function pregen_delta1(numVars, prime, restricted = false)
    return deserialize_pregen(numVars, prime, restricted)
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
        primeArray = UInt.([2033101286932481, 2033107326730241, 2033107863601153, 2033108266254337])
    elseif (numVars, prime) == (4, 13)
        primeArray = UInt.([4089429488566273, 4089440225984513, 4089445594693633, 4089451231838209, 4089452037144577])
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
    coeffs = get_coeffs(intermediate)
    degrees = copy(get_exps(intermediate))
    bits = intermediate.poly.data.bits
    mask = (one(eltype(coeffs)) << bits) - 1

    keyPowers = [key ^ i for i in 0:nvars(intermediate) - 2]

    indices = zeros(Int, length(coeffs))
    subtract = zeros(eltype(coeffs), length(coeffs))
    degrees .*= p
    for i in eachindex(degrees)
        resultidx = 1
        deg = degrees[i]
        for i in 1:nvars(intermediate) - 1
            resultidx += (deg & mask) * keyPowers[i]
            deg >>= bits
        end
        indices[i] = resultidx
        subtract[i] = coeffs[i] ^ p
    end

    return indices, subtract
end

function remove_pth_power_terms(intermediate::HomogeneousPolynomial, key::Int, vec, p, primeArray)
    removeindices, subtract = generate_remove_indices(intermediate, key, p)
    if vec isa CuArray
        gpu_remove_pth_power_terms(CuArray(removeindices), CuArray(subtract), vec, primeArray)
    else
        cpu_remove_pth_power_terms(removeindices, subtract, vec, Array(primeArray))
    end

    return nothing
end

@inline function sub_mod(x::Signed, y::Signed, m::Signed)
    return mod(x - y, m)
end

@inline function sub_mod(x::Unsigned, y::Unsigned, m::Unsigned)
    if y > x
        return m - mod(y - x, m)
    else
        return mod(x - y, m)
    end
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

function memorysafe_delta1(intermediate::HomogeneousPolynomial, prime::Int; pregen::Delta1Pregen)
    numVars = nvars(intermediate)
    if pregen === nothing
        pregen = pregen_delta1(numVars, prime)
    end
    @assert (pregen.numVars, pregen.prime) == (numVars, prime)

    vect = kronecker_substitution(intermediate, pregen.key, pregen.encodedLen, eltype(pregen.gpupregen.primeArray))

    multimoddenseresult = memorysafe_gpu_ntt_pow(vect, prime; pregen = pregen.gpupregen, docrt = false)

    remove_pth_power_terms(intermediate, pregen.key, multimoddenseresult, prime, pregen.gpupregen.primeArray)

    multimodresultCoefs, encodedDegs = sparsify(multimoddenseresult)

    resultCoefs = build_result(multimodresultCoefs, pregen.gpupregen.crtpregen, pregen.gpupregen.resultType)

    resultCoefs .÷= eltype(resultCoefs)(prime)
    resultCoefs .%= eltype(resultCoefs)(prime)
    resultCoefs = UInt64.(resultCoefs)

    resultDegs = decode_kronecker_substitution(encodedDegs, pregen.key, nvars(intermediate), intermediate.homogDegree * prime)
    resultCoefs = Array(resultCoefs)

    result = zero(intermediate.poly.parent)
    result.data = Oscar.fpMPolyRingElem(intermediate.poly.parent.data, resultCoefs, resultDegs)

    return HomogeneousPolynomial(result)
end

function memoryunsafe_delta1(intermediate::HomogeneousPolynomial, prime::Int; pregen::Delta1Pregen)
    numVars = nvars(intermediate)
    if pregen === nothing
        pregen = pregen_delta1(numVars, prime)
    end
    @assert (pregen.numVars, pregen.prime) == (numVars, prime)

    vect = kronecker_substitution(intermediate, pregen.key, pregen.encodedLen, eltype(pregen.gpupregen.primeArray))

    multimoddenseresult = gpu_ntt_pow(vect, prime; pregen = pregen.gpupregen, docrt = false)

    remove_pth_power_terms(intermediate, pregen.key, multimoddenseresult, prime, pregen.gpupregen.primeArray)
    
    multimodresultCoefs, encodedDegs = sparsify(multimoddenseresult)

    resultCoefs = build_result(multimodresultCoefs, pregen.gpupregen.crtpregen, pregen.gpupregen.resultType)
    resultCoefs .÷= eltype(resultCoefs)(prime)
    resultCoefs .%= eltype(resultCoefs)(prime)
    resultCoefs = UInt64.(resultCoefs)

    resultDegs = decode_kronecker_substitution(encodedDegs, pregen.key, nvars(intermediate), intermediate.homogDegree * prime)
    resultCoefs = Array(resultCoefs)

    result = zero(intermediate.poly.parent)
    result.data = Oscar.fpMPolyRingElem(intermediate.poly.parent.data, resultCoefs, resultDegs)

    return HomogeneousPolynomial(result)
end

function delta1(intermediate::HomogeneousPolynomial, prime::Int; pregen::Delta1Pregen)
    if nvars(intermediate) == 4 
        if prime in [2, 3, 5, 7]
            return memoryunsafe_delta1(intermediate, prime; pregen = pregen)
        else
            return memorysafe_delta1(intermediate, prime; pregen = pregen)
        end
    else 
        throw(ArgumentError("Haven't implemented this yet"))
    end
end

function decode_kronecker_substitution(encodedDegs, key, numVars, totalDegree)
    result = CUDA.zeros(UInt64, numVars, length(encodedDegs))

    kernel = @cuda launch=false decode_kronecker_substitution_kernel!(encodedDegs, key, numVars, totalDegree, result)
    config = launch_configuration(kernel.fun)
    threads = min(length(encodedDegs), config.threads)
    blocks = cld(length(encodedDegs), threads)

    kernel(encodedDegs, key, numVars, totalDegree, result; threads = threads, blocks = blocks)

    return Array(result)
end

function decode_kronecker_substitution_kernel!(encodedDegs::CuDeviceVector, key::Int, numVars::Int, totalDegree::Int, dest::CuDeviceMatrix)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if idx <= length(encodedDegs)
        num = encodedDegs[idx] - 1
        for i in numVars:-1:2
            num, r = divrem(num, key)
            dest[i, idx] = r
            totalDegree -= r
        end
        dest[1, idx] = totalDegree
    end

    return nothing
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

function cpu_kronecker_substitution(hp::HomogeneousPolynomial, key::Int, len::Int, nttType::DataType)
    result = zeros(nttType, len)
    coeffs = get_coeffs(hp)
    exps = get_exps(hp)
    keyPowers = [key ^ i for i in 0:nvars(hp) - 2]
    bits = hp.poly.data.bits
    mask = (one(eltype(exps)) << bits) - 1

    for term in eachindex(exps)
        resultIdx = 1
        deg = exps[term]
        for i in 1:nvars(hp) - 1
            resultIdx += (deg & mask) * keyPowers[i]
            deg >>= bits
        end

        result[resultIdx] = nttType(coeffs[term])
    end

    return result
end

function kronecker_substitution(hp::HomogeneousPolynomial, key::Int, len::Int, nttType::DataType)
    result = CUDA.zeros(nttType, len)
    coeffs = CuArray(get_coeffs(hp))
    exps = CuArray(get_exps(hp))
    keyPowers = CuArray([key ^ i for i in 0:nvars(hp) - 2])
    bits = hp.poly.data.bits
    mask = (one(eltype(exps)) << bits) - 1

    kernel = CUDA.@sync @cuda launch=false kronecker_substitution_kernel!(coeffs, exps, keyPowers, mask, bits, result, nttType)
    config = launch_configuration(kernel.fun)
    threads = min(length(coeffs), config.threads)
    blocks = cld(length(coeffs), threads)

    kernel(coeffs, exps, keyPowers, mask, bits, result, nttType; threads = threads, blocks = blocks)

    return result
end

function kronecker_substitution_kernel!(coeffs, exps, keyPowers, mask, bits::Int, result::CuDeviceVector, nttType::DataType)
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    @inbounds if idx <= length(coeffs)
        resultIdx = 1
        deg = exps[idx]
        for i in eachindex(keyPowers)
            resultIdx += (deg & mask) * keyPowers[i]
            deg >>= bits
        end

        result[resultIdx] = nttType(coeffs[idx])
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