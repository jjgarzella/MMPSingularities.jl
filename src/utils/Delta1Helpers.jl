function change_encoding(small::CuVector{T}, large::CuVector{T}, smallKey::Int, largeKey::Int, nvars::Int) where T<:Integer
    kernel = @cuda launch=false change_encoding_kernel!(small, large, smallKey, largeKey, nvars)
    config = launch_configuration(kernel.fun)
    threads = min(length(small), Base._prevpow2(config.threads))
    blocks = div(length(small), threads)

    kernel(small, large, smallKey, largeKey, nvars; threads = threads, blocks = blocks)
end

@inbounds function change_encoding_kernel!(small::CuDeviceVector{T}, large::CuDeviceVector{T}, smallKey::Int, largeKey::Int, nvars::Int) where T<:Integer
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if small[idx] != 0
        resultIdx = 1
        currLargeKey = 1
        smallIdx = idx - 1

        for i in 0:nvars - 2
            smallIdx, r = divrem(smallIdx, smallKey)
            resultIdx += r * currLargeKey
            currLargeKey *= largeKey
        end

        @cuassert resultIdx < length(large)
        large[resultIdx] = small[idx]
    end

    return
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
            arr[removeidx, p] = sub_mod(arr[removeidx, p], subtract[i] % primeArray[p], primeArray[p])
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
            vec[removeidx, p] = sub_mod(vec[removeidx, p], subtract[idx] % primeArray[p], primeArray[p])
        end
    end

    return nothing
end

function broadcast_mul!(v1::CuVector{T}, v2::CuVector{T}, m::CudaNTTs.Reducer{T}) where T<:Integer
    kernel = @cuda launch=false broadcast_mul_kernel!(v1, v2, m)
    config = launch_configuration(kernel.fun)
    threads = min(length(v1), Base._prevpow2(config.threads))
    blocks = div(length(v1), threads)

    kernel(v1, v2, m; threads = threads, blocks = blocks)

    return nothing
end

@inbounds function broadcast_mul_kernel!(v1::CuDeviceVector{T}, v2::CuDeviceVector{T}, m::CudaNTTs.Reducer{T}) where T<:Integer
    idx = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    
    v1[idx] = CudaNTTs.mul_mod(v1[idx], v2[idx], m)

    return nothing
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