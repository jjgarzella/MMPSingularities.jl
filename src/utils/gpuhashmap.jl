struct MyMap{A, B}
    buckets::B
    keys::A
    values::B
end

function find_biggest_prime_lt(num)
    if num % 2 == 0
        temp = num - 1
    else
        temp = num
    end

    while true
        if isprime(temp)
            return temp
        end
        temp -= 2
    end
end

function Base.getindex(mymap::MyMap, key)
    bucket = (key % length(mymap.buckets)) + 1
    startidx = mymap.buckets[bucket]

    while key != mymap.keys[startidx]
        startidx += 1
    end

    return mymap.values[startidx]
end

Adapt.@adapt_structure MyMap

function make_dict(encodedMons)
    num = find_biggest_prime_lt(length(encodedMons) >> 1)
    mods = [UInt[] for i in 1:num]
    values = [Int[] for i in 1:num]

    for i in eachindex(encodedMons)
        idx = (encodedMons[i] % num) + 1
        push!(mods[idx], encodedMons[i])
        push!(values[idx], i)
    end

    startindices = zeros(Int, num)
    curridx = 1

    hashedMons = zeros(eltype(encodedMons), length(encodedMons))
    originalIndices = zeros(Int, length(encodedMons))
    for i in eachindex(mods)
        startindices[i] = curridx
        for j in eachindex(mods[i])
            hashedMons[curridx] = mods[i][j]
            originalIndices[curridx] = values[i][j]
            curridx += 1
        end
    end

    reverseDict = MyMap(CuArray(startindices), CuArray(hashedMons), CuArray(originalIndices))

    return reverseDict
end