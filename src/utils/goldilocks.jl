const goldilocks = 0xffffffff00000001

struct GoldilocksReducer{UInt64} <: CudaNTTs.Reducer{UInt64}
    function GoldilocksReducer()
        return new{UInt64}()
    end
end

@inline function CudaNTTs.add_mod(x::UInt64, y::UInt64, m::GoldilocksReducer{UInt64})
    result = x + y

    return (result >= goldilocks || result < x) ? result - goldilocks : result
end

@inline function CudaNTTs.sub_mod(x::UInt64, y::UInt64, m::GoldilocksReducer{UInt64})
    if y > x
        return (goldilocks - y) + x
    else
        return x - y
    end
end

@inline function mul_wide(a::UInt64, b::UInt64)
    return Base.llvmcall(
        """
        %3 = zext i64 %0 to i128
        %4 = zext i64 %1 to i128
        %5 = mul i128 %3, %4
        %6 = trunc i128 %5 to i64
        %7 = lshr i128 %5, 64
        %8 = trunc i128 %7 to i64
        %9 = insertvalue [2 x i64] undef, i64 %6, 0
        %10 = insertvalue [2 x i64] %9, i64 %8, 1
        ret [2 x i64] %10
        """,
        Tuple{UInt64, UInt64},
        Tuple{UInt64, UInt64},
        a, b
    )
end

@inline function CudaNTTs.mul_mod(x::UInt64, y::UInt64, m::GoldilocksReducer{UInt64})
    lo, hi = mul_wide(x, y)

    middle = hi & 0xffffffff
    high = hi >>> 32
    
    low2 = lo - high
    if high > lo
        low2 += goldilocks
    end
    
    product = middle << 32
    product -= (product >>> 32)
    
    result = low2 + product
    
    if (result < product) || (result >= goldilocks)
        result -= goldilocks
    end
    
    return result
end