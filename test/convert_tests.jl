using Oscar

function get_exps(poly::FqMPolyRingElem)
    expsDataType = UInt
    expsPtr = Base.unsafe_convert(Ptr{expsDataType}, poly.data.exps)
    expsVec = unsafe_wrap(Vector{expsDataType}, expsPtr, poly.data.length)

    return expsVec
end


function run()
    R, (x, y, z, w) = polynomial_ring(GF(5), 4)

    f = x^4 + y^3 + z^2 + w

    display(get_exps(f))
end

run()