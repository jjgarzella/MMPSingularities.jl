

function frobeniusNu(p,e,f)
    n = 0

    ff = parent(f)(1)

    while !inPowerOfVariableIdeal(p,p^e,ff)
        n = n + 1
        ff *= f
    end

    n
end

"""
Gives the estimate in page 9 of arXiv:1906.09491
for a given nu and e, i.e.
nu_e(f) = nu
"""
function estimateFPureThreshold(p,e,nu)
    #[nu // (p^e - 1), (nu + 1) // p^e]
    nu / p^e
end

function estimateFPureThreshold(f,N)

end

function frobeniusNu_CI(p,e,fs)
    n = 0

    ffs = [parent(fs[1])(1)]

    function ideal_in_variable_ideal(generators)
        for g in generators
            if !inPowerOfVariableIdeal(p,p^e,g)
                return false
            end
        end
        true
    end

    while !ideal_in_variable_ideal(ffs)
        n = n+1
        println("Calculating I^$n...")
        new_ffs = []
        for g in ffs
            for h in fs
                push!(new_ffs,g*h)
            end
        end
        ffs = new_ffs
    end
    n
end

# complete intersections
function estimateFPureThreshold_CI(fs,N)

end
