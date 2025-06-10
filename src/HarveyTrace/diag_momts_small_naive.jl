
function diag_momts_naive_little(f,p)
    g = f^(p-1)

    d = total_degree(f)

    ws = wics(d*(p-1),n)

    nMons = 35
    diag = zeros(base_ring(parent(f)),nMons)
    evs = zeros(Int,nMons,4)

    j = 1 
    for i in 1:size(ws,2)
        wic = ws[:,i]
        
        # the ones that match and are in the diagonal are the
        # things of the form pu - u, aka (p-1)*u
        match = true
        for i in 1:n
            if !divides(wic[i],p-1)[1]
                match = false
            end
        end

        if match
            v = divexact.(wic,p-1)
            evs[j,:] .= v
            #push!(evs,v)
            diag[j] = coeff(g,wic)
            #push!(diag,coeff(g,wic))
            j += 1
        end
    end

    (diag,evs)
end

