

function nverts(edgearray)
    if !(eltype(edgearray) <: Array)
        edgearray = collect.(edgearray)
    end

    maximum(reduce(hcat,edgearray))
end

"""
creates the binomial edge ideal over F_p corresponding to the graph
(v,e)

p - prime power
v - set of vertices
e - set of edges, expected to be a collection of 2-element collections

"""
function binomial_edge_ideal(p,v,e)
    R, xs, ys = polynomial_ring(GF(p),:x => 1:length(v),:y => 1:length(v))

    gens = zeros(R,0)
    for i in 1:length(e)
        edge = e[i]
        f = xs[edge[1]]*ys[edge[2]] - ys[edge[1]]*xs[edge[2]]
        push!(gens,f)
    end
    ideal(R,gens)
end

function bracketpower(I,m)
    R = base_ring(I)

    res = zeros(R,0)
    for g in gens(I)
        new_gen = g^m
        push!(res,new_gen)
    end

    ideal(R,res)
end

function fedderColonIdeal(I)
    p = characteristic(base_ring(I))
    bracketpower(I,p):I
end

function fedder_quotient(I)
    p = characteristic(base_ring(I))
    (bracketpower(I,p), I)
end

function ideal_quotient_example_text(name,I,J)
    result = "$name\n\n"

    result *= "IDEAL QUOTIENT I:J\n"

    R = base_ring(I)
    vars = gens(R)
    result *= "$vars\n"
    p = characteristic(R)
    result *= "$p\n"

    result *= "I = (\n"

    for g in gens(I)
        result *= "$g,\n"
    end

    result *= ")\n"

    result *= "J = (\n"

    for g in gens(J)
        result *= "$g,\n"
    end

    result *= ")\n"
end

function fedder_quotient_example_text(name,I)
    ideal_quotient_example_text(name,fedder_quotient(I)...)
end

function append_fedder(filename,name,I)
    open(filename, "a") do file
        write(file,fedder_quotient_example_text(name,I))
    end
end


