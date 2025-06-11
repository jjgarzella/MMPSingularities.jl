using GaloisFields

function affine_line_points(p,r)
    if r == 1
        if 2^24 < p
            throw(ArgumentError("p is too big to fit inside a Float32"))
        end
        convert.(Float32,collect(0:p-1))
    else
        F = @GaloisField! p^r β
        collect(F)
    end
end

function projective_line_points(p,r)
   aff = affine_line_points(p,r)
   etype = eltype(aff)
   os = ones(etype,length(aff))

   [os aff;
    zero(etype) one(etype)]
end

function affine_space_points(p,r,n)
    aff = affine_line_points(p,r)
    if n == 1
        return aff
    end
    tuples = collect(Iterators.product([aff for i in 1:n]...))
    vectors = collect.(tuples)
    reduce(vcat,transpose.(vectors))
end

"""
Outputs a matrix whose rows contain
(representatives of) all of the projective points in 
\\mathbb{P}^n over \\mathbb{F}_p^r
"""
function projective_space_points(p,r,n)
    if n == 1
        return projective_line_points(p,r)
    end

    aff_n = affine_space_points(p,r,n)
    etype = eltype(aff_n)
    os = ones(etype,size(aff_n,1))

    proj_nminus1 = projective_space_points(p,r,n-1)
    zs = zeros(etype,size(proj_nminus1,1))

    [os aff_n;
     zs proj_nminus1]
end
