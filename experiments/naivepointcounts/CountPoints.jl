
"""
TODO: get this into GPUPolynomials.jl

point - vector of the entries of the point
coef - coefficients of the polynomial
exp_vecs - length(point) x length(coefs) matrix describing the
   expoenent vectors
   (in a kernel this could be in shared memory)
"""
function naive_eval(point,coefs,exp_vecs)
    n = length(point) #number of variables

    res = zero(eltype(point))
    for t in 1:length(coefs)
        val = coefs[t]
        for i in 1:n
            val *= point[i] ^ exp_vecs[i,t]
        end
        res += val
    end

    res
end

using AcceleratedKernels

"""
Here, coefs is a vector whose entries
are the coefficients of the single polynomial
"""
function pointcount(points,coefs,exp_vecs)
    fq_points = similar(exp_vecs,size(points,1))

    AcceleratedKernels.foraxes(points,1) do i
        point = points[i,:]
        val = naive_eval(point,coefs,exp_vecs)
        fq_points[i] = val == 0 ? 1 : 0
    end
   
    AcceleratedKernels.reduce(+, fq_points; init=zero(eltype(fq_points)))
end

"""
Here, coefs is a vector whose entries
are the coefficients of the single polynomial
"""
function pointcount(points,coefs,exp_vecs,p)
    fq_points = similar(exp_vecs,size(points,1))

    AcceleratedKernels.foraxes(points,1) do i
        point = points[i,:]
        val = naive_eval(point,coefs,exp_vecs)
        val %= p
        fq_points[i] = val == 0 ? 1 : 0
    end

    AcceleratedKernels.reduce(+, fq_points; init=zero(eltype(fq_points)))
end

"""
Here, coefs is a matrix whose columns
are the coefficients of each polynomial
"""
function pointcounts(points,coefs,exp_vecs)
    fq_points = similar(exp_vecs,size(points,1),size(coefs,2))

    AcceleratedKernels.foraxes(points,1) do i
        point = points[i,:]
        for j in 1:size(coefs,2)
            val = naive_eval(point,coefs[:,j],exp_vecs)
            fq_points[i,j] = val == 0 ? 1 : 0
        end
    end

    AcceleratedKernels.reduce(+,fq_points; init=zero(eltype(fq_points)),dims=1)
end

"""
Here, coefs is a matrix whose columns
are the coefficients of each polynomial

manually mods by p

also assumes that the polynomial is small enough that this 
won't overflow 2^24
"""
function pointcounts(points,coefs,exp_vecs,p)
    fq_points = similar(exp_vecs,size(points,1),size(coefs,2))

    AcceleratedKernels.foraxes(points,1) do i
        point = points[i,:]
        for j in 1:size(coefs,2)
            val = naive_eval(point,coefs[:,j],exp_vecs)
            val %= p
            fq_points[i,j] = val == 0 ? 1 : 0
        end
    end

    AcceleratedKernels.reduce(+,fq_points; init=zero(eltype(fq_points)),dims=1)
end

function main_cpu(n)
    exp_vec_ints = [4 3 2 1 0 3 2 1 0 2 1 0 1 0 0 3 2 1 0 2 1 0 1 0 0 2 1 0 1 0 0 1 0 0 0; 0 1 2 3 4 0 1 2 3 0 1 2 0 1 0 0 1 2 3 0 1 2 0 1 0 0 1 2 0 1 0 0 1 0 0; 0 0 0 0 0 1 1 1 1 2 2 2 3 3 4 0 0 0 0 1 1 1 2 2 3 0 0 0 1 1 2 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 4]

    rand_k3 = rand(0:2,35)
    rand_k3_list = rand(0:2,35,n)
    println("surface: $rand_k3")

    points_f3 = projective_space_points(3,1,3)
    exp_vecs_f3 = convert.(Float32,exp_vec_ints)
    rand_k3_f3 = convert.(Float32,rand_k3)
    rand_k3_f3_list = convert.(Float32,rand_k3_list)

    points_f9 = projective_space_points(3,2,3)
    F9 = eltype(points_f9)
    exp_vecs_f9 = Int32.(exp_vec_ints)
    rand_k3_f9 = F9.(rand_k3)
    rand_k3_f9_list = F9.(rand_k3_list)

    points_f27 = projective_space_points(3,3,3)
    exp_vecs_f27 = Int32.(exp_vec_ints)
    F27 = eltype(points_f27)
    rand_k3_f27 = F27.(rand_k3)
    rand_k3_f27_list = F27.(rand_k3_list)

    println("testing one surface over F3")
    @time pc = pointcount(points_f3,rand_k3_f3,exp_vecs_f3,3)
    println(pc)

    println("testing $n surfaces over F3")
    @time pcs = pointcounts(points_f3,rand_k3_f3_list,exp_vecs_f3,3)
    println(pcs)

    println("testing one surface over F9")
    @time pc = pointcount(points_f9,rand_k3_f9,exp_vecs_f9)
    println(pc)

    println("testing $n surfaces over F9")
    @time pcs = pointcounts(points_f9,rand_k3_f9_list,exp_vecs_f9)
    println(pcs)

    println("testing one surface over F27")
    @time pc = pointcount(points_f27,rand_k3_f27,exp_vecs_f27)
    println(pc)
end

using CUDA

function main_cuda(n)
    exp_vec_ints = [4 3 2 1 0 3 2 1 0 2 1 0 1 0 0 3 2 1 0 2 1 0 1 0 0 2 1 0 1 0 0 1 0 0 0; 0 1 2 3 4 0 1 2 3 0 1 2 0 1 0 0 1 2 3 0 1 2 0 1 0 0 1 2 0 1 0 0 1 0 0; 0 0 0 0 0 1 1 1 1 2 2 2 3 3 4 0 0 0 0 1 1 1 2 2 3 0 0 0 1 1 2 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 4]

    rand_k3 = rand(0:2,35)
    rand_k3_list = rand(0:2,35,n)
    println("surface: $rand_k3")

    points_f3 = projective_space_points(3,1,3)
    exp_vecs_f3 = convert.(Float32,exp_vec_ints)
    rand_k3_f3 = convert.(Float32,rand_k3)
    rand_k3_f3_list = convert.(Float32,rand_k3_list)

    exp_vecs = CuArray(Int32.(exp_vec_ints))

    points_f9 = projective_space_points(3,2,3)
    F9 = eltype(points_f9)
    rand_k3_f9 = CuArray(F9.(rand_k3))
    rand_k3_f9_list = CuArray(F9.(rand_k3_list))

    points_f27 = projective_space_points(3,3,3)
    F27 = eltype(points_f27)
    rand_k3_f27 = CuArray(F27.(rand_k3))
    rand_k3_f27_list = CuArray(F27.(rand_k3_list))

    #println("testing one surface over F3")
    #@time pc = pointcount(points_f3,rand_k3_f3,exp_vecs_f3,3)
    #println(pc)

    #println("testing $n surfaces over F3")
    #@time pcs = pointcounts(points_f3,rand_k3_f3_list,exp_vecs_f3,3)
    #println(pcs)

    println("testing one surface over F9")
    @time pc = pointcount(points_f9,rand_k3_f9,exp_vecs)
    println(pc)

    println("testing $n surfaces over F9")
    @time pcs = pointcounts(points_f9,rand_k3_f9_list,exp_vecs)
    println(pcs)

    println("testing one surface over F27")
    @time pc = pointcount(points_f27,rand_k3_f27,exp_vecs_f27)
    println(pc)
end
p

