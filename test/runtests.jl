
include("../src/MMPSingularities.jl")

using Oscar
using CUDA

using Test


function test_K3_5()

    best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(5,f,10)

    R, (x1,x2,x3,x4) = polynomial_ring(GF(5),4)

    ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4

    ftwo2 = x1^4 + x1^3*x2 + 2*x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + x1^2*x4^2 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 3*x1*x2^2*x4 + 2*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x4^3 + x2^3*x3 + x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 4*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + 2*x2*x4^3 + x3^2*x4^2 + x3*x4^3

    fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4

    ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3


    @test best_qfs_height(ftwo1) == 2
    @test best_qfs_height(ftwo2) == 2
    @test best_qfs_height(fthree1) == 3
    @test best_qfs_height(ffour1) == 4
end

function test_time_K3_5()
    
    qfs_height_fn = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu

    pregen = MMPSingularities.pregen_delta1(4,5)


    R, (x1,x2,x3,x4) = polynomial_ring(GF(5),4)

    ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

    # prime the jitter
    println("PRIMING THE JITTER:")
    CUDA.@time qfs_height_fn(5,ffour1,10,pregen)

    CUDA.@time qfs_height_fn(5,ffour1,10,pregen)

end


@testset "K3 surfaces" begin
    #test_K3_5()
    test_time_K3_5()
end
