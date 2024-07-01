using RowEchelon
using Test
using CSV, DataFrames

# As = Input Matrix, Rs = Expected RREF, Ps = Pivot Vectors (tentative)

As = Vector{Matrix{Int}}()
Rs = Vector{Matrix{Int}}()
Ps = Vector{Vector{Int}}()
# Issue #518 of Julia Base
push!(As, [ 1  2  0   3;
            2  4  0   7])
push!(Rs, [ 1  2  0   0;
            0  0  0   1])
push!(Ps, [ 1; 4])
# Example from Wikipedia
push!(As, [ 1  3  -1;
            0  1   7])
push!(Rs, [ 1  0 -22;
            0  1   7])
push!(Ps, [ 1; 2])
# Example from Wikibooks
push!(As, [ 0  3 -6   6  4  -5;
            3 -7  8  -5  8   9;
            3 -9 12  -9  6  15])
push!(Rs, [ 1  0 -2   3  0 -24;
            0  1 -2   2  0  -7;
            0  0  0   0  1   4])
push!(Ps, [ 1; 2; 5])
# Example from Rosetta Code
push!(As, [ 1  2 -1  -4;
            2  3 -1 -11;
           -2  0 -3  22])
push!(Rs, [ 1  0  0  -8;
            0  1  0   1;
            0  0  1  -2])
push!(Ps, [ 1; 2; 3])
# Magic Square
push!(As, [16  2  3  13;
            5 11 10   8;
            9  7  6  12;
            4 14 15   1])
push!(Rs, [ 1  0  0   1;
            0  1  0   3;
            0  0  1  -3;
            0  0  0   0])
push!(Ps, [ 1; 2; 3])

Aq = Vector{Matrix{Int}}()
Rq = Vector{Matrix{Int}}()

push!(Aq, [38962 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11726 0 0 0 3523;
0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10978 0 0;
0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10978 0 0 0;
0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10978 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 10978 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0;
0 0 0 0 0 0 0 0 0 0 0 14597 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 14597 0 0 0 0 0 0 0 0 0 0 0;
101761 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11396 0 0 0 4954;
0 0 0 0 0 0 14597 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 14597 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0;
96921 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5819 0 0 0 12192;
])

push!(Rq, [1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1])

@testset "Matrices of integers (treated as Float64)" begin
    for (A,R) in zip(As,Rs)
        C = rref(A)
        println(C)
        @test C isa Matrix{Float64}
        @test C ≈ R
    end
end

@testset "Matrices of rationals" begin
    for (A,R) in zip(As,Rs)
        B = Matrix{Rational{Int}}(A)
        C = rref(B)
        @test C isa Matrix{Rational{Int}}
        @test C == R
    end
end

@testset "Matrix of Complex numbers" begin
    A = [1+ im 1-im  4;
         3+2im 2-im  2]
    R = [1     0    -2- im;
         0     1     1+4im]
    for conv in [false, true]
        if conv
            A = Matrix{ComplexF64}(A)
        end
        for M in [A, Matrix(A')']
            C = rref(M)
            @test C isa Matrix{ComplexF64}
            @test C ≈ R
        end
    end
end

@testset "Matrices of integers (treated as Float64) with Pivots" begin
    for (A, X) in zip(As, zip(Rs,Ps))
        R=X[1]
        P=X[2]
        C = rref_with_pivots(A)
        @test C[1] isa Matrix{Float64}
        @test C[1] ≈ R
        @test C[2] == P
    end
end

@testset "Matrices of rationals with Pivots" begin
    for (A, X) in zip(As, zip(Rs,Ps))
        R=X[1]
        P=X[2]
        B = Matrix{Rational{Int}}(A)
        C = rref_with_pivots(B)
        @test C[1] isa Matrix{Rational{Int}}
        @test C[1] == R
        @test C[2] == P
    end
end

@testset "Matrix of Complex numbers with Pivots" begin
    A = [1+ im 1-im  4;
         3+2im 2-im  2]
    R = [1     0    -2- im;
         0     1     1+4im]
    for conv in [false, true]
        if conv
            A = Matrix{ComplexF64}(A)
        end
        for M in [A, Matrix(A')']
            C = rref_with_pivots(M)
            @test C[1] isa Matrix{ComplexF64}
            @test C[1] ≈ R
            @test C[2] == [1; 2]
        end
    end
end


# @testset "Quartic K3 test" begin
#     for (A,R) in zip(Aq,Rq)
#         C = rref(A)
#         println(C)
#         # @test C isa Matrix{Int64}
#         @test C ≈ R
#     end
# end

@testset "Quartic K3 test" begin
    for (A,R) in zip(Aq,Rq)
        # Measure the time taken by rref function
        t = @elapsed begin
            C = rref(A)
        end
        #println("Time taken by rref function: $t seconds")
        println("Time:",t)

        # Check if the result is of the correct type and approximately equal to the expected result
        println(C)
        @test C isa Matrix{Float64}
        @test C ≈ R
    end
end

using DelimitedFiles
A = readdlm("verylarge.csv", ',', Int64)

# df = DataFrame()
# CSV.read("output_file_1706742160023102800.csv", df)
# BigMatA = convert(Matrix, df)
# println(dump(BigMatA))
# println(BigMatA)

@testset "BIG Matrix" begin
    t = @elapsed begin
        C = rref(A)
    end
    #println("Time taken by rref function: $t seconds")
    println("Time:",t)
    writedlm( "resultB.csv",  C, ',')
end