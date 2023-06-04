#!/usr/bin/env julia

include("mps.jl")

using .MatrixProductStates
using Test
using LinearAlgebra

# @testset "Left" verbose=true begin
#     @test "Fusing" begin
        
#     end
# end
M = create_seq_matrix_set(Complex{Float64}, 3, 4, 10)
display(size(M))
Q = fuse_left(M)
display(size(Q))
println(M[2,3,4])
println(Q[10,3])

