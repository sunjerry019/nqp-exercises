#!/usr/bin/env julia

include("mps.jl")

using .MatrixProductStates
using Test
using LinearAlgebra

@testset "Creation Functions" verbose=true begin
    @testset "Sequential Matrix" begin
        alpha = abs(rand(Int, 1)[1] % 10) + 1
        beta  = abs(rand(Int, 1)[1] % 10) + 1
        d     = abs(rand(Int, 1)[1] % 10) + 1

        M = create_seq_matrix_set(Complex{Float64}, alpha, beta, d)
        
        @test size(M) == (alpha, beta, d)

        flat = reshape(M, :)
        @test flat[1] == 0.0 + 0.0im
        @test all(diff(flat) .== (1.0+0.0im)) # https://stackoverflow.com/a/56799012
    end
end

@testset "Left" verbose=true begin
    @testset "Fusing" begin
        m = abs(rand(Int, 1)[1] % 10) + 1
        n = abs(rand(Int, 1)[1] % 10) + 1
        d = abs(rand(Int, 1)[1] % 10) + 1

        M = create_seq_matrix_set(Integer, m, n, d)
        Q = fuse_left(M)
        
        @test size(Q) == (m*d, n)

        # Random matrix element
        _alpha = abs(rand(Int, 1)[1] % m) + 1
        _beta  = abs(rand(Int, 1)[1] % n) + 1
        _k     = abs(rand(Int, 1)[1] % d) + 1

        # Since Julia is 1-indexed, (k,alpha) = (k-1)*m + (a-1) + 1 = (k-1)*m + a
        @test M[_alpha, _beta, _k] == Q[(_k-1)*m + _alpha, _beta]
    end
end

@testset "Right" verbose=true begin
    @testset "Fusing" begin
        m = abs(rand(Int, 1)[1] % 10) + 1
        n = abs(rand(Int, 1)[1] % 10) + 1
        d = abs(rand(Int, 1)[1] % 10) + 1

        M = create_seq_matrix_set(Integer, m, n, d)
        Q = fuse_right(M)
        
        @test size(Q) == (m, n*d)

        # Random matrix element
        _alpha = abs(rand(Int, 1)[1] % m) + 1
        _beta  = abs(rand(Int, 1)[1] % n) + 1
        _k     = abs(rand(Int, 1)[1] % d) + 1

        # Since Julia is 1-indexed, (k,b) = (k-1)*n + (b-1) + 1 = (k-1)*n + b
        @test M[_alpha, _beta, _k] == Q[_alpha, (_k-1)*n +  _beta]
    end
end