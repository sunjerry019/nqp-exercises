#!/usr/bin/env julia

include("qm.jl")

using .QuantumMechanics
using Test
using Random
using LinearAlgebra

@testset "State Generation" verbose=true begin
    @testset "Random State Generation" begin
        L = abs(rand(Int, 1)[1] % 10)
        x = GetRandomState(L)
        @test x.L == L
        @test size(x.state) == (2^L,)
        @test norm(x.state) ≈ 1
    end
    @testset "Ferromagnetic Z" begin
        x = GetFerromagneticStateZ(1)
        @test x.L == 1
        @test x.state == [0.0+0.0im, 1.0+0.0im] || x.state == [1.0+0.0im, 0.0+0.0im] 
        @test size(x.state) == (2,)
        @test norm(x.state) ≈ 1
    end
    @testset "Ferromagnetic XX" begin
        x = GetFerromagneticStateX(2)
        twox = 2*x

        xplus  = ones(Complex{Float64}, 4)
        xminus = ones(Complex{Float64}, 4)
        xminus[2] = -1
        xminus[3] = -1

        @test x.L == 2
        @test norm(x.state) ≈ 1
        @test twox.state == xplus || twox.state == xminus 
    end
end