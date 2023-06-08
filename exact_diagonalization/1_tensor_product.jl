#!/usr/bin/env julia

push!(LOAD_PATH, "./")

using QM
using Test
using Random
using LinearAlgebra

@testset "State and Operator Generation" verbose=true begin
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
    @testset "Single Site Z > Full" begin
        matrix = [1 0; 0 -1]
        x = OperatorSingleSite(2, 2, matrix)
        y = ExpandToFullHilbertSpace(x)

        @test x.L == 2
        @test x.site == 2
        @test y.L == 2
        @test x.matrix == [1 0; 0 -1]
        @test y.matrix == [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
    end
    @testset "Single Site > Full Mismatch" begin
        matrix = [0 1; 1 0]
        x = OperatorSingleSite(1, 2, matrix)
        @test x.matrix == matrix

        try
            y = ExpandToFullHilbertSpace(x)
        catch e
            @test typeof(e) == DimensionMismatch
        end
    end
    @testset "Pauli-Matrices" begin
        x = X(2, 1); y = Y(2, 1); z = Z(2, 1)
        
        @test x.matrix == [0 1; 1 0]
        @test y.matrix == [0 -1im; 1im 0]
        @test z.matrix == [1 0; 0 -1]
    end
    @testset "Identity" begin
        q = IdentityOp(5) 

        @test q.matrix == I
    end
end

@testset "Operator on States" verbose=true begin
    @testset "Single Site spin flip" begin
        _state = GetFerromagneticStateZ(3)
        _spinflip_1 = X(3, 1)
        _newstate = _spinflip_1 * _state

        _correctstate = zeros(Complex{Float64}, 2^3)
        _correctstate[_state.state[1] == 1 ? 5 : 4] = 1
        _correctstate[_state.state[1] == 1 ? 1 : end] = 0

        # If spin up,   then after X, |100> = [0 0 0 0 1 0 0 0]
        # If spin down, then after X, |011> = [0 0 0 1 0 0 0 0]

        @test  _newstate.state == _correctstate
    end
end

@testset "Scalar Product on H" verbose=true begin
    @testset "Norm" begin
        x = GetRandomState(5)
        @test dot(x,x) ≈ norm(x.state)
    end
    @testset "Orthogonal States" begin
        x = GetFerromagneticStateZ(3, true)
        y = GetFerromagneticStateZ(3, false)

        @test dot(x,y) == 0
    end
end





