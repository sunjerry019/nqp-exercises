#!/usr/bin/env julia

include("qm.jl")

using .QuantumMechanics
using LinearAlgebra
using Test

# (2.a)
function H(L :: Integer, h :: Real) :: Operator
    _H = ZeroOp(L)

    for i in 1:L
        # Periodic boundary condition
        nn = (i+1) % L
        nn = nn == 0 ? L : nn

        _H -= IdentityOp(L)*Z(L, i)*Z(L, nn)
    end

    for i in 1:L
        _H -= h*X(L, i)
    end
    
    return _H
end

@testset "2.a) Hamiltonian Tests" verbose=true begin
    @testset "2 Sites, h = 1" begin
        x = H(2, 1)
    
        correct_zpart = -[
            2  0  0 0;
            0 -2  0 0;
            0  0 -2 0;
            0  0  0 2
        ]
        correct_xpart = -[
            0 1 1 0;
            1 0 0 1;
            1 0 0 1;
            0 1 1 0
        ]
    
        @test x.L == 2
        @test x.matrix == (correct_zpart + correct_xpart)
    end
    
    @testset "1 Site, h = 1" begin
        x = H(1,1)

        @test x.L == 1
        @test x.matrix == -[1 1; 1 1]
    end
end

# (2.b)