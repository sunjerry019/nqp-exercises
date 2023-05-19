#!/usr/bin/env julia

include("qm.jl")
include("transversefieldising.jl")

using .QuantumMechanics
using .TransverseFieldIsing
using Test

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

