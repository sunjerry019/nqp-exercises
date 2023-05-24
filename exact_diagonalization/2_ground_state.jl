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

# _h_array = collect(range(0, 100, length=100))
# L = 3

# for h_idx in eachindex(_h_array)
#     println(ground_state_energy_density(L, _h_array[h_idx]))
# end

# (2.d)
function approxequal(A :: TransverseFieldIsing.QuantumMechanics.State, B :: State) :: Bool
    return (A.L == B.L) ? (isapprox(A.state, B.state, rtol = 0.0000001)) : false
end

@testset "Limit h → 0" begin
    L = 4
    gs_0 = ground_state(L, 0)

    @test approxequal(gs_0, GetFerromagneticStateZ(L, true)) || approxequal(gs_0, GetFerromagneticStateZ(L, false))
end

@testset "Limit h → ∞" begin
    L = 2
    gs_0 = -1 * ground_state(L, 10000000)

    @test approxequal(gs_0, GetFerromagneticStateX(L, true)) || approxequal(gs_0, GetFerromagneticStateX(L, false))
end


