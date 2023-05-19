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
function ground_state_energy(L :: Integer, h :: Real) :: Real
    _H = H(L, h)
    evals = eigvals(_H.matrix)
    return minimum(evals)
end
function ground_state_energy_density(L :: Integer, h :: Real) :: Real
    return ground_state_energy(L, h) / L
end
# print(ground_state_energy_density(2,1)) gives -1.4142135623730947

function ground_state(L :: Integer, h :: Real) :: State
    _H = H(L, h)
    evecs = eigvecs(_H.matrix)
    evecs_array = [evecs[:,x] for x in axes(evecs,1)]
    # cite: https://discourse.julialang.org/t/how-do-i-create-vectors-from-matrix-columns/4139/2

    # display(evecs_array)
    # for _evec in evecs_array
    #     _S = State(2, _evec)
    #     println(expval(_H, _S))
    # end

    _S = State(L, evecs_array[1]) # the first one will have the lowest energy

    return normalize(_S)
end

# display(ground_state(2, 1).state)


function average_ground_state_magnetization(L :: Integer, h :: Real) :: Real
    m = 0

    _gs = ground_state(L, h) 

    for j in 1:L
        _x_j = IdentityOp(L) * Z(L, j)
        m += expval(_x_j, _gs)
    end

    return m/(2*L)
end

print(average_ground_state_magnetization(2, 1))