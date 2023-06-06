#!/usr/bin/env julia

"""
[Question] How do I contract the Matrices from the left and test that it all accumulates to Identity?
I can't seem to figure it out :(

I even tried to make a fuse middle, that fuses alpha and beta, but that doesn't seem to help...
"""

include("mps.jl")

using .MatrixProductStates
using LinearAlgebra
using Test

function sweep(S :: MPS, j :: Integer) :: MPS
    new_S = construct_MPS(deepcopy(S.tensor_sets))
    sweep!(new_S, j)
    return new_S
end

function sweep!(S :: MPS, j :: Integer) :: MPS
    # in-place sweep function
    
    for i in 1:(j-1)
        make_orthogonal_left!(S, i)
    end
    for i in (j+1):S.L
        make_orthogonal_right!(S, i)
    end

    return S
end

@testset "Sweep" begin
    # create_random_state(L,d,m)
    m = 3
    S = create_random_state(7,2,m)

    j = 5

    new_S = sweep(S, j)

    q = fuse_left(new_S.tensor_sets[1])
    accumulator = (q' * q)
    println(size(accumulator))
    println(size(new_S.tensor_sets[2]))

    @test accumulator ≈ I atol=10e-6
    # for i in 2:(j-1)
    #     _q = fuse_left(new_S.tensor_sets[i])
    #     accumulator = accumulator * 
    # end

    # _alpha, _beta, _k = size(new_S.tensor_sets[2])
    # _q  = fuse_middle(new_S.tensor_sets[2])
    # # Contract over k
    # _qq = split_middle((_q * _q'), _alpha)
    # println(size(_qq))
    # # accumulator *= (_q' * _q) 
    # # display(accumulator)
end
