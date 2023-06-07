#!/usr/bin/env julia

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
    for i in reverse((j+1):S.L)
        make_orthogonal_right!(S, i)
    end

    return S
end

@testset "Sweep" begin
    # create_random_state(L,d,m)
    m = 3
    d = 2
    S = create_random_state(7,d)

    j = 5

    new_S = sweep(S, j)

    q = fuse_left(new_S.tensor_sets[1])
    accumulator = (q' * q)
    # println(size(accumulator))
    # println(size(new_S.tensor_sets[2]))

    @test accumulator ≈ I atol=10e-6
    for i in 2:(j-1)
        # Absorb accumulator into the right matrix first
        _right = split_right(accumulator * fuse_right(new_S.tensor_sets[i]), d)

        # Then do the A'A
        _q = fuse_left(_right)
        accumulator = _q' *  _q

        @test accumulator ≈ I atol=10e-6
    end

    # Test Right Orthogonal
end
