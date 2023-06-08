#!/usr/bin/env julia

include("mps.jl")
include("testsetskip.jl")

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

@testset "Sweep" verbose=true begin
    # create_random_state(L,d) // Exact Repr
    d = abs(rand(Int, 1)[1] % 5) + 2 # min 2
    L = abs(rand(Int, 1)[1] % 7) + 2 # min 2
    j = abs(rand(Int, 1)[1] % L) + 1 # min 1
    # println("$L sites, $d, Center $j")

    S = create_random_state(L,d)    
    new_S = sweep(S, j)

    @testset "Test Left Orthogonal" begin
        # Run only if there are left orthogonal sites
        left_orthos = length(1:(j-1))
        
        @test j == 1 ? left_orthos == 0 : left_orthos >= 1

        if left_orthos > 0
            q = fuse_left(new_S.tensor_sets[1])
            # First do the left most
            accumulator_L = (q' * q)
            @test accumulator_L ≈ I atol=10e-6
            for i in 2:(j-1)
                # Absorb accumulator into the right matrix first
                _right = split_right(accumulator_L * fuse_right(new_S.tensor_sets[i]), d)

                # Then do the A'A
                _q = fuse_left(_right)
                accumulator_L = _q' *  _q

                @test accumulator_L ≈ I atol=10e-6
            end
        end
    end
    
    @testset "Test Right Orthogonal" begin
        # Run only if there are left orthogonal sites
        right_orthos = length((j+1):L)
        
        @test j == L ? right_orthos == 0 : right_orthos >= 1

        if right_orthos > 0
            q = fuse_right(new_S.tensor_sets[end])
            println(size(q))

            # First do the right most
            accumulator_R = (q * q')
            @test accumulator_R ≈ I atol=10e-6
            for i in reverse((j+1):length(new_S.tensor_sets)-1)
                # Absorb accumulator into the left matrix first
                _left = split_left(fuse_left(new_S.tensor_sets[i]) * accumulator_R, d)

                # Then do the AA'
                _q = fuse_right(_left)
                accumulator_R = _q *  _q'

                @test accumulator_R ≈ I atol=10e-6
            end
        end
    end
end




