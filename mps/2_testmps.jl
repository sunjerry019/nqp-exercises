include("mps.jl")

using .MatrixProductStates
using Test
using LinearAlgebra

@testset "Creation of MPS" verbose=true begin
    @testset "Reject invalid dimensions" begin
        arr = Vector{Array}()
        for i in 1:3
            push!(arr, create_seq_matrix_set(Complex{Float64}, 3,3,2))
        end
        @test_throws DimensionMismatch mp_state = construct_MPS(arr)

        arr2 = Vector{Array}()
        push!(arr2, create_seq_matrix_set(Complex{Float64}, 1,3,2))
        for i in 1:3
            push!(arr2, create_seq_matrix_set(Complex{Float64}, 3,4,2))
        end
        push!(arr2, create_seq_matrix_set(Complex{Float64}, 4,1,2))
        @test_throws DimensionMismatch mp_state = construct_MPS(arr2)
    end

    @testset "Properties" begin
        m = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        d = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        L = abs(rand(Int, 1)[1] % 5) + 2 # min 2

        arr = Vector{Array}()

        push!(arr, create_seq_matrix_set(Complex{Float64}, 1,m,d))
        if L > 2
            for i in 1:(L-2)
                # We simplify by creating square matrices
                push!(arr, create_seq_matrix_set(Complex{Float64}, m,m,d))
            end
        end
        push!(arr, create_seq_matrix_set(Complex{Float64}, m,1,d))

        mp_state = construct_MPS(arr)

        @test mp_state.m == m
        @test mp_state.L == L
        @test all(mp_state.ds .== d)
    end

    @testset "Random State" begin
        m = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        d = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        L = abs(rand(Int, 1)[1] % 5) + 1 # min 1

        mp_state = create_random_state(L, d, m)

        k = abs(rand(Int, 1)[1] % L) + 1

        @test mp_state.m == m
        @test mp_state.L == L
        @test (k == 1) ? size(mp_state.tensor_sets[k]) == (1,m,d) : ((k == L) ? size(mp_state.tensor_sets[k]) == (m,1,d) : size(mp_state.tensor_sets[k]) == (m,m,d))
        @test length(mp_state.tensor_sets) == L
    end
end