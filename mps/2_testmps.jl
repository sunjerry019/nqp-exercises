#!/usr/bin/env julia

push!(LOAD_PATH, "./")

include("testsetskip.jl")

using MatrixProductStates
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

        @test mp_state.m == (L == 1 ? 1 : m) # if L == 1, then m = 1 regardless
        @test mp_state.L == L
        @test (k == 1) ? size(mp_state.tensor_sets[k]) == (1,(L == 1 ? 1 : m),d) : ((k == L) ? size(mp_state.tensor_sets[k]) == (m,1,d) : size(mp_state.tensor_sets[k]) == (m,m,d))
        @test length(mp_state.tensor_sets) == L
    end

    @testset "Random State (Exact representation)" begin
        function check_size(ist_tensor_sets :: Vector{Array}, soll_sizes :: Vector{Tuple{Int64, Int64, Int64}})
            L = length(ist_tensor_sets)
            if L != length(soll_sizes)
                return false
            end
            
            for i in LinearIndices(soll_sizes)
                if size(ist_tensor_sets[i]) != soll_sizes[i]
                    return false
                end
            end
            return true
        end

        s = create_random_state(6, 2) # L, local_dof(k)
        @test s.L == 6
        @test s.m == 8
        @test check_size(s.tensor_sets, [(1,2,2),(2,4,2),(4,8,2),(8,4,2),(4,2,2),(2,1,2)])

        s = create_random_state(7, 2) # L, local_dof(k)
        @test s.L == 7
        @test s.m == 8
        @test check_size(s.tensor_sets, [(1,2,2),(2,4,2),(4,8,2),(8,8,2),(8,4,2),(4,2,2),(2,1,2)])

        s = create_random_state(2, 2) # L, local_dof(k)
        @test s.L == 2
        @test s.m == 2
        @test check_size(s.tensor_sets, [(1,2,2),(2,1,2)])

        s = create_random_state(1, 2) # L, local_dof(k)
        @test s.L == 1
        @test s.m == 1
        @test check_size(s.tensor_sets, [(1,1,2)])
    end
end

@testset "Collapse Singular Dimension" begin
    ns = create_seq_matrix_set(Int16, 2,3,4)
    s  = create_seq_matrix_set(Int16, 2,1,3)
    kk = create_seq_matrix_set(Int16, 2,3,1)

    @test_throws DimensionMismatch collapse_singular_dimension(ns)
    @test_throws DimensionMismatch collapse_singular_dimension(kk)

    ss = collapse_singular_dimension(s)
    @test size(ss) == (2,3)
end

@testset "Canonical Form" verbose=true begin
    @testset "Left Orthogonal" begin
        # create_random_state(L,d,m)
        m = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        d = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        L = abs(rand(Int, 1)[1] % 5) + 2 # min 2

        s = create_random_state(L, d, m)

        @test_throws DomainError make_orthogonal_left!(s, L)
    
        make_orthogonal_left!(s, 1)
        new = fuse_left(s.tensor_sets[1])

        @test (new' * new) ≈ I atol=10e-6
    end
    @testset "Right Orthogonal" begin
        # create_random_state(L,d,m)
        m = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        d = abs(rand(Int, 1)[1] % 5) + 1 # min 1
        L = abs(rand(Int, 1)[1] % 5) + 2 # min 2

        s = create_random_state(L, 3, 2)

        @test_throws DomainError make_orthogonal_right!(s, 1)
    
        make_orthogonal_right!(s, L)
        new = fuse_right(s.tensor_sets[L])

        @test (new * new') ≈ I atol=10e-6
    end
end
# Would also be good to check if the final matrices still multiply to give the same tensor
# Not done due to time-constraints, and not explicitly in the problem sheet