#!/usr/bin/env julia

include("mps.jl")

using .MatrixProductStates
using LinearAlgebra

function sweep_two_states_LTR(ϕ :: MPS, ϕ_tilde :: MPS) :: Tuple{MPS, MPS, Dict{Integer, Union{Array, UniformScaling}}}
    if ϕ.L != ϕ_tilde.L
        throw(DimensionMismatch("The two states must have the same num of sites."))
    end

    ϕ       = construct_MPS(deepcopy(ϕ.tensor_sets))
    ϕ_tilde = construct_MPS(deepcopy(ϕ_tilde.tensor_sets))

    omega_L = Dict{Integer, Union{Array, UniformScaling}}()
    omega_L[0] = I # 0 is -1 in the question

    # Sweep from 1 to (L-1)
    for j in 1:(ϕ.L-1)
        make_orthogonal_left!(ϕ, j)
        make_orthogonal_left!(ϕ_tilde, j)

        _left = split_right(omega_L[j-1] * fuse_right(ϕ.tensor_sets[j]), ϕ.ds[j])
        omega_L[j] = fuse_left(ϕ_tilde.tensor_sets[j])' * fuse_left(_left)
    end

    return ϕ, ϕ_tilde, omega_L
end

ψ       = create_random_state(5,2)
ψ_tilde = create_random_state(5,2)
x, x_tilde, omega_L = sweep_two_states_LTR(ψ, ψ_tilde)

display(omega_L)

# Prepare a state with gauge center at j = 0
# run sweep_two_states_LTR