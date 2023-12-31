#!/usr/bin/env julia

include("testsetskip.jl")

push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../exact_diagonalization/")

using QM
using TransverseFieldIsing
using MatrixProductStates
using LinearAlgebra

S = QM.GetFerromagneticStateX(7)
Q = create_MPS_from_State(S)
display(Q)

