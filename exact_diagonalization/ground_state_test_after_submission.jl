#!/usr/bin/env julia

push!(LOAD_PATH, "./")

using QM
using TransverseFieldIsing

print(ground_state_energy(14, 10))
