#!/usr/bin/env julia

include("qm.jl")
include("transversefieldising.jl")

using .QuantumMechanics
using .TransverseFieldIsing

print(ground_state_energy(14, 10))
