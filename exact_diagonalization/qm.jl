#!/usr/bin/env julia

module QuantumMechanics
    using Random
    using LinearAlgebra

    import Base: * # To extend an operator, you must first import it

    export GetRandomState, GetFerromagneticStateZ, GetFerromagneticStateX, *

    abstract type HilbertSpace end

    mutable struct State <: HilbertSpace
        L     :: Integer  # Lattice Sites
        state :: Vector{Complex{Float64}}  
        # The actual state vector, quantization axis = z
        # We use the tensor product basis
    end

    function *(A :: Number, B :: State) 
        C = State(B.L, A * B.state)
        return C
    end

    function GetRandomState(L :: Integer) :: State
        vec = rand(Complex{Float64}, 2^L)
        N   = norm(vec)
        return State(L, vec/N)
    end

    """
    Returns a random Ferromagnetic state along z  

    Ferromagnatic state = all spins align.  
    In the tensor product basis, with z as the quantization basis,  
    this is represented with [1 0 .. 0] (spin up) or [0 .. 0 1] (spin down)
    """
    function GetFerromagneticStateZ(L :: Integer) :: State
        vec = zeros(Complex{Float64}, 2^L)
        vec[rand(Bool) ? 1 : end] = 1
        return State(L, vec)
    end

    """
    Returns a random Ferromagnetic state along x

    We use the kronecker (tensor) product this time.
    """
    function GetFerromagneticStateX(L :: Integer) :: State
        # Create a non-normalized eigenstate first
        vec = ones(Complex{Float64}, 2)
        vec[2] = rand(Bool) ? 1 : -1

        # Tensor product it L times
        args = ntuple(x -> vec, L)
        fullvec = kron(args...) # Splatting

        return State(L, fullvec/norm(fullvec))
    end
end