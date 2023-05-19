#!/usr/bin/env julia

module QuantumMechanics
    using Random
    using LinearAlgebra

    import Base: * # To extend an operator, you must first import it

    # States
    export State
    export GetRandomState, GetFerromagneticStateZ, GetFerromagneticStateX, *
    # Operators
    export Operator, OperatorSingleSite
    export IdentityOp, ExpandToFullHilbertSpace
    export X, Y, Z

    abstract type HilbertSpace end

    # STATES
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

    # OPERATORS
    mutable struct Operator <: HilbertSpace
        L      :: Integer # Full Number of Lattice Sites
        matrix :: Matrix{Complex{Float64}}
    end
    mutable struct OperatorSingleSite <: HilbertSpace
        L      :: Integer # Full Number of Lattice Sites
        site   :: Integer # index of the lattice site its supposed to act on, 1 - indexed!
        matrix :: Matrix{Complex{Float64}}
    end

    function IdentityOp(L :: Integer) :: Operator
        return Operator(L, Diagonal(ones(2^L))) # Alternative to just using I
    end
    function X(L :: Integer, site :: Integer) :: OperatorSingleSite
        if (site > L)
            throw(DimensionMismatch(string("L larger than site, operator not possible")))
        end

        _x = [0 1; 1 0]
        return OperatorSingleSite(L, site, _x)
    end
    function Y(L :: Integer, site :: Integer) :: OperatorSingleSite
        if (site > L)
            throw(DimensionMismatch(string("L larger than site, operator not possible")))
        end
        _y = [0 -1im; 1im 0]
        return OperatorSingleSite(L, site, _y)
    end
    function Z(L :: Integer, site :: Integer) :: OperatorSingleSite
        if (site > L)
            throw(DimensionMismatch(string("L larger than site, operator not possible")))
        end
        _z = [1 0; 0 -1]
        return OperatorSingleSite(L, site, _z)
    end

    function *(A :: Operator, B :: Operator) :: Operator
        if (A.L != B.L)
            throw(DimensionMismatch(string("Dimensions do not match: ", A.L, " != " , B.L)))
        end

        return Operator(A.L, A.matrix * B.matrix)
    end
    
    function *(A :: Operator, B :: State) :: State
        if (A.L != B.L)
            throw(DimensionMismatch(string("Dimensions do not match: ", A.L, " != " , B.L)))
        end

        return State(A.L, A.matrix * B.state)
    end

    function *(A :: OperatorSingleSite, B :: OperatorSingleSite) :: Operator
        if (A.L != B.L)
            throw(DimensionMismatch(string("Dimensions do not match: ", A.L, " != " , B.L)))
        elseif (A.site != B.site)
            throw(DomainError("A and B act on different sites, unable to multiply"))
        end

        return OperatorSingleSite(A.L, A.site, A.matrix * B.matrix)
    end

    function ExpandToFullHilbertSpace(ssp :: OperatorSingleSite) :: Operator
        if (ssp.site > ssp.L)
            throw(DimensionMismatch(string("L larger than site, operator invalid")))
        end

        sites_before = ntuple(x -> Diagonal(ones(2)), ssp.site - 1)
        sites_after  = ntuple(x -> Diagonal(ones(2)), ssp.L - ssp.site)
        all_matrices = tuple(sites_before..., ssp.matrix, sites_after...)
        return Operator(ssp.L, kron(all_matrices...))
    end
end