#!/usr/bin/env julia

module MatrixProductStates
    using LinearAlgebra

    export create_random_matrix_set, create_seq_matrix_set
    export fuse_left, fuse_right
    export split_left, split_right

    # We just use the Array objects from Julia
    # Ways to create multidimensional arrays
    # https://stackoverflow.com/questions/25561390/declare-and-initialise-3-dimensional-array
    function create_random_matrix_set(T :: Type, alpha :: Integer, beta :: Integer, d :: Integer) :: Array
        M_ab  = rand(T, alpha, beta)
        M_kab = cat(ntuple(x -> M_ab, d)..., dims = 3)
        return M_kab
    end

    function create_seq_matrix_set(T :: Type, alpha :: Integer, beta :: Integer, d :: Integer) :: Array
        M_kab = [convert(T, z*alpha*beta + y*alpha + x) for x in 0:alpha-1, y ∈ 0:beta-1, z = 0:d-1]
        return M_kab
    end

    # FUSING
    function fuse_left(A :: Array) :: Array
        # A[alpha, beta, k]
        alpha, beta, k = size(A)
        B = permutedims(A, [1,3,2])
        B = reshape(B, (alpha*k, beta))

        # B = [A[(kx%alpha)+1, n+1, (kx÷alpha)+1] for kx in 0:(alpha*k-1), n in 0:(beta-1)]
        return B
    end

    function fuse_right(A :: Array) :: Array
        # A[alpha, beta, k]
        alpha, beta, k = size(A)
        B = permutedims(A, [1,2,3])
        B = reshape(A, (alpha, beta*k))
        return B
    end

    # SPLITTING
    function split_left(A :: Array, d :: Integer) :: Array
        # A is M(m*d, n)

        # A[alpha, beta, k]
        md, n = size(A)
        if (md % d != 0) throw(DimensionMismatch(string("d invalid: ", md, " not divisible by " , d))) end

        m = md ÷ d
        B = reshape(A, (m, d, n))
        B = permutedims(B, [1,3,2])
        return B
    end

    function split_right(A :: Array, d :: Integer) :: Array
        # A is M(m, n*d)

        # A[alpha, beta, k]
        m, nd = size(A)
        if (nd % d != 0) throw(DimensionMismatch(string("d invalid: ", nd, " not divisible by " , d))) end

        n = nd ÷ d
        B = reshape(A, (m, n, d))
        return B
    end

    export MPS
    export construct_MPS

    # MPS "Class"
    mutable struct MPS
        L           :: Integer # number of lattice sites
        m           :: Integer # max bond dimensions
        ds          :: Vector{Integer} # the dimension of the local hilbert space
        tensor_sets :: Vector{Array}
    end

    function construct_MPS(tensor_sets :: Vector{Array}) :: MPS
        maxm = 0
        ds = Vector{Integer}()

        L = length(tensor_sets)
        for (i, tensor) in enumerate(tensor_sets)
            _m, _n, _k = size(tensor)

            # Check if the first and last
            if ((i == 1) && ((_m != 1) || (_n != size(tensor_sets[i+1])[1]))) || ((i == L) && ((_n != 1) || (_m != size(tensor_sets[i-1])[2])))
                throw(DimensionMismatch("The ends of the MPS must be such that the MPS contracts into a number!"))
            elseif (i != 1) && (i != L) && (_m != _n)
                throw(DimensionMismatch("Bulk should be described by square matrices!"))
            end

            push!(ds, _k)
            maxm = max(_m, maxm)
        end

        return MPS(L, maxm, ds, tensor_sets)
    end
end