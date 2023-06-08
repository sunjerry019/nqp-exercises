#!/usr/bin/env julia

push!(LOAD_PATH, "../exact_diagonalization/")
push!(DEPOT_PATH, "../.DEPOT/")
# https://discourse.julialang.org/t/error-with-modules-sharing-types/37119/2

module MatrixProductStates
    using LinearAlgebra

    using TransverseFieldIsing
    using QM

    export create_random_matrix_set, create_seq_matrix_set
    export fuse_left, fuse_right
    export split_left, split_right
    export collapse_singular_dimension

    # INDEXING: A[alpha, beta, k]

    # We just use the Array objects from Julia
    # Ways to create multidimensional arrays
    # https://stackoverflow.com/questions/25561390/declare-and-initialise-3-dimensional-array
    function create_random_matrix_set(T :: Type, alpha :: Integer, beta :: Integer, d :: Integer) :: Array
        M_ab  = rand(T, alpha, beta)
        # M_kab = cat(ntuple(x -> M_ab, d)..., dims = 3)
        M_kab = cat([rand(T, alpha, beta) for x in 1:d]..., dims = 3)
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

    # COLLAPSE DIM 1
    function collapse_singular_dimension(A :: Array) :: Array
        # Collapse the first dimension that is equal to 1
        # This is guaranteed to be either alpha or beta, since k is the last element
        _size = collect(size(A))

        # we only look at alpha or beta
        idx = findfirst(item -> item == 1, collect(_size[1:(end-1)]))
        if idx === nothing
            throw(DimensionMismatch("At least one dimension (except k) should be 1"))
        end

        deleteat!(_size, idx)

        return reshape(A, _size...)
    end

    export MPS
    export construct_MPS, create_random_state, create_MPS_from_State

    # MPS "Class"
    mutable struct MPS
        L           :: Integer # number of lattice sites
        m           :: Integer # max bond dimensions
        ds          :: Vector{<:Integer} # the dimension of the local hilbert space
        tensor_sets :: Vector{Array} # each tensor indexed with A[alpha, beta, k]
    end

    function construct_MPS(tensor_sets :: Vector{Array}) :: MPS
        maxm = 0
        ds = Vector{Integer}()

        L = length(tensor_sets)
        for (i, tensor) in enumerate(tensor_sets)
            _m, _n, _k = size(tensor)

            # Check if the first and last
            if ((i == 1) && (_m != 1)) || (i == L) && ((_n != 1))
                throw(DimensionMismatch("The ends of the MPS must be such that the MPS contracts into a number!"))
            elseif ((i < L) && (_n != size(tensor_sets[i+1])[1])) || ((i > 1) && (_m != size(tensor_sets[i-1])[2]))
                throw(DimensionMismatch("Neighbouring dimensions should match"))
            end

            push!(ds, _k)
            maxm = max(_m, maxm)
        end

        return MPS(L, maxm, ds, tensor_sets)
    end

    function create_random_state(sites :: Integer, local_dof :: Integer) :: MPS
        # For the first site, m must be k to encode all the information (exact representation)
        # Understanding: we start from a tensor with L legs, then try to split it
        # (1..j-1) = (m) - X - Y - (n) - (j..L)
        # Middle dimension is min(m,n) (think SVD)
        # Then the dimensions is 1xk, k*k2, k2xk3 .... then goes back down to 1

        # For use with actual states

        half = (sites - 1) ÷ 2
        is_even = (sites % 2 == 0)
        middle = sites ÷ 2

        ascending_arr = [ local_dof^x for x ∈ 1:half ]
        middle_element = is_even ? [ local_dof^middle ] : []
        m_arr = vcat(ascending_arr, middle_element, reverse(ascending_arr), [1])
        m_arr = Int64.(m_arr) # broadcast to Int64 (https://stackoverflow.com/a/48818789/3211506)

        return create_random_state(sites, local_dof, m_arr)
    end

    function create_random_state(sites :: Integer, local_dof :: Integer, m :: Integer) :: MPS
        # create_random_state(L,d,m)

        # This function simplifies by using square matrices for the bulk
        arr = Vector{Array}()
        if sites == 1
            # m is ignored
            push!(arr, create_random_matrix_set(Complex{Float64}, 1, 1, local_dof))
        elseif sites >= 2
            push!(arr, create_random_matrix_set(Complex{Float64}, 1, m, local_dof))
            if sites > 2
                for i in 1:(sites - 2)
                    # We simplify by creating square matrices
                    push!(arr, create_random_matrix_set(Complex{Float64}, m, m, local_dof))
                end
            end
            push!(arr, create_random_matrix_set(Complex{Float64}, m, 1, local_dof))
        end
        
        return construct_MPS(arr)
    end

    function create_random_state(sites :: Integer, local_dof :: Integer, m :: Array{<:Integer}) :: MPS
        # create_random_state(L,d,m[])

        # We assume the alpha_(-1) i.e. the first alpha is zero 
        # NOTE: we need (L+1) alphas to describe the matrices, not L as described in the question

        # This function creates from an array of alpha_js
        if length(m) != sites
            throw(DimensionMismatch("length of alpha_js must be the same as number of sites"))
        elseif m[sites] != 1
            throw(DimensionMismatch("alpha_L must be 1 so that the MPS contracts into a number"))
        end

        arr = Vector{Array}()
        erste = true
        for i in LinearIndices(m)
            if erste
                erste = false
                push!(arr, create_seq_matrix_set(Complex{Float64}, 1, m[i], local_dof))
            else
                push!(arr, create_seq_matrix_set(Complex{Float64}, m[i-1], m[i], local_dof))
            end
        end

        return construct_MPS(arr)
    end

    function create_MPS_from_State(state :: QM.State) :: Nothing
        state_tensor = reshape(state.state, state.ds)

        println(size(state_tensor))
        display(state_tensor)
    end

    export make_orthogonal_left!, make_orthogonal_right!
    function make_orthogonal_left!(state :: MPS, site :: Integer) :: MPS
        # in-place function

        if (site < 1 || site > (state.L - 1))
            throw(DomainError("LeftOrth: Site $site not in range 1:($(state.L)-1)"))
        end

        _tensor       = state.tensor_sets[site]
        _tensor_right = state.tensor_sets[site + 1]

        _fused_matrix = fuse_left(_tensor)
        _fact = svd(_fused_matrix)
        _U, _Svals, _Vt = _fact.U, _fact.S, _fact.Vt

        _new_ML = split_left(_U, state.ds[site])
        _new_MR = split_right(Diagonal(_Svals) * _Vt * fuse_right(_tensor_right), state.ds[site + 1])

        state.tensor_sets[site]   = _new_ML
        state.tensor_sets[site+1] = _new_MR

        return state
    end

    function make_orthogonal_right!(state :: MPS, site :: Integer) :: MPS
        # in-place function

        if (site < 2 || site > state.L)
            throw(DomainError("RightOrth: Site $site not in range 2:$(state.L)"))
        end

        _tensor      = state.tensor_sets[site]
        _tensor_left = state.tensor_sets[site-1]

        _fused_matrix = fuse_right(_tensor)
        _fact = svd(_fused_matrix)
        _U, _Svals, _Vt = _fact.U, _fact.S, _fact.Vt

        _new_MR = split_right(_Vt, state.ds[site])
        _new_ML = split_left(fuse_left(_tensor_left) * _U * Diagonal(_Svals), state.ds[site - 1])

        state.tensor_sets[site-1] = _new_ML
        state.tensor_sets[site]   = _new_MR

        return state
    end
end