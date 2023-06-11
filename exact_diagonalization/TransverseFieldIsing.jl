#!/usr/bin/env julia

push!(LOAD_PATH, "./")

module TransverseFieldIsing
    using QM
    using LinearAlgebra
    using Test

    export H, ground_state_energy, ground_state_energy_density
    export ground_state, average_ground_state_magnetization

    # (2.a)
    function H(L :: Integer, h :: Real) :: Operator
        _H = ZeroOp(L)

        for i in 1:L
            # Periodic boundary condition
            nn = (i+1) % L
            nn = nn == 0 ? L : nn

            _H -= IdentityOp(L)*Z(L, i)*Z(L, nn)
        end

        for i in 1:L
            _H -= h*X(L, i)
        end
        
        return _H
    end

    # (2.b)
    function ground_state_energy(L :: Integer, h :: Real) :: Real
        _H = H(L, h)
        evals = eigvals(_H.matrix)
        return minimum(evals)
    end
    function ground_state_energy_density(L :: Integer, h :: Real) :: Real
        return ground_state_energy(L, h) / L
    end
    # print(ground_state_energy_density(2,1)) gives -1.4142135623730947

    function ground_state(L :: Integer, h :: Real) :: State
        _H = H(L, h)
        evecs = eigvecs(_H.matrix)
        evecs_array = [evecs[:,x] for x in axes(evecs,1)]
        # cite: https://discourse.julialang.org/t/how-do-i-create-vectors-from-matrix-columns/4139/2

        # display(evecs_array)
        # for _evec in evecs_array
        #     _S = State(2, _evec)
        #     println(expval(_H, _S))
        # end

        _S = State(L, evecs_array[1]) # the first one will have the lowest energy

        return normalize(_S)
    end

    # display(ground_state(2, 1).state)


    function average_ground_state_magnetization(L :: Integer, h :: Real) :: Real
        m = 0

        _gs = ground_state(L, h) 

        for j in 1:L
            _x_j = IdentityOp(L) * Z(L, j)
            m += expval(_x_j, _gs)
        end

        return m/(2*L)
    end

    # print(average_ground_state_magnetization(2, 1))
end