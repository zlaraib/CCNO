using ITensorMPS
using ITensors

# Struct of state variables
Base.@kwdef mutable struct SimulationState
    xyz::Array{Float64,2}
    p::Array{Float64, 2}

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.
    s::Vector{Index}
    
    energy_sign::Vector{Int}
    ψ::MPS
    N::Vector{Float64}
end

function sort_sites!(state::SimulationState)
    # reorder the sites in the MPS according to x coordinate
    perm = sortperm(state.xyz[:,1])

    # Create an index map: current_pos[i] = position of site i in current MPS
    current_pos = collect(1:length(state.N))

    for new_pos in 1:length(state.N)
        # Find site that should go to position `new_pos`
        site_to_move = perm[new_pos]

        # Get its current position
        pos = current_pos[site_to_move]

        if pos == new_pos
            continue  # already in the right position
        end

        # move the site
        state.ψ = movesite(state.ψ, pos=>new_pos)

        # update current_pos. Every value between new_pos and pos increases by 1
        for i in 1:length(state.N)
            if current_pos[i] >= new_pos && current_pos[i] < pos
                current_pos[i] += 1
            end
        end
        current_pos[site_to_move] = new_pos
    end

    # also reorder the other contents of the state
    state.xyz = state.xyz[perm,:]
    state.p = state.p[perm,:]
    state.energy_sign = state.energy_sign[perm]
    state.N = state.N[perm]
    state.s = state.s[perm]
end
