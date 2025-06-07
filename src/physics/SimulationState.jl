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

# Find the site number (position) in the MPS where index i appears
function find_site(i::Index, ψ::MPS)
    for site in 1:length(ψ)
        if i in inds(ψ[site])
            return site
        end
    end
    error("Index $i not found in MPS")
end

function apply_1site!(state::SimulationState, gate::ITensor)
    # get the site index for the gate
    s = siteind(gate)
    i = find_site(s, state.ψ)

    # create the gate
    T = gate * state.ψ[i]
    noprime!(T)

    # apply the gate in place
    state.ψ[i] = T
end

function apply_2site_adjacent!(state::SimulationState, gate::ITensor, parms::Parameters)
    # get the site indices for the gate
    s1, s2 = siteinds(gate)
    i1 = find_site(s1, state.ψ)
    i2 = find_site(s2, state.ψ)
    @assert i2 == i1 + 1 "Site indices must be adjacent"

    # orthogonalize to the left site
    orthogonalize!(state.ψ, i1)

    # apply the gate
    T = gate * (state.ψ[i1] * state.ψ[i2])
    noprime!(T)

    # get the left indices to pass to the SVD
    l = uniqueinds(state.ψ[i], state.ψ[i+1])

    # apply the SVD to the pair of sites
    U, S, V = ITensor.svd(T, l; maxdim=parms.maxdim, cutoff=parms.cutoff)

    # update the MPS
    state.ψ[i] = U
    state.ψ[i+1] = S * V
end