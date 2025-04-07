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
