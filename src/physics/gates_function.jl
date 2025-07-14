using ITensors
using ITensorMPS

@doc """
    Create a list of gates that need to be applied to the simulation state.
    params: CCNO.Parameters object defined in physics/Parameters.jl
    state: CCNO.SimulationSate object defined in physics/SimulationState.jl

    Return values are lists of dimensionless unitary operators that include vacuum and self-interactions
    gates_1site: list of gates applied to only one site, amenable to simple gate application method
    gates_2site_even: list of two-site gates applied to adjacent even pairs of sites. Amenable to simple two-site method.
    gates_2site_odd: list of two-site gates applied to adjacent odd pairs of sites. Amenable to simple two-site method.
    gates_2site_other: list of two site gates for sites that are not adjacent. Require the general "apply()" function.
"""
Base.@pure function create_gates(params::CCNO.Parameters, state::CCNO.SimulationState)
    
    # Define B vector that activates the vacuum oscillation term in Hamiltonian
    B = [sin(2*params.theta_nu), 0, -cos(2*params.theta_nu)] 
    B = B / norm(B) 

    # Create an empty ITensors array that will be our Trotter gates
    gates_1site = ITensor[] 
    gates_2site_even = ITensor[]
    gates_2site_odd = ITensor[]
    gates_2site_other = ITensor[]

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(state.p)  
    
    # define an array of vacuum oscillation frequencies (units of ergs)
    Δm²::Float64 = params.m2^2 - params.m1^2
    ω::Vector{Float64} = [Δm² / (2 * p_mod[i]) * state.energy_sign[i] for i in 1:params.N_sites]

    # Precompute operators for all sites
    Id::Vector{ITensor} = [op("Id", state.s[i]) for i in 1:length(state.s)]
    Sx::Vector{ITensor} = [op("Sx", state.s[i]) for i in 1:length(state.s)]
    Sy::Vector{ITensor} = [op("Sy", state.s[i]) for i in 1:length(state.s)]
    Sp::Vector{ITensor} = [op("S+", state.s[i]) for i in 1:length(state.s)]
    Sm::Vector{ITensor} = [op("S-", state.s[i]) for i in 1:length(state.s)]
    Sz::Vector{ITensor} = [op("Sz", state.s[i]) for i in 1:length(state.s)]
    
    # vacuum - one-site gate model. Faster, but order-dependent.
    for i in 1:params.N_sites
        if ω[i] != 0
            hj::ITensor = ω[i] * (B[1]*Sx[i] + B[2]*Sy[i] + B[3]*Sz[i])
            Gj::ITensor = exp(-im * params.τ * hj / hbar)
            push!(gates_1site, Gj)
        end
    end

    # loop over all pairs of sites and put gates into appropriate gate list
    for i in 1:(params.N_sites-1)
        for j in i+1:params.N_sites

            # get the shape function, multiplying all three directions together
            shape_result::Float64 = 1
            for d in 1:1
                shape_result *= shape_func(params, state, d, i, j)
            end
                
            geometric_factor::Float64 = geometric_func(params, p̂, i, j)
            interaction_strength::Float64 = 2.0 * √2 * G_F * (state.N[i] + state.N[j]) / (2*params.Δx^3) * shape_result * geometric_factor
            if interaction_strength != 0
                hj::ITensor = interaction_strength * (Sz[i]*Sz[j] + 0.5*Sp[i]*Sm[j] + 0.5*Sm[i]*Sp[j])

                # if neutrinos interacting with antineutrinos, H changes???
                #if state.energy_sign[i]*state.energy_sign[j] < 0
                #    hj *= -2

                # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors
                Gj::ITensor = exp(-im * params.τ/2 * hj / hbar)

                # The push! function adds (appends) an element to the end of an array;
                # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
                # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
                if i%2==0 && j==i+1
                    push!(gates_2site_even, Gj)
                elseif i%2==1 && j==i+1
                    push!(gates_2site_odd, Gj)
                else
                    push!(gates_2site_other, Gj)
                end
            end
        end 
    end

    return gates_1site, gates_2site_even, gates_2site_odd, gates_2site_other
end

