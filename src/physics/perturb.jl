using ITensors
using ITensorMPS

# This file generates the create_perturbation_gates function that holds ITensors Trotter gates and returns the dimensionless unitary 
# operators that will generate the perturbation via this hamiltonian which includes effects of the vacuum one-body potential for each site 
# Then, this file generates the evolve_perturbation function which utilizes the unitary operators created as perturb_gates from the 
# create_perturbation_gates function to evolve the initial ψ state in time and return the normalized perturbed state after evolution.

@doc """
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the gates function.                                                                   
    s = site index array (dimensionless and unitless)          
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    B = normalized vector that allows perurbation in different directions (dimensionless constant)
    N_sites = Total no.of sites (dimensionless and unitless)
    τ = time step (sec)
    energy_sign = array of sign of the energy (1 or -1): 1 for neutrinos and -1 for anti-neutrinos
    maxdim = max bond dimension in MPS truncation (unitless and dimensionless)
    cutoff = truncation threshold for the SVD in MPS representation (unitless and dimensionless)
"""
function create_perturbation_gates(params::CCNO.parameters, s::Vector{Index{Int64}}, k::Float64, B_pert::Vector{Float64})
    
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[] 

    for i in 1:(params.N_sites-1)
        for j in i+1:params.N_sites
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B_pert) == 1
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.

            # add perturbation via one-body oscillation term to the Hamiltonian
            hj= B_pert[1] * op("Sx", s[i]) * op("Id", s[j]) +
                B_pert[2] * op("Sy", s[i]) * op("Id", s[j]) +
                B_pert[3] * op("Sz", s[i]) * op("Id", s[j])

            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * params.α * hj) #Gj has factor of 1/hbar sinceonly Richers inhomo tests need perturbation

            # The push! function adds (appends) an element to the end of an array;
            # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
            # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
            push!(gates, Gj)
        end
    end

    # append! adds all the elements of a gates in reverse order (i.e. (N,N-1),(N-1,N-2),...) to the end of gates array.
    # appending reverse gates to create a second-order Trotter-Suzuki integration
    append!(gates, reverse(gates)) 
    return gates
end


function evolve_perturbation(params::CCNO.parameters, s,k,B_pert,ψ)

    # extract the gates array generated in the gates_function file
    perturb_gates = create_perturbation_gates(params, s,k, B_pert)

    ψ = apply(perturb_gates, ψ; params.cutoff, params.maxdim)
        
    # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
    # This is necessary to ensure that the MPS represents a valid quantum state.
    normalize!(ψ)
    return ψ
end


