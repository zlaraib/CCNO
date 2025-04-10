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
function create_perturbation_gates(s, k, B_pert, α, x, L, N_sites, energy_sign, τ)
    
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[] 

    # define an array of oscillation frequencies (units of ergs) of perturbation
    ω_pert = asin.([α * sin(k*x[i])* energy_sign[i] for i in 1:N_sites])
    println("perturb_ω = ", ω_pert)

    for i in 1:(N_sites-1)
        for j in i+1:N_sites
            #s_i, s_j are non-neighbouring spin site/indices from the s array
            s_i = s[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
            s_j = s[j]
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B_pert) == 1
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.

            # add perturbation via one-body oscillation term to the Hamiltonian
            hj= (
                (ω_pert[i] * B_pert[1] * op("Sx", s_i)* op("Id", s_j))  + (ω_pert[i] * B_pert[2] * op("Sy", s_i)* op("Id", s_j))  + (ω_pert[i] * B_pert[3] * op("Sz", s_i)* op("Id", s_j)) )

            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * τ/2 * hj*  1/hbar) #Gj has factor of 1/hbar sinceonly Richers inhomo tests need perturbation

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


function evolve_perturbation(s,k, τ_pert, B_pert, α, x, L, N_sites, ψ, cutoff, maxdim, energy_sign,ttotal)

    # extract the gates array generated in the gates_function file
    perturb_gates = create_perturbation_gates(s,k, B_pert, α, x, L, N_sites, energy_sign, τ_pert)

     # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
     for t in 0.0:τ_pert:ttotal  #perhaps not perturn till the end? stop somewhere in the mid 

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to τ_pert, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ τ_pert && break
        # apply each gate in perturb_gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(perturb_gates, ψ; cutoff, maxdim)
        

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end

    return ψ
end


