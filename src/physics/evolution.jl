using DelimitedFiles
using ITensors
using ITensorMPS

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 

"""
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the evolve function.                                                                   
    s = site index array (dimensionless and unitless)          
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    B = array of normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
    N_sites = Total no.of sites (dimensionless and unitless)
    Δx = length of the box of interacting neutrinos at a site (cm)
    Δm² = difference in mass squared (erg^2)
    p = array of momentum vectors (erg)
    x = array of positions of sites (cm)
    Δp = width of shape function (cm)
    shape_name = name of the shape function (string) ["none","triangular","flat_top"]
    τ = time step (sec)
    energy_sign = array of sign of the energy (1 or -1): 1 for neutrinos and -1 for anti-neutrinos
    maxdim = max bond dimension in MPS truncation (unitless and dimensionless)
    cutoff = truncation threshold for the SVD in MPS representation (unitless and dimensionless)
    periodic = boolean indicating whether boundary conditions should be periodic
"""
function evolve(params::CCNO.parameters, s::Vector{Index{Int64}}, N::Vector{Float64}, B::Vector{Float64}, L::Float64, Δx::Float64, Δm²::Float64, p::Array{Float64,2}, x::Vector{Float64}, ψ::MPS, energy_sign::Vector{Int})

    t_initial = 0.0
    iteration = 0
    t_recover = t_initial # Variable to store the initial recovery time 

    if params.do_recover
        
        println("Manual recovery from file ", params.recover_file)
        
        if isfile(params.recover_file)
            println("Recovering from checkpoint: $recover_file")
            # Increment t_initial by τ to ensure it starts from the next expected value
            t_initial += params.τ
            # Recover data from the specified checkpoint
            s, τ, N, B, L, Δx, p, x, ψ, energy_sign, t_initial, iteration = recover_checkpoint_hdf5(params.recover_file)
            
            s = siteinds(ψ)
        else
            error("Checkpoint file not found")
        end
     
    end    

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(p,params.N_sites) 
    p̂ₓ= [sub_array[1] for sub_array in p̂]
    
    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
    for t in t_initial:params.τ:params.ttotal
        # extract the gates array generated in the gates_function file
        gates = create_gates(params, s, ψ,N, B, Δx, Δm², p, x, L, energy_sign)

        for i in 1:params.N_sites
            x[i] += p̂ₓ[i] * c * params.τ
            if params.periodic
                # wrap around position from 0 to domain size L
                x[i] = mod(x[i],L)

                # Checking if the updated x[i] satisfies the boundary conditions
                @assert (x[i] >= 0 && x[i] <= L)
            end
        end

        println("iteration= $iteration time= $t")
        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ params.ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(gates, ψ; params.cutoff, params.maxdim)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)

        if params.save_data

            store_data(params.datadir, t, ψ, x, p)

            mkpath(params.chkptdir)
            if iteration % params.checkpoint_every == 0 
                checkpoint_filename = joinpath(params.chkptdir, "checkpoint.chkpt.it" * lpad(iteration, 6, "0") * ".h5")
                checkpoint_simulation_hdf5(params, checkpoint_filename, s, N, B, L, Δx, Δm², p, x, ψ, energy_sign, t, iteration)
            end

            iteration = iteration + 1
        end
    end
    t_array = t_initial:params.τ:params.ttotal

end


