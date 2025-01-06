using DelimitedFiles
include("gates_function.jl")  # Include the gates_functions.jl file
include("chkpt_hdf5.jl") 
include("momentum.jl")
include("constants.jl")
include("../Utilities/save_datafiles.jl")
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

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, chkptdir, checkpoint_every, do_recover, recover_file, save_data::Bool, periodic=true)

    t_initial = 0.0
    iteration = 0
    t_recover = 0.0  # Variable to store the initial recovery time 

    if do_recover
        
        println("Manual recovery from file ", recover_file)
        
        if isfile(recover_file)
            println("Recovering from checkpoint: $recover_file")
            # Recover data from the specified checkpoint
            s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, t_initial, iteration = recover_checkpoint_hdf5(recover_file)
            
            # Increment t_initial by τ to ensure it starts from the next expected value
            t_initial += τ
            s = siteinds(ψ)
        else
            error("Checkpoint file not found")
        end
     
    end    

    # Create empty array s to... 
    Sz_array = [] # to store sz values for all sites  
    Sy_array = [] # to store sy values for all sites 
    Sx_array = [] # to store sx values for all sites 
    t_array = [] # to store t values (same for all sites)
    prob_surv_array = []   # to store survival probability values for all sites 
    x_values = []  # to store x values for all sites 
    pₓ_values = [] # to store px vector values for all sites
    ρₑₑ_array = [] # to store ρₑₑ values for all sites 
    ρ_μμ_array = [] # to store ρ_μμ values for all sites 
    ρₑμ_array = [] # to store ρₑμ values for all sites 
    
    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(p,N_sites) 
    p̂ₓ= [sub_array[1] for sub_array in p̂]
    
    # Main evolution loop
    is_first_iteration = true  # Flag to mark the first iteration of the loop
    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
    for t in t_initial:τ:ttotal
        # Save `t_recover` only once during the first iteration of the loop
        if is_first_iteration
            t_recover = t
            is_first_iteration = false
        end
        # extract the gates array generated in the gates_function file
        gates = create_gates(s, ψ,N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name,L, τ, energy_sign, periodic)
        push!(x_values, copy(x))  # Record x values at each time step
        px = p[:, 1]  # Extracting the first column (which corresponds to px values)
        push!(pₓ_values, copy(px)) # Record px values at each time step

        for i in 1:N_sites
            x[i] += p̂ₓ[i] * c * τ
            if periodic
                # wrap around position from 0 to domain size L
                x[i] = mod(x[i],L)

                # Checking if the updated x[i] satisfies the boundary conditions
                @assert (x[i] >= 0 && x[i] <= L)
            end
        end

        # compute the avg expectation value of Sz at all sites
        sz_tot = expect(ψ, "Sz")  # Compute Sz for each site and store the values in sz_tot

        # compute expectation value of sy and sx (inbuilt operator in ITensors library) at all sites on the chain
        sy_tot = expect(complex(ψ), "Sy")
        sx_tot = expect(ψ, "Sx")

        push!(Sz_array, sz_tot)  # Add all elements of sz_tot to Sz_array
        push!(Sy_array, sy_tot)  # Add all elements of sz_tot to Sz_array
        push!(Sx_array, sx_tot)  # Add all elements of sz_tot to Sz_array
        
        #survival probability for all sites (neutrino) to be found in its initial flavor state
        prob_surv_tot = 0.5 * (1 .- 2 .* sz_tot)
        push!(prob_surv_array, prob_surv_tot)

        # recall that in our code sigma_z = 2*Sz so make sure these expressions are consistent with "Sz in ITensors" 
        ρₑₑ_tot = ((2 .* sz_tot) .+ 1) ./ 2
        push!(ρₑₑ_array, abs.(ρₑₑ_tot))
        
        ρ_μμ_tot = ((-2 .* sz_tot) .+ 1) ./ 2
        push!(ρ_μμ_array, abs.(ρ_μμ_tot))
        
        ρₑμ_tot = sqrt.(sx_tot.^2 .+ sy_tot.^2)
        push!(ρₑμ_array, ρₑμ_tot)

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(gates, ψ; cutoff, maxdim)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)

        println("$iteration $t")
        if save_data

            store_data(datadir, t+τ, sz_tot, sy_tot, sx_tot, prob_surv_tot,x, px, ρₑₑ_tot,ρ_μμ_tot, ρₑμ_tot)

            mkpath(chkptdir)
            if iteration % checkpoint_every == 0 
                checkpoint_filename = joinpath(chkptdir, "checkpoint.chkpt.it" * lpad(iteration, 6, "0") * ".h5")
                checkpoint_simulation_hdf5(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, gates, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t+τ, iteration)
            end

            iteration = iteration + 1
        end
    end
    t_array = t_initial:τ:ttotal



    return Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover
end


