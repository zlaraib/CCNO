using DelimitedFiles
include("gates_function.jl")  # Include the gates_functions.jl file
include("chkpt_hdf5.jl") 
include("momentum.jl")
include("constants.jl")
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

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, chkptdir, checkpoint_every, do_recover, recover_type, recover_iteration, save_data::Bool, periodic=true)

    t_initial = 0.0
    iteration = 0
    t_recover = 0.0  # Variable to store the initial recovery time 

    if do_recover
        
        if recover_type == "auto"
            println("auto recovery: recovering from last iteration")
            recover_iteration = -1
            checkpoint_files = readdir(chkptdir)
            checkpoint_files = filter(f -> endswith(f, ".h5"), checkpoint_files)
            checkpoint_files = sort(checkpoint_files)

            if !isempty(checkpoint_files)
                latest_checkpoint_file = checkpoint_files[end]
                checkpoint_filename = joinpath(chkptdir, latest_checkpoint_file)
            else
                error("No checkpoint files found in directory: $chkptdir. Unable to recover.")
            end

            s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, t_initial, iteration = recover_checkpoint_hdf5(checkpoint_filename)
            # Increment t_initial by τ to ensure it starts from the next expected value
            t_initial += τ
            s = siteinds(ψ)

        elseif recover_type == "manual"
            println("Manual recovery from iteration ", recover_iteration)
            
            # Create the checkpoint filename for the specified iteration
            checkpoint_filename = joinpath(chkptdir, "checkpoint.chkpt.it" * lpad(recover_iteration, 6, "0") * ".h5")
            
            if isfile(checkpoint_filename)
                println("Recovering from checkpoint: $checkpoint_filename")
                # Recover data from the specified checkpoint
                s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, t_initial, iteration = recover_checkpoint_hdf5(checkpoint_filename)
                
                # Increment t_initial by τ to ensure it starts from the next expected value
                t_initial += τ
                s = siteinds(ψ)
            else
                error("Checkpoint file not found for iteration $recover_iteration in directory: $chkptdir.")
            end
        end
     
    end    

    # Create empty array s to... 
    Sz_array = Float64[] # to store sz values 
    Sy_array = Float64[] # to store sy values 
    Sx_array = Float64[] # to store sx values 
    t_array = [] # to store t values 
    prob_surv_array = Float64[]   # to store survival probability values 
    x_values = []  # to store x values for all sites 
    pₓ_values = [] # to store px vector values for all sites
    ρₑₑ_array = Float64[] # to store ρₑₑ values
    ρ_μμ_array = Float64[] # to store ρ_μμ values
    ρₑμ_array = Float64[] # to store ρₑμ values
    
    ρₑμ_at_t1 = nothing  # Initialize a variable to store ρₑμ at t1
    ρₑμ_at_t2 = nothing  # Initialize a variable to store ρₑμ at t2
    Δt = t2 - t1 #time difference between growth rates 

    # H = Hamiltonian_mpo(s, N, B, N_sites, Δx, Δm², p, x, Δp, shape_name,L, τ, energy_sign, periodic)
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
        half_N = div(N_sites, 2)  # Calculate half of the number of sites

        # compute expectation value of sy and sx using S+ and S- (inbuilt operator in ITensors library) at the first site on the chain
        if p == zeros(N_sites, 3) #for rogerro's case only (b/c S+ S- needed to keep conservation of QN number)
            sy_tot = -0.5 *im * (expect(complex(ψ), "S+") - expect(complex(ψ), "S-")) 
            sx_tot = 0.5 * (expect(ψ, "S+") + expect(ψ, "S-"))
        else 
            sy_tot = expect(complex(ψ), "Sy")
            sx_tot = expect(ψ, "Sx")
        end
        println("sz_tot= ",sz_tot )
        println("sy_tot= ",sy_tot )
        println("sx_tot= ",sx_tot )
        if shape_name!=="none" #specfic to the inhomogenous case Test4
            sz = mean(abs.(sz_tot[1:half_N]))  # Take the abs value fo all enteries till half_N and then take the mean of that first half of the sz_tot array
            sy = mean(abs.(sy_tot[1:half_N])) # Take the abs value fo all enteries till half_N and then take the mean of that first half of the sy_tot array
            sx = mean(abs.(sx_tot[1:half_N]))  # Take the abs value fo all enteries till half_N and then take the mean of that first half of the sx_tot array
        else 
            #for all inhomo cases:
            sz = sz_tot[1]
            sy = sy_tot[1]
            sx = sx_tot[1]
        end
        println("sz= ",sz )
        println("sy= ",sy )
        println("sx= ",sx)

        # # seperate loop of the forst site for all noninhomo cases 
        # # compute expectation value of Sz (inbuilt operator in ITensors library) at the first site on the chain
        # sz = expect(ψ, "Sz"; sites=1)
        # # compute expectation value of sy and sx using S+ and S- (inbuilt operator in ITensors library) at the first site on the chain
        # if p == zeros(N_sites, 3) #for rogerro's case only (b/c S+ S- needed to keep conservation of QN number)
        #     sy = -0.5 *im * (expect(complex(ψ), "S+"; sites=1) - expect(complex(ψ), "S-"; sites=1)) #re-check
        #     sx = 0.5 * (expect(ψ, "S+"; sites=1) + expect(ψ, "S-"; sites=1)) #recheck
        #     println("sz= ",sz )
        #     println("sy= ",sy )
        #     println("sx= ",sx)
        # else 
        #     sy = expect(complex(ψ), "Sy"; sites=1)
        #     sx = expect(ψ, "Sx"; sites=1)
        #     println("sz= ",sz )
        #     println("sy= ",sy )
        #     println("sx= ",sx)
        # end

        # add an element sz to ... 
        push!(Sz_array, sz) # .. the end of Sz array 
        push!(Sy_array, sy) # .. the end of Sy array 
        push!(Sx_array, sx) # .. the end of Sx array 
        
        # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
        prob_surv = 0.5 * (1 - 2 * sz)

        # add an element prob_surv to the end of  prob_surv_array 
        push!(prob_surv_array, prob_surv)

        if B[1] == 1
            println("$t $sz")
        else 
            println("$t $prob_surv")
        end
        
        # recall that in our code sigma_z = 2*Sz so make sure these expressions are consistent with "Sz in ITensors" 
        ρₑₑ = ( (2 * sz) + 1)/2 
        push!(ρₑₑ_array,abs(ρₑₑ))
        ρ_μμ = ( (-2 * sz) + 1)/2 
        push!(ρ_μμ_array,abs(ρ_μμ))
        ρₑμ = sqrt(sx^2 + sy^2)
        push!(ρₑμ_array, ρₑμ)

        # Check if the current time is approximately t1
        if abs(t - t1) < τ / 2
            ρₑμ_at_t1 = ρₑμ
        end
        
        # Check if the current time is approximately t2
        if abs(t - t2) < τ / 2
            ρₑμ_at_t2 = ρₑμ
        end
        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(gates, ψ; cutoff, maxdim)
        # ψ = tdvp(H, ψ,  -im *τ;   nsweeps=1,
        # reverse_step=true,outputlevel=1)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)

        if save_data
            isdir(chkptdir) || mkpath(chkptdir)
            println("iteration: ",  iteration, "  ", iteration % checkpoint_every)
            if iteration % checkpoint_every == 0 
                println("CREATE CHECKPOINT AT ITERATION = ", iteration, " TIME = ", t)
                checkpoint_filename = joinpath(chkptdir, "checkpoint.chkpt.it" * lpad(iteration, 6, "0") * ".h5")
                checkpoint_simulation_hdf5(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, gates, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
            end

            iteration = iteration + 1
        end
    end
    t_array = t_initial:τ:ttotal

    # After the time evolution loop, calculate and print the growth rate from the ratio if both values have been captured
    if ρₑμ_at_t1 !== nothing && ρₑμ_at_t2 !== nothing
        global Im_Ω = (1/Δt ) * log(ρₑμ_at_t2 / ρₑμ_at_t1)
        println("Growth rate of flavor coherence of ρₑμ at t2 to ρₑμ at t1: $Im_Ω")
    else
        global Im_Ω = 1
        println("ρₑμ was not captured at both t1 and t2.")
    end

    if save_data && !do_recover
        save_data = isdir(datadir) || mkpath(datadir)
        # Writing data to files with corresponding headers
        fname1 = joinpath(datadir, "t_<Sz>_<Sy>_<Sx>.dat")
        writedlm(fname1, [t_array Sz_array Sy_array Sx_array])
        
        fname2 = joinpath(datadir, "t_probsurv.dat")
        writedlm(fname2, [t_array prob_surv_array])
        
        fname3 = joinpath(datadir, "t_xsiteval.dat")
        writedlm(fname3, [t_array x_values])
        
        fname4 = joinpath(datadir, "t_pxsiteval.dat")
        writedlm(fname4, [t_array pₓ_values])
        
        fname5 = joinpath(datadir, "t_ρₑₑ.dat")
        writedlm(fname5, [t_array ρₑₑ_array])
        
        fname6 = joinpath(datadir, "t_ρ_μμ.dat")
        writedlm(fname6, [t_array ρ_μμ_array])
        
        fname7 = joinpath(datadir, "t_ρₑμ.dat")
        writedlm(fname7, [t_array ρₑμ_array])

        fname8 = joinpath(datadir, "Im_Ω.dat")
        writedlm(fname8, [Im_Ω])
    
    
    elseif save_data && do_recover

        # Determine the start index for appending data
        start_index = iteration >= checkpoint_every ? checkpoint_every + 1 : 0
    
        # Function to append data as new rows to the existing file after every `checkpoint_every` values
        function append_data(filename, new_data)
            open(filename, "a") do f
                writedlm(f, new_data)
            end
        end
    
        # Append new values for each of the saved files, starting from the appropriate index
        if start_index > 0
            fname1 = joinpath(datadir, "t_<Sz>_<Sy>_<Sx>.dat")
            append_data(fname1, [t_array[start_index:end] Sz_array[start_index:end] Sy_array[start_index:end] Sx_array[start_index:end]])
    
            fname2 = joinpath(datadir, "t_probsurv.dat")
            append_data(fname2, [t_array[start_index:end] prob_surv_array[start_index:end]])
    
            fname3 = joinpath(datadir, "t_xsiteval.dat")
            append_data(fname3, [t_array[start_index:end] x_values[start_index:end]])
    
            fname4 = joinpath(datadir, "t_pxsiteval.dat")
            append_data(fname4, [t_array[start_index:end] pₓ_values[start_index:end]])
    
            fname5 = joinpath(datadir, "t_ρₑₑ.dat")
            append_data(fname5, [t_array[start_index:end] ρₑₑ_array[start_index:end]])
    
            fname6 = joinpath(datadir, "t_ρ_μμ.dat")
            append_data(fname6, [t_array[start_index:end] ρ_μμ_array[start_index:end]])
    
            fname7 = joinpath(datadir, "t_ρₑμ.dat")
            append_data(fname7, [t_array[start_index:end] ρₑμ_array[start_index:end]])
    
            # Saving just a single value here, appended as a new row
            fname8 = joinpath(datadir, "Im_Ω.dat")
            open(fname8, "a") do f
                writedlm(f, [Im_Ω])
            end
        end
    end
        
    return Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, Im_Ω, t_recover
end


