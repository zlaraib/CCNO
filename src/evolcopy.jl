using DelimitedFiles
include("gates_function.jl")  # Include the gates_functions.jl file
#include("H_MPO.jl") 
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


# function checkpoint_simulation(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
#     io = open(checkpoint_filename, "w")
#     # io = open("hello.txt", "w")
#     # println(io, s)
#     # println(io,τ)
#     # println(io,N )
#     # println(io,B)
#     # println(io,L)
#     # println(io,N_sites)
#     # println(io,Δx )
#     # println(io,Δm² )
#     # println(io,p )
#     # println(io,x )
    
#     # println(io,Δp )
#     # println(io,theta_nu)
#     # println(io,ψ)
#     # println(io,shape_name)
#     # println(io,energy_sign)
#     # println(io,cutoff)
#     # println(io,maxdim)
#     # println(io,t1)
#     # println(io,t2)
#     # println(io, ttotal)
#     println(io, t)
#     println(io, iteration)
#     println("checkpoint t: ", t)
#     close(io)

#     # Save the MPS ψ separately as a binary file
#     # write("$(checkpoint_filename)_psi.itensor", ψ)
# end

function checkpoint_simulation_hdf5(checkpoint_filename, ψ, t, iteration)
    # Open HDF5 file to write checkpoint
    f = h5open(checkpoint_filename, "w")
    
    # Save ψ (MPS) in ITensor format
    write(f, "psi", ψ)
    
    # Save other checkpoint metadata like `t` and `iteration`
    write(f, "t", t)
    write(f, "iteration", iteration)
    
    # Close the file after writing
    close(f)
    
    println("Checkpoint created at iteration $iteration, time $t")
end

function recover_checkpoint_hdf5(checkpoint_filename)
    
    # Open HDF5 file to read checkpoint
    f = h5open(checkpoint_filename, "r")
    
    # Read ψ (MPS) from ITensor format
    ψ = read(f, "psi", MPS)
    
    # Read other checkpoint metadata like `t` and `iteration`
    t_initial = read(f, "t")
    iteration = read(f, "iteration")
    
    # Close the file after reading
    close(f)
    
    println("Recovered from checkpoint at iteration $iteration, time $t_initial")
    
    return ψ, t_initial, iteration
end


# function checkpoint_simulation_hdf5(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
#     # Open an HDF5 file for writing (or create it if it doesn't exist)
#     f = h5open(checkpoint_filename, "w")
    
#     # Write scalar values
#     write(f, "τ", τ)
#     write(f, "L", L)
#     write(f, "N_sites", N_sites)
#     write(f, "Δx", Δx)
#     write(f, "Δm²", Δm²)
#     write(f, "Δp", Δp)
#     write(f, "theta_nu", theta_nu)
#     write(f, "cutoff", cutoff)
#     write(f, "maxdim", maxdim)
#     write(f, "t1", t1)
#     write(f, "t2", t2)
#     write(f, "ttotal", ttotal)
#     write(f, "t", t)
#     write(f, "iteration", iteration)
    
#     # Write arrays and vectors
#     write(f, "s", s)
#     write(f, "N", N)
#     write(f, "B", B)
#     write(f, "p", p)
#     write(f, "x", x)
#     write(f, "energy_sign", energy_sign)
    
#     # Write string values
#     write(f, "shape_name", shape_name)

#     # Write MPS ψ (ITensor type)
#     write(f, "ψ", ψ)
    
#     # Close the HDF5 file after writing
#     close(f)
    
#     println("Checkpoint created at $checkpoint_filename, iteration $iteration, time $t")
# end

# function recover_checkpoint_hdf5(checkpoint_filename)
#     # Open the HDF5 file for reading
#     f = h5open(checkpoint_filename, "r")
    
#     # Read scalar values
#     τ = read(f, "τ")
#     L = read(f, "L")
#     N_sites = read(f, "N_sites")
#     Δx = read(f, "Δx")
#     Δm² = read(f, "Δm²")
#     Δp = read(f, "Δp")
#     theta_nu = read(f, "theta_nu")
#     cutoff = read(f, "cutoff")
#     maxdim = read(f, "maxdim")
#     t1 = read(f, "t1")
#     t2 = read(f, "t2")
#     ttotal = read(f, "ttotal")
#     # t = read(f, "t")
#     t_initial = read(f, "t")
#     iteration = read(f, "iteration")
    
#     # Read arrays and vectors
#     s = read(f, "s")
#     N = read(f, "N")
#     B = read(f, "B")
#     p = read(f, "p")
#     x = read(f, "x")
#     energy_sign = read(f, "energy_sign")
    
#     # Read string values
#     shape_name = read(f, "shape_name")
    
#     # Read MPS ψ (ITensor type)
#     ψ = read(f, "ψ", MPS)
    
#     # Close the file after reading
#     close(f)
    
#     println("Recovered checkpoint from $checkpoint_filename, iteration $iteration, time $t_initial")

#     return s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t_initial, iteration
# end

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, chkptdir, checkpoint_every, do_recover, recover_type, recover_iteration, save_data::Bool, periodic=true)

    t_initial = 0.0
    iteration = 0

    # if do_recover
    #     if recover_type == "auto"
    #         println("auto recovery")
    #         f = open("/home/zohalaraib/Oscillatrino/tests/checkpoints/checkpoint.chkpt.it000004.txt", "r")
            

    #         t_initial = 0.2
    #         iteration = 2

    #         in_vals = readlines(f)
    #         println("in_vals", in_vals)
    #         t_initial = parse(Float64, in_vals[1])
    #         iteration = parse(Int, in_vals[2])
    #         # ψ = parse(Complex{Float64}, in_vals[3])  
    #         #println(readlines(f))
    #         # for lines in readlines(f)
    #         #     # print the line
    #         #     println("t_initial", lines)
    #         #     t_initial = parse(Float64, lines) 
    #         #     println("iteration", lines)
    #         #     iteration = parse(Int, lines)      
            
    #         # end
    #         close(f)
    #         # # Recover the MPS ψ from the binary file using ITensor's read function
    #         # psi_filename = "/home/zohalaraib/Oscillatrino/tests/checkpoints/checkpoint.chkpt.it000004.txt_psi.itensor"
    #         # if isfile(psi_filename)
    #         #     ψ = read(psi_filename, MPS)
    #         # else
    #         #     error("Checkpoint MPS file not found at $psi_filename. Ensure the checkpoint was saved properly.")
    #         # end
    #     elseif recover_type == "manual"
    #         println("recovery from iteration ", recover_iteration)
    #     end
    # end

    if do_recover
        if recover_type == "auto"
            println("auto recovery")
            checkpoint_filename = "/home/zohalaraib/Oscillatrino/tests/checkpoints/checkpoint.chkpt.it000004.h5"
            ψ, t_initial, iteration = recover_checkpoint_hdf5(checkpoint_filename)
            # s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t_initial, iteration = recover_checkpoint_hdf5(checkpoint_filename)
            # Increment t_initial by τ to ensure it starts from the next expected value
            t_initial += τ
            # Reinitialize necessary parameters based on recovered ψ
            @assert length(siteinds(ψ)) == N_sites "Mismatch in number of sites between recovered ψ and N_sites"
            println("Recovered ψ with length: ", length(siteinds(ψ)))

            # Normalize ψ to make sure it's valid after recovery
            normalize!(ψ)

            # Regenerate gates for the recovered ψ
            println("Generating recovery gates for recovered ψ...")
            gates_recover = create_gates(s, ψ, N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name, L, τ, energy_sign, periodic)

            # Ensure gates are properly generated
            println("Number of recovered gates: ", length(gates_recover))
            if length(gates_recover) == 0
                error("Recovered gates are empty. Please ensure the gates are correctly generated.")
            end
        elseif recover_type == "manual"
            println("recovery from iteration ", recover_iteration)
            # Optionally, add manual recovery logic for a specific iteration
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
    
    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
     for t in t_initial:τ:ttotal
        # extract the gates array generated in the gates_function file
        # gates = create_gates(s, ψ,N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name,L, τ, energy_sign, periodic)
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
        
        # if shape_name!=="none" #specfic to the inhomogenous case Test4
        #     # compute the avg expectation value of Sz at all sites
        #     sz_tot = expect(ψ, "Sz")
        #     sz = mean(sz_tot)
        # else 
            # compute expectation value of Sz (inbuilt operator in ITensors library) at the first site on the chain
            # sz = expect(ψ, "Sz"; sites=1)
        # end 

        if shape_name!=="none" #specfic to the inhomogenous case Test4
            # compute the avg expectation value of Sz at all sites
            sz_tot = expect(ψ, "Sz")  # Compute Sz for each site and store the values in sz_tot
            half_N = div(N_sites, 2)  # Calculate half of the number of sites
            sz = mean(sz_tot[1:half_N])  # Take the mean of the first half of the sz_tot array
        
        else 
            # compute expectation value of Sz (inbuilt operator in ITensors library) at the first site on the chain
            sz = expect(ψ, "Sz"; sites=1)
        end 

        # compute expectation value of sy and sx using S+ and S- (inbuilt operator in ITensors library) at the first site on the chain
        if p == zeros(N_sites, 3) #for rogerro's case only (b/c S+ S- needed to keep conservation of QN number)
            sy = -0.5 *im * (expect(complex(ψ), "S+"; sites=1) - expect(complex(ψ), "S-"; sites=1)) #re-check
            sx = 0.5 * (expect(ψ, "S+"; sites=1) + expect(ψ, "S-"; sites=1)) #recheck
        else 
            sy = expect(complex(ψ), "Sy"; sites=1)
            sx = expect(ψ, "Sx"; sites=1)
        end

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

        # Apply recovered gates if in recovery mode
        if do_recover && recover_type == "auto"
            println("Applying recovered gates...")
            @assert length(gates_recover) > 0 "Recovered gates are empty. Cannot proceed with empty gates."
            ψ = apply(gates_recover, ψ; cutoff, maxdim)
        else
            # Generate new gates for the current ψ
            println("Generating new gates for current ψ...")
            gates = create_gates(s, ψ, N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name, L, τ, energy_sign, periodic)
    
            # Ensure gates are properly generated before applying
            if length(gates) == 0
                error("Generated gates are empty. Please check the gate creation process.")
            end
    
            ψ = apply(gates, ψ; cutoff, maxdim)
        end

            
        # ψ = apply(gates, ψ; cutoff, maxdim)
        # ψ = tdvp(H, ψ,  -im *τ;   nsweeps=1,
        # reverse_step=true,outputlevel=1)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
        # println(typeof(ψ))
        println("iteration: ",  iteration, "  ", iteration % checkpoint_every)
        # if iteration % checkpoint_every == 0 
        #     println("CREATE CHECKPOINT AT ITERATION = ", iteration, " TIME = ", t)
        #     #checkpoint_filename = joinpath("checkpoint.chkpt.it.", lpad(iteration,6,"0"), ".txt")
        #     checkpoint_filename = joinpath(chkptdir, "checkpoint.chkpt.it" * lpad(iteration,6,"0") * ".txt")
        #     println("checkpoint_filename:", checkpoint_filename)
        #     checkpoint_simulation(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
        #     println("checkpoint 2 t:", t)
        #     # println(ψ)
        # end
        if iteration % checkpoint_every == 0 
            println("CREATE CHECKPOINT AT ITERATION = ", iteration, " TIME = ", t)
            checkpoint_filename = joinpath(chkptdir, "checkpoint.chkpt.it" * lpad(iteration, 6, "0") * ".h5")
            checkpoint_simulation_hdf5(checkpoint_filename, ψ, t, iteration)
            # checkpoint_simulation_hdf5(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
        end

        iteration = iteration + 1
    end
    t_array = 0.0:τ:ttotal

    # After the time evolution loop, calculate and print the growth rate from the ratio if both values have been captured
    if ρₑμ_at_t1 !== nothing && ρₑμ_at_t2 !== nothing
        global Im_Ω = (1/Δt ) * log(ρₑμ_at_t2 / ρₑμ_at_t1)
        println("Growth rate of flavor coherence of ρₑμ at t2 to ρₑμ at t1: $Im_Ω")
    else
        println("ρₑμ was not captured at both t1 and t2.")
    end

    if save_data
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
    end
    
    return Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, Im_Ω 
end


