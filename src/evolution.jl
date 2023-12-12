include("gates_function.jl")  # Include the gates_functions.jl file

"""Expected units of the quantities defined in the files in tests directory that are being used in the evolve function                                                                   
s = site index array (dimensionless and unitless) 
τ = time step (sec)      
n = no.of neutrinos (dimensionless and unitless)
ω = vacuum oscillation angular frequency (rad/s)
B = Normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
N = Total no.of sites (dimensionless and unitless)
Δx = length of the box of interacting neutrinos at a site (cm) 
cutoff = truncation threshold for the SVD in MPS (unitless, number)
ttotal = ttotal time (sec)

This test uses the evolve function to evolve the ψ state in time and compute the expectation values of Sz at each time step, along 
with their survival probabilities.The time evolution utilizes the unitary operators created as gates from the create_gates function.
The <Sz> and Survival probabilities output from this function are unitless. 
"""
mutable struct EntanglementObserver <: AbstractObserver
end

function ITensors.measure!(o::EntanglementObserver; bond, ψ, half_sweep, kwargs...)
    wf_center, other = half_sweep==1 ? (ψ[bond+1],ψ[bond]) : (ψ[bond],ψ[bond+1])
    U,S,V = svd(wf_center, uniqueinds(wf_center,other))
    SvN = 0.0
    sz=0.0
    for n=1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
      #sz =expect(ψ, "Sz"; sites=1)
    end
    println("  Entanglement across bond $bond = $SvN")
  end

  function evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal, nsweeps,
    maxdim,outputlevel, use_splitblocks, use_threaded_blocksparse)
    # Initialize shared arrays
    if outputlevel >0 && use_threaded_blocksparse==false
        Sz_array = Atomic{Float64}[]
        prob_surv_array = Atomic{Float64}[]
    else    
         # Create empty array to store sz values 
        Sz_array = Float64[]
        # Create empty array to store survival probability values 
        prob_surv_array = Float64[]
    end

    gates = create_gates(s, n, ω, B, N, Δx, τ, outputlevel, use_splitblocks)

    apply_time = 0.0

    for t in 0.0:τ:ttotal
        if outputlevel >0 && use_threaded_blocksparse==false
            # Compute initial expectation value of Sz and store as an Atomic{Float64} for thread safety.
            sz = expect(ψ, "Sz"; sites=1)
            push!(Sz_array, Atomic{Float64}(sz))

            prob_surv = 0.5 * (1 - 2 * sz)
            push!(prob_surv_array, Atomic{Float64}(prob_surv))
        
        else 
            # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site on the chain
            sz = expect(ψ, "Sz"; sites=1)
            # add an element sz to the end of Sz array 
            push!(Sz_array, sz)
            
            
            # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
            prob_surv = 0.5 * (1 - 2 * sz)
            # add an element prob_surv to the end of  prob_surv_array 
            push!(prob_surv_array, prob_surv)
        end

        if ω == fill(0, N) 
            println("$t $prob_surv")
        else println("$t $sz")
        end

        # If the current time t is approximately equal to ttotal, the loop breaks.
        t ≈ ttotal && break
        if outputlevel >0 && use_threaded_blocksparse==false
            # Use an array to hold the new MPSs from each thread
            new_psis = Array{MPS, 1}(undef, nthreads())

            # Divide gates among threads and compute in parallel
            num_threads = nthreads()
            gates_per_thread = ceil(Int, length(gates) / num_threads)

            # sync macro ensures that the enclosed block of code executes in a synchronized manner
            # It waits for all the spawned threads within its scope to complete before moving on
            Threads.@sync begin
                # The @threads macro is used to distribute the iterations of the loop across multiple threads.
                Threads.@threads for i in 1:num_threads
                    # spawn macro is used to create a new thread for executing the code block which is nested inside the @threads loop, 
                    # i.e. for each iteration of the loop, a new thread is spawned to execute a specific apply to relevant sites.
                    Threads.@spawn begin
                        # Compute the application of a subset of gates for each thread
                        # deepcopy creates an entirely independent copy of the entire structure of ψ
                        # This is important in multi-threaded applications to prevent race conditions where multiple threads are trying to read from and write to the same memory location
                        thread_psi = deepcopy(ψ)
                        start_idx = (i - 1) * gates_per_thread + 1
                        end_idx = min(i * gates_per_thread, length(gates))
                        # dividing a collection of gates into subsets. 
                        # Each thread works on a different subset of gates, determined by start_idx and end_idx. 
                        gates_subset = gates[start_idx:end_idx]

                        for gate in gates_subset
                            # Use the ITensor apply function to handle the gate application
                            apply_time += @elapsed thread_psi = apply(gate, thread_psi; cutoff,maxdim)
                        end

                        # Store the resulting MPS in the array (i.e. result for each thread (thread_psi) is stored in an array new_psis)
                        new_psis[i] = thread_psi
                    end
                end
            end

            # Combine the new MPSs from each thread
            ψ = combine_mps_states(new_psis,cutoff, nsweeps)
        end
        
        if outputlevel >0 && use_threaded_blocksparse== true
            # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
            # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
            # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
            @disable_warn_order begin
                ITensors.enable_combine_contract()
                ITensors.enable_threaded_blocksparse()
                apply_time += @elapsed ψ = apply(gates, ψ;cutoff)
                ITensors.disable_threaded_blocksparse()
                ITensors.disable_combine_contract()
            end
        end

        if outputlevel == 0
        # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
        apply_time += @elapsed ψ = apply(gates, ψ; cutoff)
        end 
        # to make sure the quantum state remains physical 
        normalize!(ψ)
        end
        
        @show apply_time

    return Sz_array, prob_surv_array, apply_time
end


function combine_mps_states(mps_array::Array{MPS, 1},cutoff, nsweeps)
    # Start with the first state in the array as the base for the superposition
    combined_state = mps_array[1]

    # Determine the maximum bond dimension we are willing to allow for the combined state
    max_bond_dim = determine_max_bond_dim(combined_state, mps_array)

    for i in 2:length(mps_array)
        # Superpose mps_array[i] onto combined_state
        # So we construct a new MPS with a probable higher bond dimension
        # that captures the superposition of combined_state and mps_array[i].
        # After creating the superposition, we truncate the MPS to the max_bond_dim
        # to keep the computation tractable.

        combined_state = superpose_and_truncate(combined_state, mps_array[i], max_bond_dim,cutoff, nsweeps)
    end

    # Normalize the combined state to ensure it represents a valid quantum state
    normalize!(combined_state)

    return combined_state
end

function determine_max_bond_dim(base_state::MPS, mps_array::Array{MPS, 1})
    # since all states are initialized with maxdim=200.
    maxdim=200
    return maxdim
end

function superpose_and_truncate(state1::MPS, state2::MPS, max_bond_dim::Int, cutoff, nsweeps)
    # First, create a new MPS that is the sum of state1 and state2
    # This operation would involve summing the tensors at each site
    # which could increase the bond dimension beyond what is tractable
    combined_mps = add_mps(state1, state2)
    
    # Then, truncate the combined MPS to the maximum bond dimension
    # with the specified cutoff to keep the computation manageable
    truncated_mps = truncate_mps(combined_mps, max_bond_dim, cutoff, nsweeps)

    return truncated_mps
end

function add_mps(state1::MPS, state2::MPS)
    # Verify that the MPS have the same length and site indices
    if length(state1) != length(state2)
        throw(DimensionMismatch("MPS must have the same length"))
    end
    
    # Retrieve site indices for each MPS
    sites1 = siteinds(state1)
    sites2 = siteinds(state2)
    
    # Check if the site indices match
    for i in 1:length(sites1)
        if sites1[i] != sites2[i]
            throw(DimensionMismatch("MPS site indices must match"))
        end
    end
    
    # Proceed with MPS addition assuming the site indices match
    # Combine the state tensors from state1 and state2
    combined_state = state1 + state2  # This is not correct as adding MPS tensors site by site is not directly possible with ITensor's MPS type
    # something correct would be like this where we create a new MPS with the correct site indices here and contract over matching indices
    # for i in 1:length(state1)
    #     combined_state[i] = state1[i] + state2[i]
    # end
    
    return combined_state
end

function truncate_mps(mps::MPS, max_bond_dim::Int, cutoff::Float64,nsweeps)
    # Apply truncation to the MPS to reduce its maximum bond dimension via cutoff and maxdim 
    sweep = Sweeps(nsweeps)  # One sweep should be enough for truncation
    maxdim!(sweep, max_bond_dim)
    cutoff!(sweep, cutoff)

    return mps
end

function eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,
    maxdim,use_splitblocks, sweeps,nsweeps,use_threaded_blocksparse)
    
    H = Hamiltonian_mpo(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
    state = Vector{String}(undef, N)
    ψ = nothing
    for i in 1:N
        # Assign "Up" or "Dn" based on even/odd index
        if i % 2 == 0
            state[i] = "Up"
        else
            state[i] = "Dn"
        end
    end
    
    #create a MPS state with half spins pointing up and half pointing down
    ψ0 = randomMPS(s, state,4)
    # ψ0 = productMPS(s, n -> isodd(n) ? "Dn" : "Up")
    #ψ0=randomMPS(s,linkdims=4)
    apply_time_parallel = 0

    #obs = EntanglementObserver()

    Sz_observer = DMRGObserver(["Sz"],s,energy_tol=1E-7)
    observer = EntanglementObserver()
    
    ITensors.enable_combine_contract()
    #apply_time_parallel += @elapsed energy, ψ =@time dmrg(H,ψ0; nsweeps, maxdim, cutoff, observer=Sz_observer, outputlevel=2)
    apply_time_parallel += @elapsed energy, ψ = @time dmrg(H, ψ0, sweeps; outputlevel)
    normalize!(ψ)
    ITensors.disable_combine_contract()
    
    for (sw,Szs) in enumerate(measurements(Sz_observer)["Sz"])
        println("Total Sz after sweep $sw = ", sum(Szs)/N)
    end

    return energy,ψ,apply_time_parallel
end

