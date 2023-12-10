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


function evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal,
    outputlevel,maxdim,use_splitblocks, use_threaded_blocksparse)

    # Create empty array to store sz values 
    Sz_array = Float64[]
    # Create empty array to store survival probability values 
    prob_surv_array = Float64[]
    # extract the gates array generated in the gates_function file
    gates = create_gates(s, n, ω, B, N, Δx, τ, outputlevel,use_splitblocks)
   
    apply_time = 0
    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
    for t in 0.0:τ:ttotal
        
        # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site on the chain
        sz = expect(ψ, "Sz"; sites=1)
        # add an element sz to the end of Sz array 
        push!(Sz_array, sz)
        
        # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
        prob_surv = 0.5 * (1 - 2 * sz)
        # add an element prob_surv to the end of  prob_surv_array 
        push!(prob_surv_array, prob_surv)

        # for i in 1:(N-1)
        #     if ω[i] != 0
        #         # Compute the expected value based on the derived analytic formula
        #         expected_sz = 0.5 * cos(ω[i] * t)
                
        #         # Checking that the value of Sz at the center spin oscillates between 0.5 and -0.5 
        #         # Compare the actual value with the expected value using a tolerance
        #         @assert abs(sz - expected_sz) < tolerance
        #     end
        # end

        if ω == fill(0, N) 
            println("$t $prob_surv")
        else println("$t $sz")
        end

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break
        
        # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
        @disable_warn_order begin
            #ITensors.enable_combine_contract()
            
            apply_time += @elapsed ψ = apply_z(gates, ψ;cutoff)
            #ITensors.disable_combine_contract()
        end
        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end
     
    return Sz_array, prob_surv_array, apply_time
end

# function apply_z(gates, ψ; cutoff)
#     num_threads = Threads.nthreads()
#     temp = [deepcopy(ψ) for _ in 1:num_threads]  # Create copies of ψ for each thread
#     shared_temp = [deepcopy(gates) for _ in 1:num_threads]  # Create copies of gates for each thread
#     results = [deepcopy(ψ) for _ in 1:num_threads]
#     Threads.@threads for i in findsites(ψ, gates)
#         thread_id = Threads.threadid()
#         temp[thread_id][i] += ψ[i]
#         shared_temp[thread_id][i] += gates[i]
#         # Perform operations using temp and shared_temp for each thread
#     end

#     # Process temp and shared_temp for each thread separately
#     for thread_id in 1:num_threads
#         results[thread_id] = apply(shared_temp[thread_id], temp[thread_id]; cutoff)
#         # Additional processing if needed
#     end

#     # Construct the final result as an MPS
#     all_results = deepcopy(ψ)  # Initialize all_results with the first thread's result
#     for i in 2:num_threads
#         # Concatenate/combine the MPS from other threads into all_results
#         # Adjust this part based on how you want to combine the MPS objects
#         all_results = vcat(all_results, results[i])  # Replace combine_mps with appropriate operation
#     end

#     return all_results
# end

function product(
    o::ITensor,
    ψ::MPS,
    ns=findsites(ψ, o);
    move_sites_back::Bool=true,
    apply_dag::Bool=false,
    kwargs...,
  )
    N = length(ns)
    ns = sort(ns)
  
    # TODO: make this smarter by minimizing
    # distance to orthogonalization.
    # For example, if ITensors.orthocenter(ψ) > ns[end],
    # set to ns[end].
    ψ = orthogonalize(ψ, ns[1])
    diff_ns = diff(ns)
    ns′ = ns
    if any(!=(1), diff_ns)
      ns′ = [ns[1] + n - 1 for n in 1:N]
      ψ = movesites(ψ, ns .=> ns′; kwargs...)
    end
    ϕ = ψ[ns′[1]]
    for n in 2:N
      ϕ *= ψ[ns′[n]]
    end
    ϕ = product(o, ϕ; apply_dag=apply_dag)
    ψ[ns′[1]:ns′[end], kwargs...] = ϕ
    if move_sites_back
      # Move the sites back to their original positions
      ψ = movesites(ψ, ns′ .=> ns; kwargs...)
    end
    return ψ
  end
  gates = create_gates(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
  apply_z(gates, ψ, kwargs...) = product(
    o::ITensor,
    ψ::AbstractMPS,
    ns=findsites(ψ, o);
    move_sites_back::Bool=true,
    apply_dag::Bool=false,
    kwargs...,
  )

# function evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal,outputlevel, 
#     maxdim,use_splitblocks,use_threaded_blocksparse)
#     Sz_array = Float64[]
#     prob_surv_array = Float64[]
#     gates = create_gates(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
#     apply_time = 0.0
#     for t in 0.0:τ:ttotal
#         sz = expect(ψ, "Sz"; sites=1)
#         push!(Sz_array, sz)
        
#         prob_surv = 0.5 * (1 - 2 * sz)
#         push!(prob_surv_array, prob_surv)
        
#         # Parallelizing the application of gates
#         threads = Task[]
#         for gate in gates
#             th = Threads.@spawn begin
#                 apply_time += @elapsed new_psi = apply(gate, deepcopy(ψ); cutoff)
#                 normalize!(new_psi)
#                 #println(new_psi)
#                 return new_psi
#             end
#             push!(threads, th)
#         end
        
#         # Waiting for all threads to finish and updating ψ
#         for thread in threads
#             ψ = fetch(thread)
#         end
#         println("$t $prob_surv")
#         #t ≈ ttotal && break
#     end
    
#     return apply_time
# end

# function apply_z(
#     As::Vector{<:ITensor}, ψ::MPS; move_sites_back::Bool=true, kwargs...
#   )
#     Aψ = ψ
#     for A in As
#       Aψ = apply_z(A, Aψ; move_sites_back=false, kwargs...)
#     end
#     if move_sites_back
#       s = siteinds(Aψ)
#       ns = 1:length(ψ)
#       ñs = [findsite(ψ, i) for i in s]
#       Aψ = movesites(Aψ, ns .=> ñs; kwargs...)
#     end
#     return Aψ
#   end


function eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,maxdim,use_splitblocks, sweeps,nsweeps,use_threaded_blocksparse)
    
    # cutoff!(sweeps, 1e-14)
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
    ψ0 = randomMPS(s, state,4)
    # ψ0 = productMPS(s, n -> isodd(n) ? "Dn" : "Up")
    #ψ0=randomMPS(s,linkdims=4)
    apply_time_parallel = 0

    

    #obs = EntanglementObserver()

    

    Sz_observer = DMRGObserver(["Sz"],s,energy_tol=1E-7)
    observer = EntanglementObserver()
    
    #ITensors.enable_combine_contract()
    #apply_time_parallel += @elapsed energy, ψ =@time dmrg(H,ψ0; nsweeps, maxdim, cutoff, observer=Sz_observer, outputlevel=2)
    apply_time_parallel += @elapsed energy, ψ = @time dmrg(H, ψ0, sweeps; outputlevel)
    # normalize!(ψ)
    #ITensors.disable_combine_contract()
    
    for (sw,Szs) in enumerate(measurements(Sz_observer)["Sz"])
        println("Total Sz after sweep $sw = ", sum(Szs)/N)
    end


     
    return energy,ψ,apply_time_parallel
end

