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


function evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal,outputlevel,maxdim,use_splitblocks, use_threaded_blocksparse)
    
    # Create empty array to store sz values 
    Sz_array = Float64[]
    # Create empty array to store survival probability values 
    prob_surv_array = Float64[]

    # extract the gates array generated in the gates_function file
    gates = create_gates(s, n, ω, B, N, Δx, τ, outputlevel,use_splitblocks)
   
    apply_time = 0
    # println(gates)
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

        # if ω == fill(0, N) 
        #     println("$t $prob_surv")
        # else println("$t $sz")
        # end

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break



                    # function apply(gates, ψ;cutoff)
            #     temp = zeros(Int, Threads.nthreads())
            #     Threads.@threads for i in eachindex(ψ)
            #         temp[Threads.threadid()] = ψ[i]
            #     end
            #     return temp
            # end
            # function apply_gate(mat, vec)
            #     num_rows, num_cols = size(mat)
            #     res = zeros(ComplexF64, num_rows)
                
            #     for i in 1:num_rows
            #         for j in 1:num_cols
            #             res[i] += mat[(i - 1) * num_cols + j] * vec[j]
            #         end
            #     end
                
            #     return res
            # end
            # function apply_g(gates, ψ; cutoff)
            #     temp = Vector{ComplexF64}(undef, length(ψ))
            #     Threads.@threads for i in eachindex(ψ)
            #         local_state = ψ[i]
            #         for gate in gates
            #             local_state = apply_gate(gate, local_state)
            #         end
            #         temp[i] = local_state
            #     end
            #     return temp
            # end

            # function product(
            #     o::ITensor,
            #     ψ::AbstractMPS,
            #     ns=findsites(ψ, o);
            #     move_sites_back::Bool=true,
            #     apply_dag::Bool=false,
            #     kwargs...,
            #   )
            #     N = length(ns)
            #     ns = sort(ns)
              
            #     # TODO: make this smarter by minimizing
            #     # distance to orthogonalization.
            #     # For example, if ITensors.orthocenter(ψ) > ns[end],
            #     # set to ns[end].
            #     ψ = orthogonalize(ψ, ns[1])
            #     diff_ns = diff(ns)
            #     ns′ = ns
            #     if any(!=(1), diff_ns)
            #       ns′ = [ns[1] + n - 1 for n in 1:N]
            #       ψ = movesites(ψ, ns .=> ns′; kwargs...)
            #     end
            #     ϕ = ψ[ns′[1]]
            #     for n in 2:N
            #       ϕ *= ψ[ns′[n]]
            #     end
            #     ϕ = product(o, ϕ; apply_dag=apply_dag)
            #     ψ[ns′[1]:ns′[end], kwargs...] = ϕ
            #     if move_sites_back
            #       # Move the sites back to their original positions
            #       ψ = movesites(ψ, ns′ .=> ns; kwargs...)
            #     end
            #     return ψ
            #   end
              
            #   apply = product
    
    # end

  #   # function mpsapply(ψ::MPS, gates::Vector{ITensor}, dt::Real, t::Real; effect! = nothing, savelast::Bool = false, kwargs...) #keyword arguments for ITensors.apply
  #     # out = [deepcopy(ψ)]
  #     # for _ in dt:dt:ttotal
  #         savelast::Bool = false
  #         if savelast
  #           ψ[1] = apply(gates, ψ[1]; normalize = true, cutoff)
  #           normalize!(ψ[1])
  #         else
  #             push!(ψ, apply(gates, ψ[end]; normalize = true, cutoff))
  #             normalize!(ψ)
  #         end
  #         !isa(effect!, Nothing) ? effect!(ψ[end]) : nothing
  #     # end
  #     # return out
  # # end


        # function apply(gates,ψ;cutoff)
        #     temp = zeros(Int, Threads.nthreads())
        #     shared_temp = zeros(Int, Threads.nthreads())
        #     Threads.@threads for i in findsites(ψ, gates)
        #         temp[Threads.threadid()]+= ψ[i]
        #         shared_temp[Threads.threadid()]+= gates[i]
        #         # push!(temp, temp[i])
        #         # normalize(temp)

        #     end
        #     # Threads.@threads for i in findsites(ψ, gates)
                
        #     #     # push!(temp, temp[i])
        #     #     # normalize(temp)

        #     # end

        #     ψ = apply(shared_temp, temp; cutoff)
        #     return ψ
        # end

        # function apply(gates, ψ; cutoff)
        #     num_threads = Threads.nthreads()
        #     temp = [deepcopy(ψ) for _ in 1:num_threads]  # Create copies of ψ for each thread
        #     shared_temp = [deepcopy(gates) for _ in 1:num_threads] 
        
        #     Threads.@threads for i in length(findsites(ψ, gates))
        #         thread_id = Threads.threadid()
        #         temp[thread_id][i] += gates[i]
        #         shared_temp[thread_id][i] += gates[i]
        #     end
        
        #     # Process temp and shared_temp for each thread separately
        #     for thread_id in 1:num_threads
        #         temp[thread_id] = apply(shared_temp[thread_id], temp[thread_id]; cutoff)
        #     end
        
        #     # Collect all results into a single array
        #     all_results = Int[]  # Initialize an empty array to collect all results
        #     for thread_temp in temp
        #         append!(all_results, thread_temp)  # Append each thread's result to all_results
        #     end

        #     return all_results
        # end

        # function apply(gates, ψ; cutoff)
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
        #         all_results = combine_mps(all_results, results[i])  # Replace combine_mps with appropriate operation
        #     end
        
        #     return all_results
        # end
        

        # function apply_gate!(psi::MPS,G::ITensor ; kwargs...)
        #     b = bond(G)
        #     orthogonalize!(psi,b)
        #     wf = psi[b]*psi[b+1]
        #     wf = noprime( G*(psi[b]*psi[b+1]) )
        #     spec = replacebond!(psi, b, wf; normalize=get(kwargs, :normalize,true), kwargs...)
        #     return spec
        # end
        
        # function apply_gates!(psi::MPS, Gs::Vector{ITensor}; kwargs...)
        #     rev = get(kwargs,:reversed,false)
        #     cb_func = get(kwargs,:cb, nothing)
        #     Gs_ordered = rev ? (@view Gs[end:-1:1]) : Gs
        #     for g in Gs_ordered
        #         spec=apply_gate!(psi,g; kwargs...)
        #         !isnothing(cb_func) && cb_func(psi; bond=bond(g),spec=spec)
        #     end
        # end
        
        # ψ = apply_gates!( ψ, gates;cutoff)
        

        # function apply(gates, ψ; cutoff)
        #     results = deepcopy(ψ)
        #     Threads.@threads for i in siteinds(commoninds, gates,ψ)
        #         result = apply_gate(gates[i],ψ[i]; cutoff)
        #         #results[i] = normalize!(result)
        #     end
        #     return results
        # end
        # function apply_gate(gates, ψ; cutoff)
        #     ψ = apply(gates, ψ; cutoff)
        #     normalize!(ψ)
        #     return ψ
        # end


        # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
        
        #ψ = mpsapply(gates, ψ; cutoff)
        @disable_warn_order begin
            ITensors.enable_threaded_blocksparse()
            ψ = apply(gates, ψ;cutoff, maxdim)
            ITensors.disable_threaded_blocksparse()
        end
        # println(ψ)
        #apply_time += @elapsed ψ = apply(gates, ψ; cutoff)
        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end
    #@show apply_time
    return Sz_array, prob_surv_array
end



# function evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal, outputlevel)
#     Sz_array = Float64[]
#     prob_surv_array = Float64[]
#     gates = create_gates(s, n, ω, B, N, Δx, τ, outputlevel)
    
#     for t in 0.0:τ:ttotal
#         sz = expect(ψ, "Sz"; sites=1)
#         push!(Sz_array, sz)
        
#         prob_surv = 0.5 * (1 - 2 * sz)
#         push!(prob_surv_array, prob_surv)
        
#         # Parallelizing the application of gates
#         threads = Task[]
#         for gate in gates
#             th = Threads.@spawn begin
#                 new_psi = apply([gate], deepcopy(ψ); cutoff)
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
    
#     return Sz_array, prob_surv_array
# end

