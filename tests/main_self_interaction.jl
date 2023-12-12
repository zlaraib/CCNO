using ITensors
using Plots
using Measures
using Base.Threads
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")


# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main(; use_splitblocks = true, nsweeps, blas_num_threads,
    strided_num_threads, use_threaded_blocksparse, outputlevel)
    
    N = 4 # number of sites 
    maxdim = 200 #bondimension for the SVD in MPS representation
    cutoff = 1E-14 # also specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.05 # time step (NEED TO BE 0.05 for Rog_results)
    ttotal = 10 # total time of evolution (NEED TO GO TILL 50 for Rog_results)
    tolerance  = 5E-3 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 

    ITensors.Strided.set_num_threads(strided_num_threads)
    BLAS.set_num_threads(blas_num_threads)

    if outputlevel > 0

        if use_threaded_blocksparse
            ITensors.enable_threaded_blocksparse()
          else
            ITensors.disable_threaded_blocksparse()
        end      
        @show Threads.nthreads()
        @show Sys.CPU_THREADS
        @show BLAS.get_num_threads()
        @show ITensors.Strided.get_num_threads()
        @show ITensors.using_threaded_blocksparse()
        println()
    end

    sweeps = Sweeps(nsweeps)
    maxdims = min.([100, 200, 400, 800, 2000, 3000, maxdim], maxdim)
    maxdim!(sweeps, maxdims...)
    noise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
    @show sweeps

    # Constants for Rogerro's fit (only self-interaction term)
    a_t = 0
    b_t = 2.105
    c_t = 0
    
    # Initialize an array of ones for all N sites
    mu = ones(N) # erg
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    n = mu .* fill((Δx)^3/(sqrt(2) * G_F), N)
    
    # Create a B vector which would be same for all N particles 
    B = [0, 0, 1]

    # Create an array ω with N elements. Each element of the array is zero.
    ω = fill(0, N) 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N; conserve_qns=true)  

    # Initialize psi to be a product state (alternating down and up)
    ψ = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, prob_surv_array, apply_time=  evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal, nsweeps,
                                            maxdim,outputlevel, use_splitblocks, use_threaded_blocksparse)
    energy,ψ,apply_time_parallel  = eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,
                                    maxdim,use_splitblocks, sweeps,nsweeps,use_threaded_blocksparse)
    if outputlevel > 0
        @show flux(ψ)
        @show maxlinkdim(ψ)
    end

    #Serial code:
    if outputlevel ==0  
        # This function scans through the array, compares each element with its neighbors, 
        # and returns the index of the first local minimum it encounters. 
        # If no local minimum is found, it returns -1 to indicate that.
        function find_first_local_minima_index_serial(arr)
            n = length(arr)
            for i in 2:(n-1)
                if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                    return i
                end
            end
            return -1  
        end
        
        # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
        i_first_local_min = find_first_local_minima_index_serial(prob_surv_array)
        
        # Writing if_else statement to communicate if local minima (not) found
        if i_first_local_min != -1
            println("Index of the first local minimum: ", i_first_local_min)
        else
            println("No local minimum found in the array.")
        end

        # Time at which the first mimimum survival probability is reached
        t_min = τ * i_first_local_min - τ
        println("Corresponding time of first minimum index= ", t_min)

        # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
        t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
        println("t_p_Rog= ",t_p_Rog)

        # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
        @assert abs(t_min - t_p_Rog) <  τ + tolerance  # Commented out because assert condition dooesnt pass, hence test produces incorrect results.

        # Plotting P_surv vs t
        plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

        # Save the plot as a PDF file
        savefig("Survival probability vs t (serial self-interaction term plot).pdf")
    end

    #My attempt at parallelizing manually
    if outputlevel >0 && use_threaded_blocksparse==false
        function find_first_local_minima_index(arr::Vector{Atomic{Float64}})
            n = length(arr)
            if n < 3
                return -1  # Can't have a local minimum with less than 3 data points
            end
        
            # Iterate over the array, comparing each element to its neighbors
            for i in 2:(n-1)
                # Extract the float values for comparison
                prev_val = arr[i-1][]
                current_val = arr[i][]
                next_val = arr[i+1][]
                # Debugging print statements
                # println("Comparing: ", prev_val, " ", current_val, " ", next_val)
                if current_val < prev_val && current_val < next_val
                    return i
                end
            end
            return -1
        end
        # Call the function to find the first local minimum
        i_first_local_min = find_first_local_minima_index(prob_surv_array)
    
        # Output the results
        if i_first_local_min != -1
            println("Index of the first local minimum: ", i_first_local_min)
            # Calculate the corresponding time of the first minimum
            t_min = τ * (i_first_local_min - 1)  # Adjusting index to time
            println("Corresponding time of the first minimum: ", t_min)
            # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
            t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
            println("t_p_Rog= ",t_p_Rog)

            # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
            @assert abs(t_min - t_p_Rog) <  τ + tolerance 
        else
            println("No local minimum found in the array.")
        end

        # First, extract the Float64 values from the Atomic{Float64} array
        prob_surv_values = [val[] for val in prob_surv_array]

        # Generate the time steps for the x-axis
        time_steps = 0.0:τ:(τ*(length(prob_surv_values)-1))

        # Now plot the data
        plot(time_steps, prob_surv_values, xlabel = "t", 
            ylabel = "Survival Probability p(t)", legend = false, 
            size=(800, 600), aspect_ratio=:auto, margin=10mm)
        # Save the plot as a PDF file
        savefig("Survival probability vs t (my_MT self-interaction term plot).pdf")
    end

    # Multithreading via ITensors parallelization
    if outputlevel >0 && use_threaded_blocksparse == true
        # This function scans through the array, compares each element with its neighbors, 
        # and returns the index of the first local minimum it encounters. 
        # If no local minimum is found, it returns -1 to indicate that.
        function find_first_local_minima_index_threaded(arr)
            n = length(arr)
            for i in 2:(n-1)
                if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                    return i
                end
            end
            return -1  
        end
        
        # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
        i_first_local_min = find_first_local_minima_index_threaded(prob_surv_array)
        
        # Writing if_else statement to communicate if local minima (not) found
        if i_first_local_min != -1
            println("Index of the first local minimum: ", i_first_local_min)
        else
            println("No local minimum found in the array.")
        end

        # Time at which the first mimimum survival probability is reached
        t_min = τ * i_first_local_min - τ
        println("Corresponding time of first minimum index= ", t_min)

        # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
        t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
        println("t_p_Rog= ",t_p_Rog)

        # # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
        @assert abs(t_min - t_p_Rog) <  τ + tolerance
        # Plotting P_surv vs t
        plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", 
            ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

        # Save the plot as a PDF file
        savefig("Survival probability vs t (ITensors_MT self-interaction term plot).pdf")
    end

end 


println("Serial:\n")
@time main(;use_splitblocks = true,nsweeps=10, blas_num_threads=1,
strided_num_threads=1, use_threaded_blocksparse=false, outputlevel=0)
println("Multi-threading without threaded block sparse(my_effort):\n")
@time main(;use_splitblocks = true,nsweeps=10, blas_num_threads=128,
strided_num_threads=128, use_threaded_blocksparse=false, outputlevel=1)
println("Multi-threading with threaded block sparse(ITensors):\n")
@time main(;use_splitblocks = true,nsweeps=10, blas_num_threads=1,
strided_num_threads=1, use_threaded_blocksparse=true, outputlevel=1)
