using ITensors
using Plots
using Measures
using Base.Threads
using LinearAlgebra

include("../src/my_MT_evolution.jl")
include("../src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

# function main(N; use_splitblocks = true,nsweeps = 10, blas_num_threads=1,
#         strided_num_threads=1, use_threaded_blocksparse=true, outputlevel=1)
    function main(N; use_splitblocks = true,nsweeps = 10, blas_num_threads=128,
        strided_num_threads=128, use_threaded_blocksparse=false, outputlevel=1)
    #N = 4 # number of sites 
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.05 # time step 
    ttotal = 10 # total time of evolution
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
    maxdim = 200 #bondimension

    ITensors.Strided.set_num_threads(strided_num_threads)
    BLAS.set_num_threads(blas_num_threads)
    

    ITensors.disable_threaded_blocksparse()


    if outputlevel > 0
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

    if outputlevel > 0
        @show sweeps
    end

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
    Sz_array, prob_surv_array, apply_time= evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, tolerance, ttotal,outputlevel, 
                                                                        use_splitblocks)
    energy,ψ  = eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,maxdim,use_splitblocks, 
                                                                sweeps,nsweeps,use_threaded_blocksparse)
    # if outputlevel > 0
    #     @show flux(ψ)
    #     @show maxlinkdim(ψ)
    # end


    #index of minimum of the prob_surv_array (containing survival probability values at each time step)
    i_min = argmin(prob_surv_array)
    # time at which the mimimum survival probability is reached
    t_min = τ * i_min - τ
    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)
    println("i_min= ", i_min)
    println("t_min= ", t_min)
    # # Check that our time of minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    #@assert abs(t_min - t_p_Rog) <  τ + tolerance 

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", 
           ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("Survival probability vs t (only self-interaction term plot).pdf")
    #return apply_time
end 
N= 4
println("Without threaded block sparse:\n")
@time main(N;nsweeps = 10, use_threaded_blocksparse=false)
println()