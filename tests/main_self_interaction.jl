using ITensors
using Plots
using Measures
include("src/expect.jl")
include("src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N = 4 # number of sites (NEED TO GO TILL 96 for Rog_results)
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.05 # time step (NEED TO BE 0.05 for Rog_results)
    ttotal = 5 # total time of evolution (NEED TO GO TILL 50 for Rog_results)
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    ∇x = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N; conserve_qns=true)  

    # Constants for Rogerro's fit (only self-interaction term)
    a_t = 0
    b_t = 2.105
    c_t = 0
    
    # Initialize an array of ones for all N sites
    mu = ones(N) # erg
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos 
    n = mu .* fill((∇x)^3/(sqrt(2) * G_F), N)
    
    # Create an array B with N elements. Each element of the array is a vector [0, 0, 1]
    B = fill([0, 0, 1], N)

    # Create an array ω with N elements. Each element of the array is zero.
    ω = fill(0, N) 

    # Initialize psi to be a product state (alternating down and up)
    ψ = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, n, ω, B, N, ∇x, ψ, cutoff, tolerance, ttotal)

    #index of minimum of the prob_surv_array (containing survival probability values at each time step)
    i_min = argmin(prob_surv_array)
    # time at which the mimimum survival probability is reached
    t_min = τ * i_min - τ
    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)
    println("i_min= ", i_min)
    println("t_min= ", t_min)
    # Check that our time of minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    @assert abs(t_min - t_p_Rog) <  τ + tolerance 

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("Survival probability vs t (only self-interaction term plot).pdf")
end 

@time main()

