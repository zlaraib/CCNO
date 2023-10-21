using ITensors
using Plots
using Measures
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")


# BTW this test is not finding the first minima(as was proposed in Rog paper), just the global minima. Discuss optimal conditions.
# So, need to put appropriate assert conditions.

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N = 4 # number of sites # variable
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    τ = 0.05 # time step # variable
    ttotal = 150 # total time of evolution #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm  #variable
    Δp = 5 # width of shape function #variable
    del_m2 = 0 # fixed for 'only' self interactions
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N; conserve_qns=true) #fixed
    
    # Initialize an array of ones for all N sites
    mu = ones(N) # erg #variable
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    n = mu .* fill((Δx)^3/(sqrt(2) * G_F), N)
    
    # Create a B vector which would be same for all N particles 
    B = [0, 0, -1] #variable

    x = fill(rand(), N) #variable
    y = fill(rand(), N) #variable
    z = fill(rand(), N) #variable

    p = rand(N, 3) # p array with random numbers for all components (x, y, z)

    # Initialize psi to be a product state (alternating down and up)
    ψ = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p_hat = momentum(p,N)

    #extract output for the survival probability values at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, n, B, N, Δx,del_m2, p, p_mod, p_hat, x, Δp, ψ, shape_name, cutoff, tolerance, ttotal)

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("Survival probability vs t (self-interactions w geo+shape).pdf")
end 

@time main()

