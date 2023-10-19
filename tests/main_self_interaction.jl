using ITensors
using Plots
using Measures
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")


# BTW this test is not finding the first minima(as was proposed in Rog paper), just the global minima. Discuss optimal conditions.
# So, need to put appropriate assert conditions and define and use proper unit vectors. This is not done yet and will produce wrong results. 
# Also write ω in terms of norm(p) and remove dependence on \omega in functions that contribute through the vacuum oscillations. But in this test is a problem b/c for the self inetractions part we consider the ω to be zero for all sites which would mean that norm(p) should blow up for all sites. Weird. Discuss.


# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N = 4 # number of sites (NEED TO GO TILL 96 for Rog_results)
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.05 # time step (NEED TO BE 0.05 for Rog_results)
    ttotal = 5 # total time of evolution (NEED TO GO TILL 50 for Rog_results)
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N; conserve_qns=true)  

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
    x = fill(rand(), N)
    y = fill(rand(), N)
    z = fill(rand(), N)

    Δp = 5 #width of shape function
    p = zeros(N, 3) # Initialize the p array with zeros for all components (x, y, z)
    for i in 1:N
        p[i, 1] = rand() #p_x array
        p[i, 2] = rand() #p_y array
        p[i, 3] = rand() #p_z array
    end
    println(p) 
    unit_p = p/norm(p)
    println(norm(p))
    println(unit_p) 
    del_m2 = 2.0
    #norm(p) = [del_m2 / (2 * ω[i]) for i in 1:N]
    # The norm function returns the magnitude of the vector (by default):
    # A unit vector is then found by scaling by the reciprocal of the magnitude: e.g array/vector v  can have unit vector v to be  defined as: v / norm(v)

    # Initialize psi to be a product state (alternating down and up)
    ψ = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, n, ω, B, N, Δx, p, x, Δp, ψ, cutoff, tolerance, ttotal)

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
    #@assert abs(t_min - t_p_Rog) <  τ + tolerance 

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("Survival probability vs t (only self-interaction term plot).pdf")
end 

@time main()

