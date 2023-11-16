using ITensors
using Plots
using Measures
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")
include("../src/momentum.jl")

# BTW this test is not finding the first minima(as was proposed in Rog paper), just the global minima. Discuss optimal conditions.
# So, need to put appropriate assert conditions.

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N_sites = 4 # number of sites # variable
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    τ = 0.05 # time step # variable
    ttotal = 150 # total time of evolution #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm  #variable
    Δp = 5 # width of shape function #variable
    del_m2 = 0 # fixed for 'only' self interactions

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 

    # # n = sum N_i /V # total number density
    # # N = n * V 
    # # n_i = N_/N_sites

    # # Sherwood's initial conditions: # throughout this cdoe i am assuming each site is occupied by a particlei.e. each site containsa umber of neutrinos of all sma eflvor 
    # # so all neutrinos are electron flavored while all muon flavored are anti neutrinos
    L = 1 # cm # domain size
    # n_mu_e =  4.891290848285061e+32 # cm^-3 # number density of electron flavor of neutrino
    # N_sites= 50 # total sites/ particles # particles evenly spaced for each flavor 


    # # big box for all 'sites'/particles
    # # total vol = V 
    # # length of each side = L 
    # V = L^3 
    # # N =  total no.of 'particles'
    # # n = total no.density of 'particles'
    # n = N / V

    # # small box (for each 'site')
    # # total vol = V_i 
    # # length of each side = Δx 
    # Δx^3 = V_i
    # # N_i =  total no.of 'particles' at site i 
    # # n_i  = total no.density of 'particles' at site i 
    # N_i = n_i  * V_i  

    # # Combining 
    # # sum n_i = n # total no.density 
    # # sum N_i = N # total no.particles
    # # sum V_i = V # total volume 
    # # sum Δx_i = L # domain size


    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    s = siteinds("S=1/2", N_sites; conserve_qns=true) #fixed
    
    # Initialize an array of ones for all N sites
    mu = ones(N_sites) # erg #variable
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    N = mu .* fill((Δx)^3/(sqrt(2) * G_F), N_sites) #fixed
    
    # Create a B vector which would be same for all N particles 
    B = [0, 0, -1] #variable

    function generate_x_array(N_sites, L)
        return [i * L / (2N_sites) for i in 1:N_sites]
    end
    
    x = generate_x_array(N_sites, L)
    #x = fill(rand(), N) #variable
    y = fill(rand(), N_sites) #variable
    z = fill(rand(), N_sites) #variable

    function generate_p_array(N_sites)
        half_N_sites = div(N_sites, 2)
        return [fill(10.0, half_N_sites); fill(-10.0, half_N_sites)]
    end

    # p matrix all with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites), generate_p_array(N_sites), generate_p_array(N_sites))
    #p = rand(N_sites, 3) # p array with random numbers for all components (x, y, z)

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Up" : "Dn") # Fixed to produce consistent results for the test assert conditions 

    #extract output for the survival probability values at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, N, B, N_sites, Δx, del_m2, p, x, Δp, ψ, shape_name, cutoff, tolerance, ttotal)

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("Survival probability vs t (self-interactions w geo+shape).pdf")
end 

@time main()

             