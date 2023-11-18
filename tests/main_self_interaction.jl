using ITensors
using Plots
using Measures
using LinearAlgebra
#using NDTensors

include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")
include("../src/momentum.jl")

# BTW When you apply gates to an MPS, it will in general increase the bond dimension 
# (one exception is that single-site gates don’t change the bond dimension), 
# and if you continue applying gates without truncation the bond dimension will in general grow exponentially.
# There isn’t a way to set the truncation when you initialize the state
# this is all from julia itensors discourse answered by matthew fisherman


# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.
# throughout this code the words particle and sites are used interchangeably.
# Where each site is occupied by either some neutrinos or some antineutrinos. 

function main()
    N_sites = 8 # number of sites # variable
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    τ = 0.05 # time step # sec # variable
    ttotal = 10 # total time of evolution # sec #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    Δx = 1E-3 # length of the box of interacting neutrinos at a site in cm  #variable
    Δp = 5 # width of shape function  # cm #variable
    del_m2 = 0 # fixed for 'only' self interactions # (erg^2)
    maxdim = 5 # max bond dimension in MPS truncation

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 

    """
    PIC SETUP: 

    big box for all 'sites' or particles 
    total vol = V 
    length of each side = L 
    V = L^3 
    N_sites = total no.of sites or particles
    N =  total no.of neutrinos at all sites
    n = total no.density of neutrinos at all sites
    n = N / V 

    small box (for each 'site')
    total vol = V_i 
    length of each side = Δx 
    Δx^3 = V_i
    N_i =  total no.of neutrinos at site i 
    n_i  = total no.density of neutrinos at site i 
    N_i = n_i  * V_i  

    Combining 
    sum n_i = n # total no.density of neutrinos
    sum N_i = N # total no.of neutrinos
    sum V_i = V # total volume of the grid
    sum Δx_i = L # domain size

    """

    # Richers(2021) initial conditions:
    # throughout this code I am assuming each site is occupied by a particle i.e. each site contains some number of neutrinos all of same flavor 
    # so all neutrinos are electron flavored (at a site) which interact with electron flavored anti neutrinos (at a different site) in the opposing beam.
    L = 1 # cm # domain size
    n_mu_e =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    n_mu_e_bar =  4.891290848285061e+32 # cm^-3 # number density of electron flavor antineutrino
    # N_sites= 50 # total sites/particles that evenly spaced "for each (electron) flavor" 
    #N_sites = 100 # total particles/sites for all neutrino and anti neutrino electron flavored

    V = L^3 

    # Create an array of dimension N and fill it half with values of sites containing all electron neutrinos 
    # and other half with sites conatining electron anti-neutrino. 
    N_mu_e  = n_mu_e * V 
    N_1 = fill(N_mu_e / (N_sites ÷ 2), N_sites ÷ 2)
    N_mu_e_bar  = n_mu_e_bar * V 
    N_2 = fill(N_mu_e_bar / (N_sites ÷ 2), N_sites ÷ 2)
    N = vcat(N_1, N_2) # This is the total number of neutrinos. 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    s = siteinds("S=1/2", N_sites; conserve_qns=true) #fixed
    
    # Create a B vector which would be same for all N particles 
    B = [0, 0, -1] #variable

    function generate_x_array(N_sites, L)
        return [i * L / (2N_sites) for i in 1:N_sites]
    end
    
    x = generate_x_array(N_sites, L)
    y = fill(rand(), N_sites) #variable
    z = fill(rand(), N_sites) #variable

    function generate_p_array(N_sites)
        half_N_sites = div(N_sites, 2)
        return [fill(10.0e6, half_N_sites); fill(-10.0e6, half_N_sites)]
    end

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites), generate_p_array(N_sites), generate_p_array(N_sites))
    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= N_sites ÷ 2 ? 1 : -1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Up" : "Dn") # Fixed to produce consistent results for the test assert conditions 
    
    # # # METHOD 1 to perturb initial state
    # d = 2 # dimension for the MPS with N_site site indices
    # # Define a small perturbation strength
    # epsilon = randn(d^N_sites) # Adjust this as needed
    # ϵ_MPS = MPS(epsilon,s;cutoff,maxdim) # converting array to MPS 

    # # METHOD 2 to perturb initial state
    # # Initialize the state array
    # state = [n <= N_sites ÷ 2 ? "Dn" : "Up" for n in 1:N_sites]
    # ϵ_MPS = randomMPS(s, state, maxdim)

    # # Perturb the state by adding random noise to the amplitudes
    # perturbed_ψ = ψ .+ ϵ_MPS
    # ψ = perturbed_ψ 

    #extract output for the survival probability values at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, N, B, N_sites, Δx,del_m2, p, x, Δp, ψ, shape_name, energy_sign, cutoff, maxdim, tolerance, ttotal)

    # rho = outer(ψ', ψ)
    # println(rho)

    rho_ee = ( (2* Sz_array) + 1)/2
    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(rho_ee)-1), rho_ee, xlabel = "t", ylabel = "<rho_ee>", legend = false, size=(800, 600), aspect_ratio=:auto,margin= 10mm) 

    # Save the plot as a PDF file
    savefig("<rho_ee> vs t (self-interactions w geo+shape)_MF_FFI.pdf")
end 

@time main()
