using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")
include("../src/momentum.jl")
include("../src/perturb.jl")
"""
PIC SETUP: 

big box (for 'all' sites/particles)
total vol = V 
length of each side = L 
V = L^3 
N_sites = total no.of sites or particles
N =  total no.of neutrinos at all sites
n = total no.density of neutrinos at all sites
n = N / V 

small box (for 'each' site/particle)
total vol = V_i 
length of each side = Δx 
Δx^3 = V_i
N_i =  total no.of neutrinos at site i 
n_i  = total no.density of neutrinos at site i 
N_i = n_i  * V_i  

Combining (where index i represent a site and runs from 1:N_sites)
sum n_i = n # total no.density of neutrinos
sum N_i = N # total no.of neutrinos
sum V_i = V # total volume of the grid
sum Δx_i = L # domain size

"""

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

# throughout this code I am assuming each site is occupied by a particle i.e. each site contains some number of neutrinos all of same flavor 
# so all neutrinos are electron flavored (at a site) which interact with electron flavored anti neutrinos (at a different site) in the opposing beam.
function main()
 
    # my initial conditions 
    N_sites =2  # number of sites # variable
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    τ = 1e-10 # time step # sec # variable
    ttotal = 1e-9 # total time of evolution # sec #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    Δp = 1/4 # width of shape function  # cm #variable
    del_m2 = 0 # fixed for 'only' self interactions # (erg^2)
    maxdim = 1 # max bond dimension in MPS truncation
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 
    periodic = true  # true = imposes periodic boundary conditions while false doesn't
    function generate_p_array(N_sites)
        half_N_sites = div(N_sites, 2)
        return [fill(50.0e6, half_N_sites); fill(-50.0e6, half_N_sites)]
    end
    L = 1 # cm # domain size # (aka big box length)
    n_mu_e =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    n_mu_e_bar =  4.891290848285061e+32 # cm^-3 # number density of electron flavor antineutrino

    # #Richers(2021) initial conditions:
    # L = 1 # cm # domain size # (aka big box length)
    # n_mu_e =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    # n_mu_e_bar =  4.891290848285061e+32 # cm^-3 # number density of electron flavor antineutrino
    # N_sites_eachflavor= 50 # total sites/particles that evenly spaced "for each (electron) flavor" 
    # N_sites = 100 #(2* N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
    # τ = 6.5E-13 # time step # sec # variable
    # ttotal = 9E-11 # total time of evolution # sec #variable
    # function generate_p_array(N_sites)
    #     half_N_sites = div(N_sites, 2)
    #     return [fill(50.0e6, half_N_sites); fill(-50.0e6, half_N_sites)]
    # end

    V = L^3 
    Δx = L/N_sites # length of the box of interacting neutrinos at a site in cm  #variable
    # Create an array of dimension N and fill it half with values of sites containing all electron neutrinos 
    # and other half with sites conatining electron anti-neutrino. 
    N_mu_e  = n_mu_e * V 
    N_1 = fill(N_mu_e / (N_sites ÷ 2), N_sites ÷ 2)
    N_mu_e_bar  = n_mu_e_bar * V 
    N_2 = fill(N_mu_e_bar / (N_sites ÷ 2), N_sites ÷ 2)
    N = vcat(N_1, N_2) # This is the total number of neutrinos. 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.
    s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function
    
    # Create a B vector that allows for perturbation to inital state in different directions
    B = [0.02, -0.02, 1] #variable
    # Normalize B to have a norm of 1
    B = B / norm(B)

    # generate x_array such that the first particle is at position L/(2*N_sites) while subsequent particles are at a position incremental by L/N_sites. # grid style
    function generate_x_array(N_sites, L)
        return [(i - 0.5) * L / N_sites for i in 1:N_sites]
    end
    
    x = generate_x_array(N_sites, L)
    y = fill(rand(), N_sites) #variable
    z = fill(rand(), N_sites) #variable


    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites), generate_p_array(N_sites), generate_p_array(N_sites))
    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= N_sites ÷ 2 ? 1 : -1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Up" : "Dn") # Fixed to produce consistent results for the test assert conditions 
    
    # Perturb the state via one-body Hamiltonian
    ψ_0 = evolve_perturbation(s, τ, B, N_sites, ψ, cutoff, maxdim, ttotal)

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "..","misc","datafiles","FFI", "par_"*string(N_sites), "tt_"*string(ttotal))

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, px_values, ρ_ee_array= evolve(s, τ, N, B,L, N_sites, 
                    Δx,del_m2, p, x, Δp, ψ_0, shape_name, energy_sign, cutoff, maxdim, datadir, ttotal,periodic)

    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","FFI", "par_"*string(N_sites), "tt_"*string(ttotal))
    
    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)
    # Plotting ρ_ee vs t
    plot(0.0:τ:τ*(length(ρ_ee_array)-1), ρ_ee_array, xlabel = "t", ylabel = "<ρ_ee>", legend = false, 
    left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 

    # Save the plot as a PDF file
    savefig(joinpath(plotdir, "<ρ_ee>_vs_t_self-interactions_w_geo+shape_MF_FFI.pdf"))

    # Plotting Sz vs t
    plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false,
        left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm,margin= 10mm,
        ylims = (minimum(Sz_array), maximum(Sz_array) + 0.1 * abs(maximum(Sz_array) - minimum(Sz_array)))
    )

    # Save the plot as a PDF file
    savefig(joinpath(plotdir, "<Sz>_vs_t_self-interactions_w_geo+shape_MF_FFI.pdf"))

    plot(0.0:τ:τ*(length(Sy_array)-1), Sy_array, xlabel = "t", ylabel = "<Sy_array>", legend = false,
    left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sy> vs t (self-interactions w geo+shape)_MF_FFI.pdf"))

    plot(0.0:τ:τ*(length(Sx_array)-1), Sx_array, xlabel = "t", ylabel = "<Sx_array>", legend = false,
    left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sx> vs t (self-interactions w geo+shape)_MF_FFI.pdf"))

    plot(title="Particle Position Evolution", xlabel= "Position (x)",ylabel="Time")
    for site in 1:N_sites
        site_positions = [(x_values[t][site]) for t in 1:length(x_values)]
        plot!(site_positions, 0.0:τ:ttotal, label="Site $site",
        left_margin = 30mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm, margin= 10mm)
    end

    savefig(joinpath(plotdir,"Particles evolution.pdf"))

end 

@time main()
