using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics
using Random
using HDF5

"""
For github unit tests runs: 
src_dir = ../
save_data and save_plots_flag should be false to run test files. 
"""

src_dir = "../"
save_data = false  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs


"""
For science runs: 
src_dir = /home/zohalaraib/Oscillatrino/ # should be changed to users PATH
save_data and save_plots_flag should be true to run test files. 

"""
# src_dir= "/home/zohalaraib/Oscillatrino/" # should be changed to users PATH
# save_data = true  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
# save_plots_flag = true # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    

include(src_dir * "src/evolution.jl")
include(src_dir * "src/constants.jl")
include(src_dir * "src/shape_func.jl")
include(src_dir * "src/momentum.jl")
include(src_dir * "src/perturb.jl")
include(src_dir * "src/chkpt_hdf5.jl")
include(src_dir * "Initializations/initial_cond.jl")
include(src_dir * "Utilities/gen_input_file.jl")
include(src_dir * "Utilities/save_plots.jl")

function main()
    """ Richers(2021) Test 3 initial conditions for many-body dynamics: """
    N_sites_eachflavor= 1 # total sites/particles that evenly spaced "for each (electron) flavor" 
    N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
    τ = 1.666e-7 # time step from Richers test # sec # variable
    ttotal = 1.666e-2 # total time of evolution # sec #variable
    tolerance  = 5E-3 # acceptable level of error or deviation from the exact value or solution #variable
    m1 = -0.008596511*eV #ergs #1st mass eigenstate of neutrino in Richers(2021)
    m2 = 0*eV   #ergs #2nd mass eigenstate of neutrino in Richers(2021)
    Δm² = (m2^2-m1^2) # mass square difference # (erg^2)
    maxdim = 2 # max bond dimension in MPS truncation
    cutoff = 1e-100 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    L = 1e7 # cm # domain size # (aka big box length)
    n_νₑ =  2.92e24 # cm^-3 # number density of electron flavor neutrino
    n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
    Eνₑ =  50*MeV # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
    Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
    Δx = L # length of the box of interacting neutrinos at a site in cm

    theta_nu = 1.74532925E-8  #1e-6 degrees # mixing_angle # = 1.74532925E-8 radians 
    B = [-sin(2 * theta_nu), 0, cos(2 * theta_nu)]  # actual b vector that activates the vacuum oscillation term in Hamiltonian
    B = B / norm(B) 
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "none"  # Change this to the desired shape name #variable 
    Δp = L # width of shape function  # cm #variable
    t1 = 0.0084003052 #choose initial time for growth rate calculation
    t2 = 0.011700318 #choose final time for growth rate calculation
    periodic = true  # true = imposes periodic boundary conditions while false doesn't
    analytic_growth_rate=  (abs(m2^2 - m1^2)/ (2*hbar* Eνₑ)) # analytic growth rate 

    checkpoint_every = 4
    do_recover = false
    recover_type = "auto" 
    recover_iteration = 80 # change it to the iteration you want to recover from, for manual iteration. Currently auto recovery already recovers from last iteration (i.e. recover_iteration = -1 for auto recovery). 

    x = generate_x_array(N_sites_eachflavor, L)
    y = generate_x_array(N_sites_eachflavor, L)
    z = generate_x_array(N_sites_eachflavor, L)


    # p matrix with numbers generated from the p_array for all components (x, y, z) #sherood has 
    p = hcat(generate_px_array(N_sites, Eνₑ, Eνₑ̄), generate_py_array(N_sites), generate_pz_array(N_sites))

    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= N_sites ÷ 2 ? -1 : 1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

    s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

    # Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up in Richers notation which is equivalently half spin up and half chain spin down in my TN notation)
    ψ₀ = productMPS(s, n -> n <= N_sites/2 ? "Up" : "Dn")

    N = Neutrino_number(s, τ, B,L, N_sites, N_sites_eachflavor, tolerance,
                    n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, ttotal,periodic)

    # Specify the relative directory path
    datadir = joinpath(@__DIR__,"datafiles","FFI_MB", "par_"*string(N_sites))
    chkptdir = joinpath(@__DIR__, "checkpoints","FFI_MB", "par_"*string(N_sites))

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, Im_Ω, t_recover = evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ₀, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,chkptdir, checkpoint_every,  do_recover, recover_type, recover_iteration, save_data , periodic)

    # @assert abs((Im_Ω - analytic_growth_rate)/  analytic_growth_rate) < tolerance #this assert condition doesnt pass with 'any' non-1 maxdim

    if save_data
        # Generate input data
        input_data = extract_initial_conditions(N_sites,N_sites_eachflavor,τ, ttotal,tolerance,
        Δm², maxdim, cutoff, p,ψ₀,L, Δx,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,B, N, shape_name,Δp,periodic)
        # Call the function to generate the inputs file in the specified directory
        generate_inputs_file(datadir, "inputs.txt", input_data)
    end


    if save_plots_flag 
        # Specify the relative directory path
        plotdir = joinpath(@__DIR__, "plots","FFI_MB", "par_"*string(N_sites))
            
        save_plots(τ, N_sites,L,tolerance, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,datadir, plotdir, save_plots_flag)
        
        # Call the function to generate the inputs file in the specified directory
        generate_inputs_file(plotdir, "inputs.txt", input_data)
    end
    if !save_plots_flag 
        # Plotting ρₑμ vs t # for jenkins file 
        plot(0.0:τ:τ*(length(ρₑμ_array)-1), ρₑμ_array, xlabel = "t", ylabel = "<ρₑμ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig( "Homo_MB_<ρₑμ>_vs_t for $N_sites particles.pdf")
    end
end
@time main()