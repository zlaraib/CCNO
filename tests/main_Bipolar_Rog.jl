using ITensors
using Plots
using Measures
using DelimitedFiles
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
include(src_dir * "src/chkpt_hdf5.jl")
include(src_dir * "Utilities/save_plots.jl")
include(src_dir * "Initializations/initial_cond.jl")
include(src_dir * "Utilities/save_datafiles.jl")

# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N_sites = 4 # number of sites # make it 24 to produce Rog results. #reduced for unit test passing
    N_sites = 4 # number of sites # make it 24 to produce Rog results. #reduced for unit test passing
    cutoff = 1E-10 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.25 # time step (NEED TO BE 0.05 for Rog_ main_text results)
    ttotal = 50 # total time of evolution
    τ = 0.25 # time step (NEED TO BE 0.05 for Rog_ main_text results)
    ttotal = 50 # total time of evolution
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
    Δm²= 0.2 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.
    Δm²= 0.2 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.
    maxdim = 1 #bond dimension
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent. 
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test
    periodic = true  # true = imposes periodic boundary conditions while false doesn't
   
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent. 
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test
    periodic = true  # true = imposes periodic boundary conditions while false doesn't
   
    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N_sites.
    # conserve_qns=false doesnt conserve the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N_sites; conserve_qns=false)  
    
    # Initialize an array of ones for all N_sites sites
    mu = ones(N_sites) # erg
    
    # Create an array of dimension N_sites and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    N = mu .* fill(((Δx)^3 )/(√2 * G_F * N_sites), N_sites)

    N = mu .* fill(((Δx)^3 )/(√2 * G_F * N_sites), N_sites)

    # Create a B vector which would be same for all N_sites particles 
    theta_nu= 0.1 #rad # =34.3 degrees
    theta_nu= 0.1 #rad # =34.3 degrees
    B = [sin(2 *theta_nu), 0, -cos(2*theta_nu)]
    B = B / norm(B) 

    checkpoint_every = 4
    do_recover = false
    recover_type = "auto" 
    recover_iteration = 80 # change it to the iteration you want to recover from, for manual iteration. Currently auto recovery already recovers from last iteration (i.e. recover_iteration = -1 for auto recovery). 
    
    x = fill(rand(), N_sites) # variable.
    y = fill(rand(), N_sites) # variable.
    z = fill(rand(), N_sites) # variable.

    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "none"  # Change this to the desired shape name # variable.

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites),fill(0, N_sites), fill(0, N_sites))
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites] # all of the sites are neutrinos

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "none"  # Change this to the desired shape name # variable.

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites),fill(0, N_sites), fill(0, N_sites))
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites] # all of the sites are neutrinos

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "datafiles","Roggero", "par_"*string(N_sites), "τ_"*string(τ))
    chkptdir = joinpath(@__DIR__, "checkpoints","Roggero", "par_"*string(N_sites), "τ_"*string(τ))
    
    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover = evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,chkptdir, checkpoint_every,  do_recover, recover_type, recover_iteration, save_data , periodic)

    # extract the prob_surv on the first site 
    prob_surv_array_site1= [row[1] for row in prob_surv_array]
    # Defining Δω as in Rogerro(2021)
    Δω = vcat((ω_a - ω_b)/2, (ω_a - ω_b)/2)
    @assert all(Δω./mu .== 0.1)

    if save_plots_flag
        # Specify the relative directory path
        plotdir = joinpath(@__DIR__, "plots","Roggero", "par_"*string(N_sites), "τ_"*string(τ))

        # Read the data files
        t_Sz_tot = readdlm(joinpath(datadir, "t_<Sz>.dat"))
        t_Sy_tot = readdlm(joinpath(datadir, "t_<Sy>.dat"))
        t_Sx_tot = readdlm(joinpath(datadir, "t_<Sx>.dat"))
        t_probsurv_tot = readdlm(joinpath(datadir, "t_probsurv.dat"))
        t_xsiteval = readdlm(joinpath(datadir, "t_xsiteval.dat"))
        t_pxsiteval = readdlm(joinpath(datadir, "t_pxsiteval.dat"))
        t_ρₑₑ_tot = readdlm(joinpath(datadir, "t_ρₑₑ.dat"))
        t_ρ_μμ_tot = readdlm(joinpath(datadir, "t_ρ_μμ.dat"))
        t_ρₑμ_tot = readdlm(joinpath(datadir, "t_ρₑμ.dat"))

        # Extract time array and corresponding values for plotting
        t_array = t_Sz_tot[:, 1]  

        #Extract the array for first site only
        Sz_array = t_Sz_tot[:,2]
        Sy_array = t_Sy_tot[:,2]
        Sx_array= t_Sx_tot[:,2] 
        prob_surv_array = t_probsurv_tot[:, 2]
        ρₑₑ_array = t_ρₑₑ_tot[:, 2]
        ρ_μμ_array = t_ρ_μμ_tot[:, 2]
        ρₑμ_array = t_ρₑμ_tot[:, 2]
        
        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

        # # Parsing arrays containing strings like "[1.0,", into a numeric array suitable for plotting
        Sz_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sz_array]
        Sy_array = [parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sy_array]
        Sx_array =[parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sx_array]
        prob_surv_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in prob_surv_array]
        ρₑₑ_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρₑₑ_array]
        ρ_μμ_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρ_μμ_array]
        ρₑμ_array =[parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρₑμ_array]

        save_plots(τ, N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,datadir, plotdir, save_plots_flag)
    end
    if !save_plots_flag 
        # Plotting P_surv vs t 
        plot(t_array, prob_surv_array_site1, xlabel = "t", ylabel = "Survival Probability p(t)",title = "Running main_Bipolar_Rog script \n for N_sites$(N_sites) with maxdim=1 and cutoff for τ$(τ)", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N$(N_sites)"]) 
        # Save the plot as a PDF file # for jenkins archive 
        savefig("Survival probability_site1 vs t (Rog_bipolar)for N_sites$(N_sites) with maxdim=1 and cutoff for τ$(τ).pdf")
    end
end 

@time main()