push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5


# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    Δm²= 0.2 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.

    params = CCNO.Parameters(
        N_sites = 4, # number of sites # make it 24 to produce Rog results. #reduced for unit test passing
        cutoff = 1E-10, # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
        τ = 0.25*CCNO.hbar, # time step (NEED TO BE 0.05 for Rog_ main_text results)
        ttotal = 50*CCNO.hbar, # total time of evolution
        tolerance  = 5E-1, # acceptable level of error or deviation from the exact value or solution
        m1 = 0,
        m2 = sqrt(Δm²),
        L=L,
        theta_nu= 0.1, #rad # =34.3 degrees
        Δx = 1E-3, # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
        Δp = L, # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent. 
        maxdim = 1, #bond dimension
        periodic = true,  # true = imposes periodic boundary conditions while false doesn't
        checkpoint_every = 100,
        do_recover = false,
        recover_file = "",
        geometric_name = "none",
        shape_name = "none",  # Change this to the desired shape name # variable.
        datadir = joinpath(@__DIR__, "datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        save_plots_flag = false, # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
        α = 0
    )

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N_sites.
    # conserve_qns=false doesnt conserve the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  
    
    # Initialize an array of ones for all N_sites sites
    mu = ones(params.N_sites) # erg
    
    # Create an array of dimension N_sites and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    N = mu .* fill(((params.Δx)^3 )/(√2 * CCNO.G_F * params.N_sites), params.N_sites)

    x = fill(rand(), params.N_sites) # variable.
    y = fill(rand(), params.N_sites) # variable.
    z = fill(rand(), params.N_sites) # variable.

    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Dn" : "Up")

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(CCNO.generate_p_array(params.N_sites),fill(0, params.N_sites), fill(0, params.N_sites))
    energy_sign = [i <= params.N_sites ÷ 2 ? 1 : 1 for i in 1:params.N_sites] # all of the sites are neutrinos
    
    state = CCNO.SimulationState(ψ=ψ,
                                 s=s,
                                 p=p,
                                 energy_sign = energy_sign,
                                 N=N,
                                 xyz = hcat(x,y,z))

    #extract output for the survival probability values at each timestep
    CCNO.evolve(params, state)

    # Read the data files
    t_Sz_tot = readdlm(joinpath(params.datadir, "t_<Sz>.dat"))
    t_Sy_tot = readdlm(joinpath(params.datadir, "t_<Sy>.dat"))
    t_Sx_tot = readdlm(joinpath(params.datadir, "t_<Sx>.dat"))
    t_xsiteval = readdlm(joinpath(params.datadir, "t_xsiteval.dat"))
    t_pxsiteval = readdlm(joinpath(params.datadir, "t_pxsiteval.dat"))
    t_ρₑₑ_tot = readdlm(joinpath(params.datadir, "t_ρₑₑ.dat"))
    t_ρ_μμ_tot = readdlm(joinpath(params.datadir, "t_ρ_μμ.dat"))
    t_ρₑμ_tot = readdlm(joinpath(params.datadir, "t_ρₑμ.dat"))
    
    # Extract time array and corresponding values for plotting
    t_array = t_Sz_tot[:, 1]  
    
    #Extract the array for first site only
    Sz_array = t_Sz_tot[:,2]
    Sy_array = t_Sy_tot[:,2]
    Sx_array= t_Sx_tot[:,2] 
    ρₑₑ_array = t_ρₑₑ_tot[:, 2]
    ρ_μμ_array = t_ρ_μμ_tot[:, 2]
    ρₑμ_array = t_ρₑμ_tot[:, 2]
        
    if params.save_plots_flag

        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

        CCNO.save_plots(τ, params.N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,datadir, plotdir, save_plots_flag)
    end
    if !params.save_plots_flag 
        # Plotting P_surv vs t 
        plot(t_array/CCNO.hbar, Sz_array, xlabel = "t", ylabel = "S_z(t)",title = "Running main_Bipolar_Rog script \n for N_sites$(params.N_sites) with maxdim=1 and cutoff for τ$(params.τ)", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N$(params.N_sites)"]) 
        # Save the plot as a PDF file # for jenkins archive 
        savefig("S_z vs t (Rog_bipolar)for N_sites$(params.N_sites) with maxdim=$(params.maxdim) and cutoff for τ$(params.τ).pdf")
    end
end 

@time main()
