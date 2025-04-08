push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    
    params = CCNO.Parameters(
        N_sites = 4, # number of sites (NEED TO GO TILL 96 for Rog_results) # variable.
        cutoff = 1E-14, # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) # variable.
        τ = 0.05*CCNO.hbar, # time step #sec #fixed for Rogerros result
        ttotal = 5*CCNO.hbar, # total time of evolution #sec # variable.
        tolerance = 5E-1, # acceptable level of error or deviation from the exact value or solution # variable.
        Δx = 1E-3, # length of the box of interacting neutrinos at a site in cm  # variable.
        L=L,
        m1 = 0.0,
        m2 = 0.0,
        maxdim = 1000, # max bond dimension in MPS truncation
        Δp = L, # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.  
        periodic = false,  # true = imposes periodic boundary conditions while false doesn't
        save_plots_flag = false, # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
        checkpoint_every = 4,
        do_recover = false,
        recover_file = "" ,
        theta_nu = 0, # mixing_angle #rad
        shape_name = "none",  # Change this to the desired shape name # variable.
        geometric_name = "none",
        datadir = joinpath(@__DIR__, "datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        α = 0
    )

    

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  #fixed

    # Fixed Constants for Rogerro's fit (only self-interaction term)
    a_t = 0
    b_t = 2.105
    c_t = 0
    
    # Initialize an array of ones for all N sites
    mu = ones(params.N_sites) # erg #fixed
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the total number of neutrinos. 
    N = mu .* fill(((params.Δx)^3 )/(√2 * CCNO.G_F * params.N_sites), params.N_sites)

    x = fill(rand(), params.N_sites) # variable.
    y = fill(rand(), params.N_sites) # variable.
    z = fill(rand(), params.N_sites) # variable.

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Dn" : "Up") #fixed for Rog case

    # array p with N rows and 3 columns, all initialized to 0.0 with colums representing components and rows representing sites
    p = ones(params.N_sites, 3) #fixed for Rogerro's case
    energy_sign = fill(1, params.N_sites) # all of the sites are neutrinos

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

    # This function scans through the array, compares each element with its neighbors, 
    # and returns the index of the first local minimum it encounters. 
    # If no local minimum is found, it returns -1 to indicate that.
    function find_first_local_minima_index(arr)
        n = length(arr)
        for i in 2:(n-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
    
    tmin_ifirstlocalmin_file = joinpath(params.datadir, "tmin_ifirstlocalmin.dat")
    # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
    i_first_local_min = find_first_local_minima_index(-Sz_array)
    
    # Writing if_else statement to communicate if local minima (not) found
    t_min = nothing
    if i_first_local_min != -1
        println("Index of the first local minimum: ", i_first_local_min)
        t_min = params.τ * i_first_local_min - params.τ
    else
        println("No local minimum found in the array.")
    end

    println("Corresponding time of first minimum index= ", t_min)

    # Save t_min and i_first_local_min to file
    writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
    println("Saved t_min to file: ", tmin_ifirstlocalmin_file)


    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(params.N_sites) + b_t * sqrt(params.N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

    if params.save_plots_flag
        
        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

        CCNO.save_plots(τ, params.N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,params.datadir, plotdir, save_plots_flag)
    end 

    if !params.save_plots_flag
        # Plotting P_surv vs t
        plot(t_array/CCNO.hbar, -Sz_array, xlabel = "t", ylabel = "-Sz(t)",
        title = "Running main_self_interaction_Rog script", legend = true, size=(800, 600), aspect_ratio=:auto,margin= 10mm, 
        label= ["My_plot_for_N_sites$(params.N_sites)"]) 
        scatter!([t_p_Rog],[-Sz_array[i_first_local_min]], label= ["t_p_Rog"])
        scatter!([t_min/CCNO.hbar],[-Sz_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
        # Save the plot as a PDF file # for jenkins archive 
        savefig("main_self_interaction_Rog.pdf")
    end

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    @assert abs(t_min/CCNO.hbar - t_p_Rog) < params.τ/CCNO.hbar + params.tolerance
end 

@time main()


