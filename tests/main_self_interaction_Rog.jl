using ITensors
using Plots
using Measures
using LinearAlgebra
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
include(src_dir * "Utilities/save_plots.jl")
include(src_dir * "src/chkpt_hdf5.jl")
include(src_dir * "Utilities/save_datafiles.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N_sites = 4 # number of sites (NEED TO GO TILL 96 for Rog_results) # variable.
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) # variable.
    τ = 0.05 # time step #sec #fixed for Rogerros result
    ttotal = 5 # total time of evolution #sec # variable.
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution # variable.
    Δx = 1E-3 # length of the box of interacting neutrinos at a site in cm  # variable.
    Δm²= 0.0 # erg^2 # Fixed for rog only self-int test case. Please dont play with it. 
    maxdim = 1000 # max bond dimension in MPS truncation
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.  
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test
    periodic = false  # true = imposes periodic boundary conditions while false doesn't
    
    checkpoint_every = 4
    do_recover = false
    recover_file = "" 

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N_sites; conserve_qns=false)  #fixed

    # Fixed Constants for Rogerro's fit (only self-interaction term)
    a_t = 0
    b_t = 2.105
    c_t = 0
    
    # Initialize an array of ones for all N sites
    mu = ones(N_sites) # erg #fixed
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the total number of neutrinos. 
    N = mu .* fill(((Δx)^3 )/(√2 * G_F * N_sites), N_sites)

    # Create a B vector which would be same for all N particles 
    theta_nu = 0 # mixing_angle #rad 
    B = [sin(2*theta_nu), 0, -cos(2*theta_nu)] # is equivalent to B = [0, 0, -1] # fixed for Rogerro's case
    B = B / norm(B)

    x = fill(rand(), N_sites) # variable.
    y = fill(rand(), N_sites) # variable.
    z = fill(rand(), N_sites) # variable.

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up") #fixed for Rog case

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "none"  # Change this to the desired shape name # variable.

    # array p with N rows and 3 columns, all initialized to 0.0 with colums representing components and rows representing sites
    p = zeros(N_sites, 3) #fixed for Rogerro's case
    energy_sign = fill(1, N_sites) # all of the sites are neutrinos

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "datafiles")
    chkptdir = joinpath(@__DIR__, "checkpoints")

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover = evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,chkptdir, checkpoint_every,  do_recover, recover_file, save_data , periodic)

    # extract the prob_surv on the first site 
    prob_surv_array_site1= [row[1] for row in prob_surv_array]
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
    
    tmin_ifirstlocalmin_file = joinpath(datadir, "tmin_ifirstlocalmin.dat")
    if !do_recover 
        # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
        i_first_local_min = find_first_local_minima_index(prob_surv_array_site1)
        
        # Writing if_else statement to communicate if local minima (not) found
        if i_first_local_min != -1
            println("Index of the first local minimum: ", i_first_local_min)
            t_min = τ * i_first_local_min - τ
            println(t_min)
        else
            t_min = nothing
            println(t_min)
            println("No local minimum found in the array.")
        end
        println("Corresponding time of first minimum index= ", t_min)

        if save_data
            # Save t_min and i_first_local_min to file
            writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
            println("Saved t_min to file: ", tmin_ifirstlocalmin_file)
        else t_min = t_min
        end
    elseif do_recover 

        # Read the saved t_min and i_first_local_min from the file
        tmin_data = readdlm(tmin_ifirstlocalmin_file)

        # Extracting t_min and i_first_local_min from the read data
        t_min_value = tmin_data[1]  # First value
        i_first_local_min_value = tmin_data[2]  # Second value

        # Convert `t_min` and `i_first_local_min` to appropriate types if needed
        t_min = t_min_value isa AbstractString ? tryparse(Float64, t_min_value) : t_min_value
        i_first_local_min = i_first_local_min_value isa AbstractString ? tryparse(Int, i_first_local_min_value) : i_first_local_min_value

        println("Recover t_min=",t_min)
        println("Recover i_first_local_min=",i_first_local_min)

        if t_min !== nothing && i_first_local_min != -1 
            # t_min = readdlm(tmin_ifirstlocalmin_file)[1]  # assign the previous t_min
            println("Recovered first local t_min from previous run: ", t_min)
            println("Recovered index of first local t_min from previous run: ", i_first_local_min)
        end
        if t_min === nothing && i_first_local_min == -1 
            i_first_local_min = find_first_local_minima_index(prob_surv_array_site1)
            t_min = (τ * i_first_local_min) - τ + t_recover
            println("Recalculated t_min=",t_min)
            println("Recalculated i_first_local_min=",i_first_local_min)
            if i_first_local_min != -1
                println("Index of the first local minimum: ", i_first_local_min)
                println("Final t_min=",t_min)
            else
                t_min = nothing
                i_first_local_min = -1
                println(t_min, i_first_local_min)
                println("No local minimum found in the array.")
            end
            
        end
        # Save t_min and i_first_local_min to file (again, updated if required)
        writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
        println("Saved t_min to file: ", tmin_ifirstlocalmin_file)
    end 

    println("Corresponding time of first minimum index= ", t_min)


    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(N_sites) + b_t * sqrt(N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

   
    if save_data
        t_min_value = readdlm(tmin_ifirstlocalmin_file)[1]  # read the previous t_min
        t_min = t_min_value isa AbstractString ? tryparse(Float64, t_min_value) : t_min_value
        println(t_min)
    end

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    if t_min !== nothing
        @assert abs(t_min - t_p_Rog) < τ + tolerance
    elseif t_min === nothing
        if do_recover
            t_min = (τ * i_first_local_min) - τ + t_recover
        else
            t_min = τ * i_first_local_min - τ
        end
        @assert abs(t_min - t_p_Rog) <  τ + tolerance
    end
    if save_plots_flag
        # Specify the relative directory path
        plotdir = joinpath(@__DIR__, "plots")
        
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

        save_plots(τ, N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,datadir, plotdir, save_plots_flag)
    end 

    if !save_plots_flag
        # Plotting P_surv vs t
        plot(t_array, prob_surv_array_site1, xlabel = "t", ylabel = "Survival Probabillity p(t)",
        title = "Running main_self_interaction_Rog script", legend = true, size=(800, 600), aspect_ratio=:auto,margin= 10mm, 
        label= ["My_plot_for_N_sites$(N_sites)"]) 
        scatter!([t_p_Rog],[prob_surv_array_site1[i_first_local_min]], label= ["t_p_Rog"])
        scatter!([t_min],[prob_surv_array_site1[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
        # Save the plot as a PDF file # for jenkins archive 
        savefig("Survival probability_site1 vs t for N_sites$(N_sites).pdf")
    end

end 

@time main()


