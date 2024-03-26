using ITensors
using Plots
using Measures
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")

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
   
    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N_sites; conserve_qns=true)  #fixed

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
    datadir = joinpath(@__DIR__, "..","misc","datafiles","Rog_self_int", "par_"*string(N_sites), "tt_"*string(ttotal))

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, Im_Ω = evolve(s, τ, N, B,L, N_sites, 
                    Δx,Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,periodic)


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
    
    # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
    i_first_local_min = find_first_local_minima_index(prob_surv_array)
    
    # Writing if_else statement to communicate if local minima (not) found
    if i_first_local_min != -1
        println("Index of the first local minimum: ", i_first_local_min)
    else
        println("No local minimum found in the array.")
    end

    # Time at which the first mimimum survival probability is reached
    t_min = τ * i_first_local_min - τ
    println("Corresponding time of first minimum index= ", t_min)

    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(N_sites) + b_t * sqrt(N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    @assert abs(t_min - t_p_Rog) <  τ + tolerance 

    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","Rog_self_int", "par_"*string(N_sites), "tt_"*string(ttotal))
    
    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)",
    title = "Running main_self_interaction_Rog script", legend = true, size=(800, 600), aspect_ratio=:auto,margin= 10mm, 
    label= ["My_plot_for_N_sites$(N_sites)"]) 
    scatter!([t_p_Rog],[prob_surv_array[i_first_local_min]], label= ["t_p_Rog"])
    scatter!([t_min],[prob_surv_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
    # Save the plot as a PDF file
    savefig(joinpath(plotdir,"Survival probability vs t (only self-interaction term plot)_Rog for N_sites$(N_sites).pdf"))

    plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false, size=(800, 600), 
    left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, aspect_ratio=:auto,margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sz> vs t (Rog_self-int).pdf"))

    plot(0.0:τ:τ*(length(Sy_array)-1), Sy_array, xlabel = "t", ylabel = "<Sy>", legend = false, size=(800, 600), 
    left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, aspect_ratio=:auto,margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sy> vs t (Rog_self_int).pdf"))

    plot(0.0:τ:τ*(length(Sx_array)-1), Sx_array, xlabel = "t", ylabel = "<Sx>", legend = false, size=(800, 600), 
    left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, aspect_ratio=:auto,margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sx> vs t (Rog_self_int).pdf"))
end 

@time main()

