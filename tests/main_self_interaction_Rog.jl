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
    τ = 0.05 # time step #fixed for Rogerros result
    ttotal = 5 # total time of evolution (NEED TO GO TILL 50 for Rog_results) # variable.
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution # variable.
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm  # variable.
    Δp = rand() # shape function width # variable.
    del_m2 = 0.0 # Fixed for rog test case. Please dont play with it. 


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
    N = mu .* fill((Δx)^3/(sqrt(2) * G_F), N_sites)
    
    # Create a B vector which would be same for all N particles 
    B = [0, 0, -1] # fixed for Rogerro's case

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

    #extract output for the survival probability values at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, N, B, N_sites, Δx,del_m2, p, x, Δp, ψ, shape_name, energy_sign, cutoff, tolerance, ttotal)

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

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)",title = "Running main_self_interaction_Rog script", legend = true, size=(800, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N_sites$(N_sites)"]) 
    scatter!([t_p_Rog],[prob_surv_array[i_first_local_min]], label= ["t_p_Rog"])
    scatter!([t_min],[prob_surv_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
    # Save the plot as a PDF file
    savefig("Survival probability vs t (only self-interaction term plot)_Rog.pdf")
end 

@time main()

