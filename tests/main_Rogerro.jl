using ITensors
using Plots
using Measures
using ITensorTDVP
# using TimeEvoMPS
include("../src/evolution.jl")
include("../src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N = 4 # number of sites (NEED TO GO TILL 96 for Rog_results)
    cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.25 # time step (NEED TO BE 0.05 for Rog_results)
    ttotal = 50 # total time of evolution (NEED TO GO TILL 50 for Rog_results)
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
    # maxdim = 1

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=false doesnt conserve the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", N; conserve_qns=false)  

    # Constants for Rogerro's fit (only self-interaction term)
    a_t = 1.224
    b_t = 0
    c_t = 1.62
    
    # Initialize an array of ones for all N sites
    mu = ones(N) # erg
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    n = mu .* fill((Δx)^3/(sqrt(2) * G_F), N)
    
    # Create a B vector which would be same for all N particles 
    # B = [0, 0, -1]
    theta_nu= 0.1 #0.5986 #rad # =34.3 degrees
    B = [sin(2 *theta_nu), 0, -cos(2*theta_nu)]
    B = B / norm(B) 
    # Create arrays ω_a and ω_b
    ω_a = fill(0.2, div(N, 2))
    ω_b = fill(0, div(N, 2))

    # Defining Δω as in Rogerro(2021)
    Δω = (ω_a - ω_b)/2
    
    # Concatenate ω_a and ω_b to form ω
    ω = vcat(ω_a, ω_b)

    ψ = productMPS(s, n -> n <= N/2 ? "Dn" : "Up")
    energy_sign = [i <= N ÷ 2 ? 1 : 1 for i in 1:N]

    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, n, ω, B, N, Δx, ψ, energy_sign, cutoff,ttotal)

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
    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    # @assert abs(t_min - t_p_Rog) <  τ + tolerance 

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)",title = "Running main_Rogerro script", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N$(N)"]) 
    scatter!([t_p_Rog],[prob_surv_array[i_first_local_min]], label= ["t_p_Rog"])
    scatter!([t_min],[prob_surv_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
    # Save the plot as a PDF file
    savefig("Survival probability vs t (Rog)for N$(N).pdf")
end 

@time main()

