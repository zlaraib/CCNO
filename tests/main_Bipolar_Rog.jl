using ITensors
using Plots
using Measures
using ITensorTDVP
using DelimitedFiles
# using TimeEvoMPS
include("../src/evolution.jl")
include("../src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    N_sites = 4 # number of sites # make it 24 to produce Rog results. #reduced for unit test passing
    cutoff = 1E-10 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    τ = 0.25 # time step (NEED TO BE 0.05 for Rog_ main_text results)
    ttotal = 50 # total time of evolution
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
    Δm²= 0.2 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.
    maxdim = 1 #bond dimension
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

    # Create a B vector which would be same for all N_sites particles 
    theta_nu= 0.1 #rad # =34.3 degrees
    B = [sin(2 *theta_nu), 0, -cos(2*theta_nu)]
    B = B / norm(B) 

    x = fill(rand(), N_sites) # variable.
    y = fill(rand(), N_sites) # variable.
    z = fill(rand(), N_sites) # variable.

    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")

    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "none"  # Change this to the desired shape name # variable.

    #generate a momentum array that depicts the energy of neutrinos and anti-neutrinos in opposing beams
    function generate_p_array(N_sites)                                                                                                                                                                                   
        half_N_sites = div(N_sites, 2)
        return [fill(1, half_N_sites); fill(1, half_N_sites)]
    end

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites),fill(0, N_sites), fill(0, N_sites))
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites] # all of the sites are neutrinos

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "..","misc","datafiles","Rog_bipolar", "par_"*string(N_sites), "tt_"*string(ttotal), "τ_"*string(τ))
    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, Im_Ω = evolve(s, τ, N, B,L, N_sites, 
                    Δx,Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,periodic)
    
    # Defining Δω as in Rogerro(2021)
    Δω = vcat((ω_a - ω_b)/2, (ω_a - ω_b)/2)
    @assert all(Δω./mu .== 0.1)

    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","Rog_bipolar", "par_"*string(N_sites), "tt_"*string(ttotal), "τ_"*string(τ))

    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)

    # Plotting P_surv vs t
    plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probability p(t)",title = "Running main_Bipolar_Rog script \n for N_sites$(N_sites) with maxdim=1 and cutoff for τ$(τ)", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N$(N_sites)"]) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir,"Survival probability vs t (Rog_bipolar)for N_sites$(N_sites) with maxdim=1 and cutoff for τ$(τ).pdf"))
end 

@time main()