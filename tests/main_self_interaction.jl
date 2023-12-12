using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics

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
V = L³ 
N_sites = total no.of sites or particles
N =  total no.of neutrinos at all sites
n = total no.density of neutrinos at all sites
n = N / V 

small box (for 'each' site/particle)
total vol = Vᵢ
length of each side = Δx 
Δx³ = Vᵢ
Nᵢ =  total no.of neutrinos at site i 
nᵢ  = total no.density of neutrinos at site i 
Nᵢ = nᵢ  * Vᵢ

Combining (where index i represent a site and runs from 1:N_sites)
∑ᵢ nᵢ = n # total no.density of neutrinos
∑ᵢ  Nᵢ = N # total no.of neutrinos
∑ᵢ Vᵢ = V # total volume of the grid
∑ᵢ Δxᵢ = L # domain size

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
function main(N_sites_eachflavor,τ,ttotal,tolerance,Δm²,maxdim,cutoff,L,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,shape_name,periodic)

    function generate_inputs_file(directory, filename, data)
        filepath = joinpath(directory, filename)
        
        # Open the file in write mode
        open(filepath, "w") do file
            # Write data to the file
            for line in data
                println(file, line)
            end
        end
    end
    
    function extract_initial_conditions()
        # Putting all my initial conditions into a dictionary
        initial_conditions = Dict(
            "N_sites" => N_sites,
            "N_sites_eachflavor" => N_sites_eachflavor,
            "τ" => τ,
            "ttotal" => ttotal,
            "tolerance" => tolerance,
            "Δm²" => Δm²,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "L" => L,
            "n_νₑ" => n_νₑ,
            "n_νₑ̄" => n_νₑ̄,
            "Eνₑ" => Eνₑ,
            "Eνₑ̄" => Eνₑ̄,
            "shape_name" => shape_name,
            "Δp" => Δp,
            "periodic" => periodic
        )

        # Convert the dictionary to a string and write it into the input file
        input_data = []
        for (key, value) in initial_conditions
            push!(input_data, "$key = $value")
        end

        return input_data
    end
    
    # Generate input data
    input_data = extract_initial_conditions()

    V = L^3 # volume of the big box containing all sites/particles
    Δx = L/N_sites # length of the box of interacting neutrinos at a site in cm  #variable

    # Create an array of dimension N and fill it half with values of sites containing all electron neutrinos 
    # and other half with sites containing electron anti-neutrino. 
    N_νₑ  = n_νₑ * V 
    N_1 = fill(N_νₑ / (N_sites ÷ 2), N_sites ÷ 2)
    N_νₑ̄  = n_νₑ̄ * V 
    N_2 = fill(N_νₑ̄/ (N_sites ÷ 2), N_sites ÷ 2)
    N = vcat(N_1, N_2) # This is the total number of neutrinos. 

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

    #generate a momentum array that depicts the energy of neutrinos and anti-neutrinos in opposing beams
    function generate_p_array(N_sites)                                                                                                                                                                                   
        half_N_sites = div(N_sites, 2)
        return [fill(Eνₑ, half_N_sites); fill(Eνₑ̄, half_N_sites)]
    end

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites), generate_p_array(N_sites), generate_p_array(N_sites))
    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= N_sites ÷ 2 ? 1 : -1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry
    
    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.
    s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    #ψ = productMPS(s, N -> N <= N_sites/2 ? "Up" : "Dn") # Fixed to produce consistent results for the test assert conditions 
    ψ = productMPS(s, N_sites -> "Up") # inital state all in spin up direction i.e. all electron flavor
    # Perturb the state via one-body Hamiltonian
    ψ₀= evolve_perturbation(s, τ, B, N_sites, ψ, cutoff, maxdim, ttotal)
    
    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "..","misc","datafiles","FFI", "par_"*string(N_sites), "tt_"*string(ttotal))

    # Call the function to generate the inputs file in the specified directory
    generate_inputs_file(datadir, "inputs.txt", input_data)

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array= evolve(s, τ, N, B,L, N_sites, 
                    Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, datadir, ttotal,periodic)

    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","FFI", "par_"*string(N_sites), "tt_"*string(ttotal))
    
    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)

    # Call the function to generate the inputs file in the specified directory
    generate_inputs_file(plotdir, "inputs.txt", input_data)
    
    # Plotting ρ_μμ vs t
    plot(0.0:τ:τ*(length(ρ_μμ_array)-1), ρ_μμ_array, xlabel = "t", ylabel = "<ρ_μμ>", legend = false, 
    left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
    # Save the plot as a PDF file
    savefig(joinpath(plotdir, "<ρ_μμ>_vs_t_self-interactions_w_geo+shape_MF_FFI.pdf"))

    # Plotting ρ_ee vs t
    plot(0.0:τ:τ*(length(ρₑₑ_array)-1), ρₑₑ_array, xlabel = "t", ylabel = "<ρₑₑ>", legend = false, 
    left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
    # Save the plot as a PDF file
    savefig(joinpath(plotdir, "<ρₑₑ>_vs_t_self-interactions_w_geo+shape_MF_FFI.pdf"))

   # Plotting P_surv vs t
   plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)",
   legend = false, left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm)
   savefig(joinpath(plotdir,"Survival probability vs t (only self-interaction term)_FFI.pdf"))

    # Plotting Sz vs t
    plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false,
        left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm,margin= 10mm,
        ylims = (minimum(Sz_array), maximum(Sz_array) + 0.1 * abs(maximum(Sz_array) - minimum(Sz_array))) )
    # Save the plot as a PDF file
    savefig(joinpath(plotdir, "<Sz>_vs_t_self-interactions_w_geo+shape_MF_FFI.pdf"))

    # Plotting Sy vs t
    plot(0.0:τ:τ*(length(Sy_array)-1), Sy_array, xlabel = "t", ylabel = "<Sy>", legend = false,
    left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sy> vs t (self-interactions w geo+shape)_MF_FFI.pdf"))

    # Plotting Sx vs t
    plot(0.0:τ:τ*(length(Sx_array)-1), Sx_array, xlabel = "t", ylabel = "<Sx>", legend = false,
    left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
    #Save the plot as a PDF file
    savefig(joinpath(plotdir,"<Sx> vs t (self-interactions w geo+shape)_MF_FFI.pdf"))

    # Plotting particles positional evolution
    plot(title="Position Evolution for $N_sites particles", xlabel= "Position (x)",ylabel="Time(s)")
    for site in 1:N_sites
        site_positions = [(x_values[t][site]) for t in 1:length(x_values)]
        plot!(site_positions, 0.0:τ:ttotal, label="Site $site",
        left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm)
    end
    savefig(joinpath(plotdir,"Particles position(x) evolution.pdf"))
    
    # Plotting particles momentum evolution
    plot(title="Particle Momentum Evolution", xlabel= "Momentum in x direction(pₓ)",ylabel="Time")
    for site in 1:N_sites
        site_momentum = [(pₓ_values[t][site]) for t in 1:length(pₓ_values)]
        plot!(site_momentum, 0.0:τ:ttotal, label="Site $site",left_margin = 25mm, right_margin = 5mm, 
        top_margin = 5mm, bottom_margin = 10mm)
    end
    savefig(joinpath(plotdir,"Particles momentum(pₓ) evolution.pdf"))

end 


