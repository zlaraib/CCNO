push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots 
using Measures 
using LinearAlgebra
using DelimitedFiles
using HDF5


# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
    
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    ttotal = 1.94e-4
    tolerance  = 1E-5 # acceptable level of error or deviation from the exact value or solution #variable
    
    params = CCNO.Parameters(
        N_sites = 6, # number of sites, #variable
        cutoff = 1E-14, # specifies a truncation threshold for the SVD in MPS representation #variable
        τ = ttotal/50.0, # time step # sec #variable
        ttotal = ttotal, # total time of evolution #sec #variable
        Δx = 1E-3, # length of the box of interacting neutrinos at a site in cm #variable
        maxdim = 1000, # max bond dimension in MPS truncation
        L=L,
        Δp = L, # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.
        periodic = false,  # true = imposes periodic boundary conditions while false doesn't
        checkpoint_every = 4,
        do_recover = false,
        recover_file = "",
        m1 = 0*CCNO.eV,
        m2 = 0.008596511*CCNO.eV,
        datadir = joinpath(@__DIR__, "datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        shape_name = "none",  # variable.
        geometric_name = "physical",
        theta_nu = π/4, # mixing_angle #rad 
        α = 0,
        save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    )

    # Make an array of 'site' indices and label as s 
    # conserve_qns=false doesnt conserve the total spin quantum number "S"(in z direction) in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  #Fixed
    
    # Initialize an array of zeros for all params.N_sites particles
    mu = zeros(params.N_sites) #Fixed
    
    # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos.
    N = mu.* fill((params.Δx)^3/(sqrt(2) * CCNO.G_F), params.N_sites) 
    
    # Create a B vector which would be same for all N particles 
    B = [sin(2*params.theta_nu), 0, -cos(2*params.theta_nu)] # is equivalent to B = [1, 0, 0] # variable. But only other case that can be tested from this file is B = [0,0,-1] for which theta_nu = π/4.
    B = B / norm(B)
    
    x = fill(rand(), params.N_sites) # variable.
    y = fill(rand(), params.N_sites) # variable.
    z = zeros(params.N_sites) # variable.
    
    # Generate an Nx3 array for p with random values
    p = ones(params.N_sites, 3)*CCNO.MeV # variable, but will need to make sure that p_vector.jl file if statment stays constsnet 
    energy_sign = fill(1, params.N_sites) # all of the sites are neutrinos
    p_mod, p_hat = CCNO.momentum(p)
    
    # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Dn" : "Up") # Fixed to produce consistent results for the test assert conditions 
    
    state = CCNO.SimulationState(ψ=ψ,
                                 s=s,
                                 s0=s,
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

    expected_sz_array = []
    Δm²::Float64 = abs(params.m2^2 - params.m1^2)
    ω::Vector{Float64} = [Δm² / (2 * p_mod[i]) * state.energy_sign[i] for i in 1:params.N_sites] / CCNO.hbar
    println("one cycle time:",2.0*π/ω[1])
    for t in 0.0:params.τ:params.ttotal
        i = 1 # change it according to the corresponding site number in the expect function
        # Compute the expected value based on the derived analytic formula
        expected_sz = -0.5 * cos(ω[i] * t)
        push!(expected_sz_array, expected_sz)
    end
    print(expected_sz_array)
    # Check if every element in Sz_array is less than tolerance away from the corresponding element in expected_sz_array
    # for B vector in x, it checks that the value of Sz at the first spin site oscillates between -0.5 and 0.5 
    # for B vector in -z, it checks that the value of Sz at the firstspin site never oscillates from -0.5 
    #@assert all(abs.(Sz_array .- expected_sz_array) .< tolerance)

  if params.save_plots_flag
    
    x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
    pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

    CCNO.save_plots(τ, params.N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,params.datadir, plotdir, save_plots_flag)
  end 
  if !params.save_plots_flag
    plot(t_array, Sz_array, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script",
    legend = true, size=(700, 600), aspect_ratio=:auto,left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, 
    bottom_margin = 10mm, label = "My_sz_site1") 
    plot!(t_array, expected_sz_array, linestyle=:dot, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script", 
    legendfontsize=8, legend=:topright, label = "Expected_sz from Sakurai_site1") 
    # Save the plot as a PDF file # for jenkins archive 
    savefig("<Sz>_site1 vs t.pdf")
  end
    
    # clean up
    rm(params.datadir, recursive=true)
    rm(params.chkptdir, recursive=true)
end

@time main()
