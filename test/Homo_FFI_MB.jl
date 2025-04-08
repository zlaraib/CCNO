push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics
using Random
using HDF5


""" Richers(2021) Test 3 initial conditions: """
function main()
    N_sites_eachflavor= 1 # total sites/particles that evenly spaced "for each (electron) flavor" 
    L = 1e7 # cm # domain size # (aka big box length)
    Δx = L # length of the box of interacting neutrinos at a site in cm

    params = CCNO.Parameters(
        save_plots_flag = false,
        N_sites = 2* (N_sites_eachflavor),
        τ = 1.666e-5,
        Δp = L,
        Δx=Δx,
        L=L,
        ttotal = 1.666e-2,
        tolerance  = 0.2,
        m1 = -0.008596511*CCNO.eV,
        m2 = 0*CCNO.eV,
        maxdim = 2,
        cutoff = 1e-100,
        theta_nu = 1.74532925E-8,  #1e-6 degrees # = 1.74532925E-8 radians 
        shape_name = "none",
        periodic = true,
        checkpoint_every = 1000,
        do_recover = false,
        recover_file = "",
        datadir = joinpath(@__DIR__,"datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        α = 0
    )

    Δm² = (params.m2^2-params.m1^2) # mass square difference # (erg^2)
    n_νₑ =  2.92e24 # cm^-3 # number density of electron flavor neutrino
    n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
    Eνₑ =  50*CCNO.MeV # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
    Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign

    B = [-sin(2 * params.theta_nu), 0, cos(2 * params.theta_nu)]  # actual b vector that activates the vacuum oscillation term in Hamiltonian
    B = B / norm(B) 
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    t1 = 0.008 #choose initial time for growth rate calculation
    t2 = 0.012 #choose final time for growth rate calculation
    analytic_growth_rate=  (abs(params.m2^2 - params.m1^2)/ (2*CCNO.hbar* Eνₑ)) # analytic growth rate 
    println("analytic_growth_rate:",analytic_growth_rate)

    x = CCNO.generate_x_array(N_sites_eachflavor, L)
    y = CCNO.generate_x_array(N_sites_eachflavor, L)
    z = CCNO.generate_x_array(N_sites_eachflavor, L)

    # p matrix with numbers generated from the p_array for all components (x, y, z) #sherood has 
    p = hcat(CCNO.generate_px_array(params.N_sites, Eνₑ, Eνₑ̄), CCNO.generate_py_array(params.N_sites), CCNO.generate_pz_array(params.N_sites))

    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= params.N_sites ÷ 2 ? -1 : 1 for i in 1:params.N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

    s = siteinds("S=1/2", params.N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

    # Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up in Richers notation which is equivalently half spin up and half chain spin down in my TN notation)
    ψ = productMPS(s, n -> n <= params.N_sites/2 ? "Up" : "Dn")

    N = CCNO.Neutrino_number(params, n_νₑ,n_νₑ̄)

    state = CCNO.SimulationState(ψ=ψ,
                                 s=s,
                                 p=p,
                                 energy_sign = energy_sign,
                                 N=N,
                                 xyz = hcat(x,y,z))

    ρₑμ_at_t1 = nothing  # Initialize a variable to store ρₑμ at t1
    ρₑμ_at_t2 = nothing  # Initialize a variable to store ρₑμ at t2
    Δt = t2 - t1 #time difference between growth rates

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
        
    # Loop over the time array to match t1 and t2
    for (i, t) in enumerate(t_array) 
        # Check if the current time is approximately t1
        if abs(t - t1) < params.τ / 2
            println("corresponding ρₑμ index from the time array =",i)
            ρₑμ_at_t1 = ρₑμ_array[i]
            println("ρₑμ_at_t1=",ρₑμ_at_t1)
        end

        # Check if the current time is approximately t2
        if abs(t - t2) < params.τ / 2
            println("corresponding ρₑμ index from the time array =",i)
            ρₑμ_at_t2 = ρₑμ_array[i]
            println("ρₑμ_at_t2=",ρₑμ_at_t2)
        end
    end

    # After the time evolution loop, calculate and print the growth rate
    if ρₑμ_at_t1 !== nothing && ρₑμ_at_t2 !== nothing
        Im_Ω = (1 / Δt) * log(ρₑμ_at_t2 / ρₑμ_at_t1)
        println("Growth rate of flavor coherence of ρₑμ at t2 to ρₑμ at t1: $Im_Ω")
    else
        println("ρₑμ was not captured at both t1 and t2.")
    end

    if params.save_plots_flag 
        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

        CCNO.save_plots(params, s,L,t_array, Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array)
    end
    if !params.save_plots_flag 
        # Plotting ρₑμ vs t # for jenkins file 
        plot(t_array, ρₑμ_array, xlabel = "t", ylabel = "<ρₑμ>_1", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig( "Homo_MB_<ρₑμ>_site1_vs_t for $(params.N_sites) particles.pdf")
    end

    # comment out assert because we don't actually have analytic solution for MB case.
    #@assert abs((Im_Ω - analytic_growth_rate)/  analytic_growth_rate) < params.tolerance 
end

@time main()
