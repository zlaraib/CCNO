push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics
using HDF5

""" Richers(2021) Test 2 initial conditions: """
function main()
    N_sites_eachflavor= 1 # total sites/particles that evenly spaced "for each (electron) flavor" 
    t_bipolar = 8.96e-4 #characteristic bipolar time #sec
    L = 1e7 # cm # domain size # (aka big box length)
    Δx = L # length of the box of interacting neutrinos at a site in cm 

    params = CCNO.Parameters(
        N_sites = 2* (N_sites_eachflavor),
        τ = 1e-7, # time step # sec/sec = unitless # variable # using this time step for faster unit testing in jenkins, actually the bipolar richers results are obtained with timestep= 1e-9/t_bipolar.
        ttotal = 0.002, # total time of evolution # sec/sec = unitless #using this total time for faster unit testing in jenkins, actually bipolar richers results are produced with ttotal = 0.01 / t_bipolar
        tolerance  = 5E-1, # acceptable level of error or deviation from the exact value or solution #variable
        m2 = 0.008596511*CCNO.eV, #ergs #1st mass eigenstate of neutrino in Richers(2021)
        m1 = 0*CCNO.eV,   #ergs #2nd mass eigenstate of neutrino in Richers(2021)
        maxdim = 1, # max bond dimension in MPS truncation
        cutoff = 1e-100, # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
        shape_name = "none",  # Change this to the desired shape name #variable
        geometric_name = "physical",
        Δp = L, # width of shape function  # cm #variable
        Δx=Δx,
        L=L,
        periodic = true,  # true = imposes periodic boundary conditions while false doesn't
        theta_nu= 0.01, #mixing angle # =34.3 degrees
        checkpoint_every = 100000,
        do_recover = false,
        recover_file = "",
        datadir = joinpath(@__DIR__,"datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        save_plots_flag = false,
        α = 0
    )
    
    Δm² = (params.m2^2-params.m1^2) # mass square difference # (erg^2) #not used in the test
    Eνₑ =  50.0*CCNO.MeV #ergs # energy of all neutrinos (P.S the its negative is energy of all antineutrinos) 
    Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
    n_νₑ =  (10* (params.m2-params.m1)^2 )/(2*(√2)*CCNO.G_F*Eνₑ) # cm^-3 # number density of electron flavor neutrino
    n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    t1 = 0.0084003052 #choose initial time for growth rate calculation
    t2 = 0.011700318 #choose final time for growth rate calculation
    B = [-sin(2 *params.theta_nu), 0, cos(2*params.theta_nu)] # for inverted mass hierarchy
    B = B / norm(B) 


    x = CCNO.generate_x_array(N_sites_eachflavor, L)
    y = fill(0, params.N_sites) #variable
    z = fill(0, params.N_sites) #variable

    p = hcat(CCNO.generate_px_array(params.N_sites, Eνₑ, Eνₑ̄), CCNO.generate_py_array(params.N_sites), CCNO.generate_pz_array(params.N_sites))

    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= params.N_sites ÷ 2 ? -1 : 1 for i in 1:params.N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

    s = siteinds("S=1/2", params.N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

    # Initialize psi to be a product state (Of first half electron flavor neutrino i.e. spin up while other half anti-neutrinos electron flavor i.e. half spin down)
    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Up" : "Dn")

    # @time main(s, τ, B,L, N_sites, N_sites_eachflavor, tolerance,
    # n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, ttotal,periodic)
    N = CCNO.Neutrino_number(params, n_νₑ,n_νₑ̄)

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
        
        CCNO.save_plots(params, L,t_array, Sz_array, Sy_array, Sx_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array)
        # Call the function to generate the inputs file in the specified directory
        CCNO.generate_inputs_file(params, "inputs.txt")
    end
    if !params.save_plots_flag
        # Plotting ρ_ee vs t # for jenkins file 
        plot(t_array/t_bipolar, ρₑₑ_array, xlabel = "t", ylabel = "<ρₑₑ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig("Bipolar Richers_site1 for $(params.N_sites) particles <ρₑₑ>_vs_t.pdf")
    end

    # clean up
    rm(params.datadir, recursive=true)
    rm(params.chkptdir, recursive=true)
end 

@time main()
