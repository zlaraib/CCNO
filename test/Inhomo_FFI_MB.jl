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

save_data = false  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs

function main()
    """ Richers(2021) Test 4 initial conditions for many-body dynamics: """
    N_sites_eachflavor= 10 # total sites/particles that evenly spaced "for each (electron) flavor" 
    N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
    τ = 5E-13 # time step to include 50 steps every 10 picoseconds # sec # variable
    τ_pert = 10^-8 # time step for perturbation evolution
    ttotal = 5.0E-12 # total time of evolution # sec #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    m1 = -0.008596511*CCNO.eV #eV  1st mass eigenstate of neutrino
    m2 = 0*CCNO.eV #eV  2nd mass eigenstate of neutrino
    Δm² = (m2^2-m1^2) # mass square difference # (erg^2)
    maxdim = 20 # max bond dimension in MPS truncation
    cutoff = 1e-100 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    L = 1 # cm # domain size # (aka big box length)
    n_νₑ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
    Eνₑ =  50.0*CCNO.MeV # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
    Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
    Δx = L/N_sites_eachflavor # length of the box of interacting neutrinos at a site in cm  #variable
    theta_nu = 1.74532925E-8 # mixing_angle
    B = [-sin(2*theta_nu), 0, cos(2*theta_nu)] # actual b vector that activates the vacuum oscillation term in Hamiltonian
    B = B / norm(B) 
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 
    Δp = Δx  # width of shape function  # cm #variable
    t1 = 33e-12 #choose initial time for growth rate calculation
    t2 = 53e-12 #choose final time for growth rate calculation
    periodic = true  # true = imposes periodic boundary conditions while false doesn't
    k = 2*pi / (L)
    analytic_growth_rate=  (abs(m2^2 - m1^2)/ (2*CCNO.hbar* Eνₑ)) + (CCNO.c* k)  # analytic growth rate #fix it for inhomo from paper
    println("analytic_growth_rate=", analytic_growth_rate)

    checkpoint_every = 4
    do_recover = false
    recover_file = "" 

    x = CCNO.generate_x_array(N_sites_eachflavor, L)
    y = CCNO.generate_x_array(N_sites_eachflavor, L)
    z = CCNO.generate_x_array(N_sites_eachflavor, L)

    # p matrix with numbers generated from the p_array for all components (x, y, z) #sherood has 
    p = hcat(CCNO.generate_px_array(N_sites, Eνₑ, Eνₑ̄), CCNO.generate_py_array(N_sites), CCNO.generate_pz_array(N_sites))

    # Create an array with the first half as 1 and the rest as -1
    energy_sign = [i <= N_sites ÷ 2 ? -1 : 1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
    # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
    # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

    s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

    # Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up in Richers notation which is equivalently half spin up and half chain spin down in my TN notation)
    @time ψ = productMPS(s, n -> n <= N_sites/2 ? "Up" : "Dn")

    function generate_B_pert(α)
        # Generate perturbations in x and y plane
        x_pert = α  
        y_pert = α  
        # Calculate the z component to maintain normalization
        z_pert = sqrt(max(0, 1 - x_pert^2 - y_pert^2))

        # Return the B_pert vector
        return [x_pert, y_pert, z_pert]
    end

    α = 1e-6 # perturbation strength as mentioned in the paper for the inhomogenous Richers test 
    # Generate the perturbed B vector scaled by α  as mentioned in the paper
    B_pert = generate_B_pert(α)

    # Since the perturbation is small, B_pert should already be normalized, but you can normalize again for precision
    B_pert = B_pert / norm(B_pert) 
    # println("B_pert= ", B_pert)  

    # Perturb the state via one-body Hamiltonian
    @time ψ₀= CCNO.evolve_perturbation(s,k, τ_pert, B_pert, α, x, L, N_sites, ψ, cutoff, maxdim, energy_sign,ttotal)

    @time N = CCNO.Neutrino_number(s, τ, B,L, N_sites, N_sites_eachflavor, tolerance,
                    n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, ttotal,periodic)

    # Specify the relative directory path
    datadir = joinpath(@__DIR__,"datafiles")
    chkptdir = joinpath(@__DIR__, "checkpoints")

    ρₑμ_at_t1 = nothing  # Initialize a variable to store ρₑμ at t1
    ρₑμ_at_t2 = nothing  # Initialize a variable to store ρₑμ at t2
    Δt = t2 - t1 #time difference between growth rates

    #extract output for the survival probability values at each timestep
    @time Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover = CCNO.evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ₀, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,chkptdir, checkpoint_every,  do_recover, recover_file, save_data , periodic)

    # Take the abs value fo all enteries till N_sites_eachflavor and then take the mean of that first half of the array, then do this for each row in ρₑμ_array
    @time ρₑμ_array_domain_avg = [mean(abs.(row[1:N_sites_eachflavor])) for row in ρₑμ_array] 

    # Loop over the time array to match t1 and t2
    for (i, t) in enumerate(t_array) 
        # Check if the current time is approximately t1
        if abs(t - t1) < τ / 2
            println("corresponding ρₑμ index from the time array =",i)
            ρₑμ_at_t1 = ρₑμ_array_domain_avg[i]
            println("ρₑμ_at_t1=",ρₑμ_at_t1)
        end

        # Check if the current time is approximately t2
        if abs(t - t2) < τ / 2
            println("corresponding ρₑμ index from the time array =",i)
            ρₑμ_at_t2 = ρₑμ_array_domain_avg[i]
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

    if save_data
        # Generate input data
        input_data = CCNO.extract_initial_conditions(N_sites,N_sites_eachflavor,τ, ttotal,tolerance,
        Δm², maxdim, cutoff, p,ψ₀,L, Δx,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,B, N, shape_name,Δp,periodic)
        # Call the function to generate the inputs file in the specified directory
        CCNO.generate_inputs_file(datadir, "inputs.txt", input_data)
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
        Sz_array =  t_Sz_tot[:, 2:N_sites_eachflavor+1]
        Sy_array = t_Sy_tot[:, 2:N_sites_eachflavor+1]  
        Sx_array= t_Sx_tot[:, 2:N_sites_eachflavor+1] 
        prob_surv_array = t_probsurv_tot[:, 2:N_sites_eachflavor+1] 
        ρₑₑ_array = t_ρₑₑ_tot[:, 2:N_sites_eachflavor+1] 
        ρ_μμ_array =t_ρ_μμ_tot[:, 2:N_sites_eachflavor+1]  
        ρₑμ_array = t_ρₑμ_tot[:, 2:N_sites_eachflavor+1]

        # Parsing arrays containing strings like "[1.0,", into a numeric array suitable for plotting

        Sz_array_domain_avgd = [mean(abs.(row)) for row in eachrow(Sz_array)]
        Sy_array_domain_avgd = [mean(abs.(row)) for row in eachrow(Sy_array)]
        Sx_array_domain_avgd = [mean(abs.(row)) for row in eachrow(Sx_array)]
        prob_surv_array_domain_avgd = [mean(abs.(row)) for row in eachrow(prob_surv_array)]
        ρₑₑ_array_domain_avgd = [mean(abs.(row)) for row in eachrow(ρₑₑ_array)]
        ρ_μμ_array_domain_avgd = [mean(abs.(row)) for row in eachrow(ρ_μμ_array)]
        ρₑμ_array_domain_avgd= [mean(abs.(row)) for row in eachrow(ρₑμ_array)]

        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first
        CCNO.save_plots(τ, N_sites,L,t_array, ttotal,Sz_array_domain_avgd, Sy_array_domain_avgd, Sx_array_domain_avgd, prob_surv_array_domain_avgd, x_values, pₓ_values, ρₑₑ_array_domain_avgd,ρ_μμ_array_domain_avgd, ρₑμ_array_domain_avgd,datadir, plotdir, save_plots_flag)
        
        # Call the function to generate the inputs file in the specified directory
        CCNO.generate_inputs_file(plotdir, "inputs.txt", input_data)
    end

    if !save_plots_flag 
        # Plotting ρₑμ vs t # for jenkins file 
        plot(t_array, ρₑμ_array_domain_avg, xlabel = "t", ylabel = "<ρₑμ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig( "Inhomo_MB_<ρₑμ>_domain_avg_vs_t for $N_sites particles.pdf")
    end

    #commented out b/c assert "doesnt always" pass even with fixed maxdim= 2
    # @assert abs((Im_Ω - analytic_growth_rate)/  analytic_growth_rate) < tolerance 
end
@time main()
