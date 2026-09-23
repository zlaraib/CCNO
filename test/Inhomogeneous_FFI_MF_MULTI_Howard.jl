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

using Plots
using LinearRegression

const hbar::Float64 = 1.05457266e-27 # erg s
const c::Float64 = 2.99792458e10 # cm/s
const eV::Float64 = 1.60218e-12 # erg
const MeV::Float64 = 1e6 * eV # erg
const GeV::Float64 = 1e9 * eV # erg
const G_F::Float64 = 1.1663787e-5 / GeV^2 * (hbar*c)^3 # erg cm^3
const kB::Float64 = 1.3807e-16 # erg/K

function analytic(samples, numerical_values)
    x_data = range(0, pi, samples)
    n_nu_e = 4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    mu = sqrt(2) * G_F * n_nu_e
    L = 1.0
    k = 2 * pi / L
    k = k * hbar * c
    
    cs = cos.(x_data) # effectively cos theta
    
    y_data = imag.(
                   sqrt.(Complex.(
                                  k .* ((cs .* mu) .- (3.0 .* mu) .+ k)
                                 )),
    )
    
    y_data = y_data ./ hbar

    plot = plot!(
        cos.(x_data),
        # [y_data, numerical_values],
        # label = ["Analytic" "Numerical"],
        # NOTE: analytical result removed as it may not be correct
        [numerical_values],
        label = ["Numerical"],
        xlabel = "cos(theta)",
        ylabel = "Growth Rate",
        title = "Analytic growth rate vs incident angle.",
    )

    return savefig("Inhomogeneous_MF_incident_vs_growth.pdf")
end

function main()
    zone_number = 16 # quick fix to set the necessary zone number for now...
    half_zone_number = 8
    qtr_zone_number = 4

    N_sites_eachflavor = half_zone_number

    L = 1.0 # cm # domain size # (aka big box length)
    Delta_x = L/N_sites_eachflavor # length of the box of interacting neutrinos at a site in cm  #variable
    tolerance = 5E-1

    Delta_x = Delta_x * 2 # Sites "take up twice as much space" with 4 beams (instead of 2)

    # NOTE: odd number of samples to include cos(theta) = 0
    # We use a particularly low number of samples in this case to reduce computational load.
    samples = 3
    numerical_values = zeros(samples)

    for idx = 1:samples
        params = CCNO.Parameters(
            N_sites = zone_number,
            τ = 5E-13,
            ttotal = 9.0E-11, 
            m1 = 0.0,
            m2 = 0.0,
            maxdim = 1,
            cutoff = 1e-100,
            theta_nu = 1.74532925E-8,
            shape_name = "triangular",
            homogeneous = [false, true, true], # homogeneous in y, z dims
            geometric_name = "physical",
            Delta_x = Delta_x,
            L = L,
            Delta_p = Delta_x,
            periodic = true,
            checkpoint_every = 20,
            do_recover = false,
            recover_file = "",
            plotdir = joinpath(@__DIR__, "plots"),
            datadir = joinpath(@__DIR__, "datafiles-" * string(idx)),
            chkptdir = joinpath(@__DIR__, "checkpoints-" * string(idx)),
            save_plots_flag = false,
            alpha = 1e-6,
        )

        Delta_m_squared = (params.m2^2-params.m1^2) # mass square difference # (erg^2)
        n_nu_e = 4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
        n_nu_ē = n_nu_e # cm^-3 # number density of electron flavor antineutrino
        Enu_e = 50.0*CCNO.MeV # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
        Enu_ē = -1 * Enu_e # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
        #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
        t1 = 33e-12 #choose initial time for growth rate calculation
        t2 = 53e-12 #choose final time for growth rate calculation
        k = 2*pi / (L)
        analytic_growth_rate =
            (abs(params.m2^2 - params.m1^2) / (2 * CCNO.hbar * Enu_e)) +
            (CCNO.c * k)  # analytic growth rate #fix it for inhomo from paper
        println("analytic_growth_rate=", analytic_growth_rate)
        
        x = [
            CCNO.generate_x_array(qtr_zone_number, L);
            CCNO.generate_x_array(qtr_zone_number, L)
        ]
        y = [
            fill(0.0, half_zone_number);
            fill(0.0, half_zone_number)
        ]
        z = [fill(0.0, half_zone_number); fill(0.0, half_zone_number)]

        xyz = hcat(x, y, z) 

        theta = range(0, pi, samples)
        
        p = hcat(
            [
                CCNO.generate_px_array(half_zone_number, Enu_e, Enu_ē);
                CCNO.generate_px_array(half_zone_number, Enu_e, Enu_ē)
            ],
            [
                fill(0, half_zone_number);
                fill(0, half_zone_number);
            ],
            [fill(0, half_zone_number); fill(0, half_zone_number)],
        )

        
            rotation = [cos(theta[idx]) -sin(theta[idx]) 0; sin(theta[idx]) cos(theta[idx]) 0; 0 0 1]
        for j in half_zone_number+1:zone_number
            p[j, :] = rotation * p[j, :]
        end

        
        energy_sign = [-1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1]
        
        # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
        # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
        # conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
        # i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
        # conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

        s = siteinds("S=1/2", params.N_sites; conserve_qns = false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

        # Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up in Richers notation which is equivalently half spin up and half chain spin down in my TN notation)
        # Psi = productMPS(s, n -> n <= params.N_sites/2 ? "Up" : "Dn")

        spin_vec = ["Up", "Up", "Up", "Up", "Dn", "Dn", "Dn", "Dn", "Dn", "Dn", "Dn", "Dn", "Up", "Up", "Up", "Up"]
        Psi = productMPS(s, spin_vec)

        N = CCNO.Neutrino_number(params, n_nu_e, n_nu_ē)

        state = CCNO.SimulationState(
            Psi = Psi,
            s = s,
            s0 = s,
            p = p,
            energy_sign = energy_sign,
            N = N,
            xyz = xyz,
        )

        # Perturb the state via one-body Hamiltonian
        CCNO.perturb(params, state, k, params.theta_nu)

        #extract output for the survival probability values at each timestep
        CCNO.evolve(params, state)

        # Read the data files #

        t_Sz_tot = readdlm(joinpath(params.datadir, "t_<Sz>.dat"))
        t_Sy_tot = readdlm(joinpath(params.datadir, "t_<Sy>.dat"))
        t_Sx_tot = readdlm(joinpath(params.datadir, "t_<Sx>.dat"))
        t_xsiteval = readdlm(joinpath(params.datadir, "t_xsiteval.dat"))
        t_pxsiteval = readdlm(joinpath(params.datadir, "t_pxsiteval.dat"))
        t_rho_e_e_tot = readdlm(joinpath(params.datadir, "t_rho_e_e.dat"))
        t_rho_mumu_tot = readdlm(joinpath(params.datadir, "t_rho_mumu.dat"))
        t_rho_emu_tot = readdlm(joinpath(params.datadir, "t_rho_emu.dat"))

        # Extract time array and corresponding values for plotting
        t_array = t_Sz_tot[:, 1]
        Sz_array = t_Sz_tot[:, 2:(N_sites_eachflavor+1)]
        Sy_array = t_Sy_tot[:, 2:(N_sites_eachflavor+1)]
        Sx_array = t_Sx_tot[:, 2:(N_sites_eachflavor+1)]
        rho_e_e_array = t_rho_e_e_tot[:, 2:(N_sites_eachflavor+1)]
        rho_mumu_array = t_rho_mumu_tot[:, 2:(N_sites_eachflavor+1)]
        rho_emu_array = t_rho_emu_tot[:, 2:(N_sites_eachflavor+1)]

        println(size(rho_emu_array))
        # Take the abs value fo all enteries till N_sites_eachflavor and then take the mean of that first half of the array, then do this for each row in rho_emu_array 
        rho_emu_array_domain_avg = mean(abs.(rho_emu_array), dims = 2)

        rho_emu_at_t1 = nothing  # Initialize a variable to store rho_emu at t1
        rho_emu_at_t2 = nothing  # Initialize a variable to store rho_emu at t2
        Delta_t = t2 - t1 #time difference between growth rates

        # NOTE: the idea is to fit a line to the last quarter of the evolution. 
        # There is surely a better way to do this, but this works for the time being.
        
        veclen = size(rho_emu_array_domain_avg)[1]
        vecqtr = (veclen - (veclen % 4)) / 4 # mod for divisi
        vecqtr = Int(vecqtr)

        t_array = t_Sz_tot[:, 1]
        x = t_array[(2*vecqtr):veclen]
        
        y = log.(rho_emu_array_domain_avg[(2*vecqtr):veclen])

        numerical_values[idx] = coef(linregress(x, y))[1]
        
        # clean up
        rm(params.datadir, recursive = true)
        rm(params.chkptdir, recursive = true)
    
        # Assert checks against a full rotation, where we expect the result to match the maxiumum growth rate case. (Directly opposed beams, maximum asymmetry)
        if idx == samples
            @assert abs((numerical_values[idx] - analytic_growth_rate)/  analytic_growth_rate) < tolerance 
        end
        
    end

    analytic(samples, numerical_values)

    return 0
end

@time main()
