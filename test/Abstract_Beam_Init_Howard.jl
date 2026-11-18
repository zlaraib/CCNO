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
const Enu::Float64 = 50.0*CCNO.MeV
const L::Float64 = 1.0
const n_nu_e::Float64 = 4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino

function run(beam_arr, num_beams, bonddim)
    
    params, state = CCNO.generate_params_from_beamarr(beam_arr, L, n_nu_e, bonddim)
    
    k = 2*pi / (L)
    analytic_growth_rate = (abs(params.m2^2 - params.m1^2)/ (2*CCNO.hbar* Enu)) + (CCNO.c* k)  # analytic growth rate #fix it for inhomo from paper
    println("analytic_growth_rate=", analytic_growth_rate)

    Delta_m_squared = (params.m2^2-params.m1^2) # mass square difference # (erg^2)


    # Perturb the state via one-body Hamiltonian
    CCNO.perturb(params, state, k, params.theta_nu)

    #extract output for the survival probability values at each timestep
    CCNO.evolve(params, state)

    # Site info for determining growth rate for the Jenkins test ONLY
    # This will not be valid for all configurations.
    N_sites_eachflavor = beam_arr[1].num_sites
    
    #=====================#
    # Read the data files #
    #=====================#
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
    Sz_array =  t_Sz_tot[:, 2:N_sites_eachflavor+1]
    Sy_array = t_Sy_tot[:, 2:N_sites_eachflavor+1]  
    Sx_array= t_Sx_tot[:, 2:N_sites_eachflavor+1] 
    rho_e_e_array = t_rho_e_e_tot[:, 2:N_sites_eachflavor+1] 
    rho_mumu_array =t_rho_mumu_tot[:, 2:N_sites_eachflavor+1]  
    rho_emu_array = t_rho_emu_tot[:, 2:N_sites_eachflavor+1]

    # Take the abs value fo all entries till N_sites_eachflavor and then take the mean of that first half of the array, then do this for each row in rho_emu_array 
    rho_emu_array_domain_avg = mean(abs.(rho_emu_array), dims=2) 
    
    t1 = 33e-12 #choose initial time for growth rate calculation
    t2 = 53e-12 #choose final time for growth rate calculation

    rho_emu_at_t1 = nothing  # Initialize a variable to store rho_emu at t1
    rho_emu_at_t2 = nothing  # Initialize a variable to store rho_emu at t2
    Delta_t = t2 - t1 #time difference between growth rates

    # Loop over the time array to match t1 and t2
    for (i, t) in enumerate(t_array) 
        # Check if the current time is approximately t1
        if abs(t - t1) < params.τ / 2
            println("corresponding rho_emu index from the time array =",i)
            rho_emu_at_t1 = rho_emu_array_domain_avg[i]
            println("rho_emu_at_t1=",rho_emu_at_t1)
        end

        # Check if the current time is approximately t2
        if abs(t - t2) < params.τ / 2
            println("corresponding rho_emu index from the time array =",i)
            rho_emu_at_t2 = rho_emu_array_domain_avg[i]
            println("rho_emu_at_t2=",rho_emu_at_t2)
        end
    end

    # After the time evolution loop, calculate and print the growth rate
    if rho_emu_at_t1 !== nothing && rho_emu_at_t2 !== nothing
        Im_Ω = (1 / Delta_t) * log(rho_emu_at_t2 / rho_emu_at_t1)
        println("Growth rate of flavor coherence of rho_emu at t2 to rho_emu at t1: $Im_Ω")
    else
        println("rho_emu was not captured at both t1 and t2.")
    end
    
    # clean up
    rm(params.datadir, recursive = true)
    rm(params.chkptdir, recursive = true)

    return 0
end

function main()

    num_beams = 2
    beam_arr = CCNO.init_beam_arr(num_beams)

    beam_arr[1] = CCNO.Beam(num_sites = 10, rotation = 0.0, energy = Enu, flavor = 1)
    beam_arr[2] = CCNO.Beam(num_sites = 10, rotation = pi, energy = Enu, flavor = -1)

    bonddim = 1
    run(beam_arr, num_beams, bonddim)
    
    return 0
end


@time main()
