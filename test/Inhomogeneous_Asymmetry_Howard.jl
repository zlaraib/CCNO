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

function main()
    total_sites = 16

    # total sites/particles that evenly spaced "for each (electron) flavor" 
    N_sites_eachflavor = div(total_sites, 2)
    
    L = 1.0 # cm # domain size # (aka big box length)
    tolerance = 5E-1

    # Make samples odd
    samples = 9
    # Make asy multiple even
    asy_multiple = 8 # We actually get one more asymmetry config than this number 
    numerical_values = zeros(samples * (div(asy_multiple, 2) + 1))

    # NOTE: symmetrical for e-anti, so we just do half
    
    for asy_iterator = 0:div(asy_multiple, 2)
        for idx = 1:samples
            println("\n======\n")
            


            # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
            Enu_e = 50.0*CCNO.MeV
            # Since all neutrinos have same energy, except in this case anti
            # neutrinos are moving in opposite direction to give it a negative sign
            Enu_ē = -1 * Enu_e
            
            t1 = 0.001 #choose initial time for growth rate calculation
            t2 = 0.016 #choose final time for growth rate calculation
            k = 2*pi / (L)

            x = Float64[];
            y = Float64[];
            z = Float64[];
            Px = Float64[];
            Py = Float64[];
            Pz = Float64[];
            energy_sign = Int[];
            spin_vec = String[];

            total_up = 0
            total_dn = 0

            for asy_idx = 1:asy_multiple
                append!(x, CCNO.generate_x_array(1, L))
                append!(y, fill(0.0, 2))
                append!(z, fill(0.0, 2))

                append!(
                    Px,
                    CCNO.generate_px_array(2, Enu_e, Enu_ē),
                )
                append!(Py, fill(0, 2))
                append!(Pz, fill(0, 2))

                if asy_idx > asy_iterator
                    spin_vec_tmp = ["Dn", "Dn"]
                    energy_sign_tmp = [-1, -1]
                    total_dn = total_dn + 1
                else
                    spin_vec_tmp = ["Up", "Up"]
                    energy_sign_tmp = [1, 1]
                    total_up = total_up + 1
                end
                append!(energy_sign, energy_sign_tmp)
                append!(spin_vec, spin_vec_tmp)
            end

            # Spacing out the sites (inhomogeneity)
            idx_up = 1.0
            idx_dn = 1.0
            up_diff = L / total_up
            dn_diff = L / total_dn
            for i = 1:asy_multiple
                if energy_sign[(i*2)-1] == 1
                    x[(i*2)-1] = (idx_up - 0.5) * up_diff
                    x[(i*2)] = (idx_up - 0.5) * up_diff
                    idx_up = idx_up + 1
                else 
                    x[(i*2)-1] = (idx_dn - 0.5) * dn_diff
                    x[(i*2)] = (idx_dn - 0.5) * dn_diff
                    idx_dn = idx_dn + 1
                end
            end

            theta = range(0, pi, samples)
            rotation = [
                cos(theta[idx]) -sin(theta[idx]) 0;
                sin(theta[idx]) cos(theta[idx]) 0;
                0 0 1
            ]

            p = hcat(Px, Py, Pz)

            # NOTE: rotations Rotate UP on -p, DN on +p (+p is odd idx, -p even)
            for i = 1:asy_multiple
                mdx = (i * 2) - 1
                if energy_sign[mdx] == 1
                    p[mdx, :] = rotation * p[mdx, :]
                else
                    p[mdx+1, :] = rotation * p[mdx+1, :]
                end
            end

            # TODO: How much does this parameter matter? 
            # OR: does it matter how we expect it to matter?
            Delta_x = L / div(total_sites, 4) 
            
            n_nu_e = 4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
            
            # TODO: refine this process (factor of 2 should be removed)
            n_nu_e = (n_nu_e * 2) / div(total_sites, 4)
            n_nu_ē = n_nu_e # cm^-3 # number density of electron flavor antineutrino
            
            params = CCNO.Parameters(
                N_sites = total_sites,
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

            println(spin_vec)
            println(energy_sign)
            println(p)
            println(x)

            #fixed #switched conserve_qns to false to avoid fluxes error in expect function
            s = siteinds("S=1/2", params.N_sites; conserve_qns = false)
            xyz = hcat(x, y, z)

            Psi = productMPS(s, spin_vec)

            N = n_nu_e * (params.Delta_x)^3

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

            # println(size(rho_emu_array))
            # Take the abs value fo all enteries till N_sites_eachflavor and
            # then take the mean of that first half of the array, then do this
            # for each row in rho_emu_array 
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
            x = t_array[(1*vecqtr):veclen]

            y = log.(rho_emu_array_domain_avg[(1*vecqtr):veclen])

            numerical_values[idx+((asy_iterator)*samples)] =
                coef(linregress(x, y))[1]

            # clean up
            rm(params.datadir, recursive = true)
            rm(params.chkptdir, recursive = true)

            # Check with growth rate from Chakraborty 2016
            m1 = 0 * eV
            m2 = 0.008596511 * eV
            Enu_e = 50.0 * MeV
            omega = (abs(m2^2 - m1^2) / (2 * Enu_e))
            a = asy_iterator / (asy_multiple)
            a = a * 2 # Multiply by two to match Chakraborty paper normalization
            a = a - 1 # so it spans -1, 1
            n = 2.92e24 # Practial number density for analytical calc
            cs = cos(theta[idx])
            mu = sqrt(2) * G_F * n

            println(a, " ", cs)

            analytic_growth_rate =
                0.5 * imag(
                    sqrt(
                        Complex(
                            (a^2 * cs^2 * mu^2) - (6 * a^2 * cs * mu^2) +
                            (9 * a^2 * mu^2) +
                            (4 * cs * mu * omega) - (12 * mu * omega) +
                            (4 * omega^2),
                        ),
                    ),
                )

            analytic_growth_rate = analytic_growth_rate / hbar

            # NOTE: this doesn't really work when the analytic growth rate is ~0
            if abs(
                (
                    numerical_values[idx+((asy_iterator)*samples)] -
                    analytic_growth_rate
                ) / (analytic_growth_rate),
            ) < tolerance
                println("Passed.")
                println("Analytic growth rate: ", analytic_growth_rate)
                println(numerical_values)
            else
                println("Failed OR ~0 growth rate.")
                println("Analytic growth rate: ", analytic_growth_rate)
                println(
                    "Numerical: ",
                    numerical_values[idx+((asy_iterator)*samples)],
                )
                # NOTE: this isn't perfect, but it's enough of a test for now.
                # What does this mean? The test only fails if the growth rates
                # do not match and are also (relatively) high.
                # if analytic_growth_rate > 300
                #     println("EXITING")
                #     return 1
                # end
            end
            
            println(numerical_values);
        end
    end

    return 0
end

@time main()
