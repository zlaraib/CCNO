push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5

# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()

    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm
    delta_omega = 0.5
    Δm²= 1.0 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.

    params = CCNO.Parameters(
        N_sites = 4, # number of sites 
        cutoff = 1E-14, # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
        τ = 0.05*CCNO.hbar, # time step 
        ttotal = 5*CCNO.hbar, # total time of evolution 
        tolerance  = 5E-1, # acceptable level of error or deviation from the exact value or solution
        Δx = Δx, # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm
        maxdim = 1000, #bond dimension
        m1 = 0.0,
        m2 = sqrt(Δm²),
        L=L,
        Δp = L, # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.  
        periodic = true,  # true = imposes periodic boundary conditions while false doesn't
        checkpoint_every = 20,
        do_recover = false,
        recover_file = "",
        shape_name = "none",  # Change this to the desired shape name # variable.
        geometric_name = "none",
        datadir = joinpath(@__DIR__, "datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        plotdir = joinpath(@__DIR__, "plots"),
        save_plots_flag = false,
        α = 0,
        theta_nu = 0 # mixing_angle #rad 
    )
 
    # s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all params.N_sites.
    # conserve_qns=false doesnt conserve the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  

    # Constants for Rogerro's fit (corresponding to Δω = 0.25)
    a_t = 1.224
    b_t = 0
    c_t = 1.62
    
    # Initialize an array of ones for all N_sites sites
    mu = ones(params.N_sites) # erg
    
    # Create an array of dimension N_sites and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    N = mu .* fill(((Δx)^3 )/(√2 * CCNO.G_F * params.N_sites), params.N_sites)
    
    B = [sin(2*params.theta_nu), 0, -cos(2*params.theta_nu)] # is equivalent to B = [0, 0, -1] # fixed for Rogerro's case
    B = B / norm(B)

    x = fill(rand(), params.N_sites) # variable.
    y = fill(rand(), params.N_sites) # variable.
    z = fill(rand(), params.N_sites) # variable.

    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Dn" : "Up")

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    # set energy such that omega_a = delta_omega and omega_b=0
    px_a = fill(Δm²/(2.0*delta_omega), div(params.N_sites,2))
    px_b = fill(10^6, div(params.N_sites,2))
    px = vcat(px_a, px_b)
    p = hcat(px,fill(0, params.N_sites), fill(0, params.N_sites))
    energy_sign = [i <= params.N_sites ÷ 2 ? 1 : 1 for i in 1:params.N_sites] # all of the sites are neutrinos

    state = CCNO.SimulationState(ψ=ψ,
                                 s=s,
                                 p=p,
                                 energy_sign = energy_sign,
                                 N=N,
                                 xyz = hcat(x,y,z))

    #extract output for the survival probability values at each timestep
    CCNO.evolve(params, state)

    # This function scans through the array, compares each element with its neighbors, 
    # and returns the index of the first local minimum it encounters. 
    # If no local minimum is found, it returns -1 to indicate that.
    function find_first_local_minima_index(arr)
        N = length(arr)
        for i in 2:(N-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
    
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
    
    # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
    tmin_ifirstlocalmin_file = joinpath(params.datadir, "tmin_ifirstlocalmin.dat")
    # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
    i_first_local_min = find_first_local_minima_index(-Sz_array)
    
    # Writing if_else statement to communicate if local minima (not) found
    t_min = nothing
    if i_first_local_min != -1
        println("Index of the first local minimum: ", i_first_local_min)
        t_min = params.τ * i_first_local_min - params.τ
    else
        println("No local minimum found in the array.")
    end
    println("Corresponding time of first minimum index= ", t_min/CCNO.hbar)

    # Save t_min and i_first_local_min to file
    writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
    println("Saved t_min to file: ", tmin_ifirstlocalmin_file)

    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(params.N_sites) + b_t * sqrt(params.N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

    if params.save_plots_flag
        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first
        CCNO.save_plots(params.τ, params.N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,params.datadir, plotdir, params.save_plots_flag)
    end
    if !params.save_plots_flag
        # Plotting P_surv vs t
        plot(t_array/CCNO.hbar, -Sz_array, xlabel = "t", ylabel = "S_z",title = "Running main_Rogerro script \n for N_sites$(params.N_sites) with maxdim=1 MF for τ$(params.τ)", legend = false, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label= ["My_plot_for_N$(params.N_sites)"]) 
        scatter!([t_p_Rog],[-Sz_array[i_first_local_min]], label= ["t_p_Rog"])
        scatter!([t_min/CCNO.hbar],[-Sz_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:bottomleft)
        # Save the plot as a PDF file
        savefig("main_Roggero.pdf")
    end

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
@assert abs(t_min/CCNO.hbar - t_p_Rog) < params.τ/CCNO.hbar + params.tolerance

end 

@time main()

