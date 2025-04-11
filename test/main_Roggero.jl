push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

# array of initial parameters
const delta_omega_list::Vector{Float64} = [-0.5, 0.0, 1.0, 1.0, 0.5, 0.25, 0.125]
const N_sites_list::Vector{Int64}       = [   4,   4,   4,   8,   4,    4,     4]

# Roggero parameters from 10.1103/PhysRevD.104.103016
const Roggero_table::Array{Float64,2} = [
    −0.5 0 0 1.8 1.0;
    0.0 0 2.10 0 0.36;
    0.01 5.4 0 −8 0.35;
    0.05 2.58 0 0 0.31;
    0.125 1.65 0 1.4 0.32;
    0.25 1.22 0 1.6 0.41;
    0.5 1.06 0 1.4 0.61;
    1.0 0.96 0 0 0.99]
const delta_omega_Rog::Vector{Float64} = Roggero_table[:,1]
const a_t::Vector{Float64} = Roggero_table[:,2]
const b_t::Vector{Float64} = Roggero_table[:,3]
const c_t::Vector{Float64} = Roggero_table[:,4]
const pmin::Vector{Float64} = Roggero_table[:,5]

function main(N_sites, delta_omega)
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δm²= delta_omega # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.

    params = CCNO.Parameters(
        N_sites = N_sites, # number of sites 
        cutoff = 1E-14, # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
        τ = 0.05*CCNO.hbar, # time step 
        ttotal = 5*CCNO.hbar, # total time of evolution 
        tolerance  = 5E-1, # acceptable level of error or deviation from the exact value or solution
        Δx = 1E-3, # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm
        maxdim = 1000, #bond dimension
        m1 = Δm²<0 ? sqrt(-Δm²) : 0.0,
        m2 = Δm²>0 ? sqrt( Δm²) : 0.0,
        L=L,
        Δp = L, # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.  
        periodic = false,  # true = imposes periodic boundary conditions while false doesn't
        checkpoint_every = 100,
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
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
    # conserve_qns=true conserves the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  
    
    # Initialize an array of ones for all N_sites sites
    mu = ones(params.N_sites) # erg
    
    # Create an array of dimension N_sites and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos. 
    N = mu .* fill(((params.Δx)^3 )/(√2 * CCNO.G_F * params.N_sites), params.N_sites)
    
    x = fill(0, params.N_sites) # variable.
    y = fill(0, params.N_sites) # variable.
    z = fill(0, params.N_sites) # variable.

    ψ = productMPS(s, N -> N <= params.N_sites/2 ? "Dn" : "Up")

    # Fixed Constants for Rogerro's fit (only self-interaction term)
    #a_t = 0
    #b_t = 2.105
    #c_t = 0
    

    # array p with N rows and 3 columns, all initialized to 0.0 with colums representing components and rows representing sites
    px_a = fill(0.5, div(params.N_sites,2))
    px_b = fill(0.5, div(params.N_sites,2)) # Δm²/(2.0*omega_b)
    px = vcat(px_a, px_b)
    p = hcat(px,fill(0, params.N_sites), fill(0, params.N_sites))

    energy_sign = [i <= params.N_sites ÷ 2 ? 1 : -1 for i in 1:params.N_sites] # all of the sites are neutrinos

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
    t_array = t_Sz_tot[:, 1]  
    Sz_array = t_Sz_tot[:,2]

    # This function scans through the array, compares each element with its neighbors, 
    # and returns the index of the first local minimum it encounters. 
    # If no local minimum is found, it returns -1 to indicate that.
    function find_first_local_minima_index(arr)
        n = length(arr)
        for i in 2:(n-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
    
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

    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    Rog_index = findfirst(x -> x == delta_omega, delta_omega_Rog)
    t_p_Rog = a_t[Rog_index]*log(params.N_sites) + b_t[Rog_index] * sqrt(params.N_sites) + c_t[Rog_index]
    println("delta_omega=",delta_omega)
    println("N_sites=",N_sites)
    println("t_p_Rog= ",t_p_Rog)
    println()

    # Plotting P_surv vs t
    plot(t_array/CCNO.hbar, -Sz_array, xlabel = "t", ylabel = "-Sz(t)",
         title = "Running main_self_interaction_Rog script", legend = true, size=(800, 600), aspect_ratio=:auto,margin= 10mm, 
         label= ["My_plot_for_N_sites$(params.N_sites)"]) 
    scatter!([t_p_Rog],[-Sz_array[i_first_local_min]], label= ["t_p_Rog"])
    scatter!([t_min/CCNO.hbar],[-Sz_array[i_first_local_min]], label= ["My_t_min)"], legendfontsize=5, legend=:topright)
    # Save the plot as a PDF file # for jenkins archive 
    savefig("main_self_interaction_Rog.pdf")

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    @assert abs(t_min/CCNO.hbar - t_p_Rog) < params.τ/CCNO.hbar + params.tolerance

    # clean up
    rm(params.datadir, recursive=true)
    rm(params.chkptdir, recursive=true)
end 

for list_index in 1:length(N_sites_list)
    main(N_sites_list[list_index], delta_omega_list[list_index])
end


