push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5


# This file evolves the system under the vaccum oscillations + self-interaction
# Hamiltonian and then plots survival probability for 
# multiple system sizes through loops. 
const N_start::Int64 = 4 
const N_step::Int64 = 4
const N_stop::Int64 = 8


function main(N_sites)
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm
    Δm²= 1.0 # erg^2 # Artifically Fixed for Rog bipolar test #change accordingly in gates_fnction too if need be.

    delta_omega = 1.0
    a_t = 0.965
    b_t = 0
    c_t = 0

    params = CCNO.Parameters(
        N_sites = N_sites, # number of sites 
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
    # We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all params.N_sites.
    # conserve_qns=false doesnt conserve the total spin quantum number "S" in the system as it evolves
    s = siteinds("S=1/2", params.N_sites; conserve_qns=false)  

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
    px_a = fill(Δm²/(2.0*2.0*delta_omega), div(params.N_sites,2))
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

    # extract the prob_surv on the first site 
    function find_first_local_minima_index(arr)
        N = length(arr)
        for i in 2:(N-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
 
    t_Sz_tot = readdlm(joinpath(params.datadir, "t_<Sz>.dat"))
    t_array = t_Sz_tot[:, 1]  
    Sz_array = t_Sz_tot[:,2]
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
    t_p_Rog = a_t*log(params.N_sites) + b_t * sqrt(params.N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

   
    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
@assert abs(t_min/CCNO.hbar - t_p_Rog) < params.τ/CCNO.hbar + params.tolerance
end

# Loop from 4 to 8 particles with an increment of 4 particles each time
for N_sites in N_start:N_step:N_stop
    println("")
    println("COMPUTING FOR N_sites=",N_sites)
    main(N_sites)
end
