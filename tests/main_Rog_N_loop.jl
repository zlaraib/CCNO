using ITensors
using Plots
using Measures
using DelimitedFiles

"""
For github unit tests runs: 
src_dir = ../
save_data and save_plots_flag should be false to run test files. 
"""

src_dir = "../"
save_data = false  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs


"""
For science runs: 
src_dir = /home/zohalaraib/Oscillatrino/ # should be changed to users PATH
save_data and save_plots_flag should be true to run test files. 

"""
# src_dir= "/home/zohalaraib/Oscillatrino/" # should be changed to users PATH
# save_data = true  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
# save_plots_flag = true # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    
include(src_dir * "src/evolution.jl")
include(src_dir * "src/constants.jl")
include(src_dir * "Utilities/save_plots.jl")
include(src_dir * "Initializations/initial_cond.jl")

# This file evolves the system under the vaccum oscillations + self-interaction
# Hamiltonian and then plots survival probability for 
# multiple system sizes through loops. 

N_start = 4 
N_step= 4
N_stop= 8
ttotal = 5

function main(N_sites)
    cutoff = 1E-14
    τ = 0.05
    tolerance  = 5E-1
    Δx = 1E-3
    Δm² = 2
    maxdim = 1000 #bond dimension
    L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.  
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test
    periodic = false  # true = imposes periodic boundary conditions while false doesn't
    s = siteinds("S=1/2", N_sites; conserve_qns=false)
    a_t = 0.965
    b_t = 0
    c_t = 0
    mu = ones(N_sites)
    N = mu .* fill(((Δx)^3 )/(√2 * G_F * N_sites), N_sites)
    theta_nu = 0 # mixing_angle #rad 
    B = [sin(2*theta_nu), 0, -cos(2*theta_nu)] # is equivalent to B = [0, 0, -1] # fixed for Rogerro's case
    B = B / norm(B)

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_p_array(N_sites),fill(0, N_sites), fill(0, N_sites))
    x = fill(rand(), N_sites) # variable.
    y = fill(rand(), N_sites) # variable.
    z = fill(rand(), N_sites) # variable.
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites]
    shape_name = "none"  # Change this to the desired shape name # variable.

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "datafiles","Rog_N_loop", "par_"*string(N_sites), "tt_"*string(ttotal))
    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, Im_Ω = evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, save_data , periodic)
    
    function find_first_local_minima_index(arr)
        N = length(arr)
        for i in 2:(N-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
    
    i_first_local_min = find_first_local_minima_index(prob_surv_array)
    if i_first_local_min != -1
        println("Index of the first local minimum: ", i_first_local_min)
    else
        println("No local minimum found in the array.")
    end

    t_min = τ * i_first_local_min - τ
    println("Corresponding time of first minimum index= ", t_min)
    t_p_Rog = a_t*log(N_sites) + b_t * sqrt(N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)
    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    @assert abs(t_min - t_p_Rog) <  τ + tolerance 


    plot!(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probability p(t)", title = "Running Rog_particle_loop script", aspect_ratio=:auto, margin= 10mm, legend= true, label= ["My_plot_for_N$(N_sites)"])
    scatter!([t_p_Rog],[prob_surv_array[i_first_local_min]], label= ["t_p_Rog_for_N$(N_sites)"])
    scatter!([t_min],[prob_surv_array[i_first_local_min]], label= ["My_t_min__for_N$(N_sites)"], legendfontsize=5, legend=:topright)
    savefig("Survival probability vs t (Rog)_loop.pdf")
end

# Loop from 4 to 8 particles with an increment of 4 particles each time
for N_sites in N_start:N_step:N_stop
    if save_plots_flag
        # Specify the relative directory path
        global plotdir = joinpath(@__DIR__, "plots","Rog_N_loop", "par_"*string(N_sites), "tt_"*string(ttotal))
        save_plot_flag = isdir(plotdir) || mkpath(plotdir)
        savefig(joinpath(plotdir,"Survival probability vs t (Rog)_loop.pdf"))
    end 
    
    @time main(N_sites)
end

