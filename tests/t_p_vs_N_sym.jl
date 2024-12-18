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
# Hamiltonian and then plots the system size N_sites on x axis while minimum time tp  
# on the y-axis. It then compares the values of our calculations with Rogerros 
# calculations in Table I.
# Here a fixed, but symmetric δω is used i.e. ω_a = - ω_b for a given δω = 0.25.

#changing variables here 
Δω = 0.25
N_start = 4 
N_step= 4
N_stop= 8
ttotal = 10

function main(N_sites, Δω)
    cutoff = 1E-14
    τ = 0.05
    tolerance  = 5E-1
    Δx = 1E-3
    Δm² = Δω
    maxdim = 1000 #bond dimension
    L = 10 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test  
    periodic = false
    s = siteinds("S=1/2", N_sites; conserve_qns=false)
    # check for Δω = 0.25
    a_t = 1.224
    b_t = 0
    c_t = 1.62
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
    shape_name = "none" 
    datadir = joinpath(@__DIR__, "datafiles","Rog_Fig_3b", "par_"*string(N_sites), "tt_"*string(ttotal))
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
    return t_p_Rog, t_min
end
datadir = joinpath(@__DIR__, "datafiles","Rog_Fig_3b", "N_start_"*string(N_start), "N_stop_"*string(N_stop))

# Arrays to store t_p_Rog and t_min for each N_sites
t_p_Rog_array = Float64[]
t_min_array = Float64[]

# Loop from N_start to N_stop particles with an increment of N_step particles each time
for N_sites in N_start: N_step:N_stop
    t_p_Rog, t_min = @time main(N_sites, Δω)
    push!(t_p_Rog_array, t_p_Rog)
    push!(t_min_array, t_min)
end
N_values = range(N_start, stop=N_stop, step=N_step)
if save_data 
    save_data = isdir(datadir) || mkpath(datadir)
    fname1 = joinpath(datadir, "Nvals_tpRog_tpmine.dat")
    writedlm(fname1, [N_values t_p_Rog_array t_min_array])
end 

if save_plots_flag
    plotdir = joinpath(@__DIR__, "plots","Rog_Fig_3b", "N_start_"*string(N_start), "N_stop_"*string(N_stop))
        
    # check if a directory exists, and if it doesn't, create it using mkpath
    save_plot_flag =  isdir(plotdir) || mkpath(plotdir)
    # Create the plot
    plot(N_values, t_p_Rog_array, label="Rog_tp", xlabel="N_sites", ylabel="Minimum Time(t_p)", title = "Table I Rogerro(2021) \n expanded for a symmetric δω=0.25", legend=:topleft, aspect_ratio=:auto,margin= 10mm)
    plot!(N_values, t_min_array, label="Our_tp")
    savefig(joinpath(plotdir,"t_p_vs_N_for_symmetric_del_w.pdf"))
end

# Create the plot
plot(N_values, t_p_Rog_array, label="Rog_tp", xlabel="N_sites", ylabel="Minimum Time(t_p)", title = "Table I Rogerro(2021) \n expanded for a symmetric δω=0.25", legend=:topleft, aspect_ratio=:auto,margin= 10mm)
plot!(N_values, t_min_array, label="Our_tp")
savefig("t_p_vs_N_for_symmetric_del_w.pdf")