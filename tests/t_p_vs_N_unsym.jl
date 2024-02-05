using ITensors
using Plots
using Measures
using DelimitedFiles
include("../src/evolution.jl")
include("../src/constants.jl")

# This file evolves the system under the vaccum oscillations + self-interaction
# Hamiltonian and then plots the system size N_sites on x axis while minimum time tp  
# on the y-axis. It then compares the values of our calculations with Rogerros 
# calculations in Table I.
# Here a fixed, but unsymmetric δω is used i.e. ω_a = - ω_b for a given δω = 0.25.

#changing variables here
N_start = 4 
N_step= 4
N_stop= 24

function main(N_sites)
    cutoff = 1E-14
    τ = 0.05
    ttotal = 10
    tolerance  = 5E-1
    Δx = 1E-3
    s = siteinds("S=1/2", N_sites; conserve_qns=false)
    # check for Δω = 0.25
    a_t = 1.224
    b_t = 0
    c_t = 1.62
    mu = ones(N_sites)
    N = mu .* fill((Δx)^3/(sqrt(2) * G_F), N_sites)
    B = [0, 0, -1]
    ω_a = fill(0.5, div(N_sites, 2))
    ω_b = fill(0, div(N_sites, 2))
    Δω = (ω_a - ω_b)/2
    ω = vcat(ω_a, ω_b)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites]
    Sz_array, prob_surv_array = evolve(s, τ, N, ω, B, N_sites, Δx, ψ, energy_sign, cutoff, ttotal)
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

datadir = joinpath(@__DIR__, "..","misc","datafiles","Rog", "N_start_"*string(N_start), "N_stop_"*string(N_stop))
isdir(datadir) || mkpath(datadir)

# Arrays to store t_p_Rog and t_min for each N_sites
t_p_Rog_array = Float64[]
t_min_array = Float64[]

# Loop from N_start to N_stop particles with an increment of N_step particles each time
for N_sites in N_start: N_step:N_stop
    t_p_Rog, t_min = @time main(N_sites)
    push!(t_p_Rog_array, t_p_Rog)
    push!(t_min_array, t_min)
end
N_values = range(N_start, stop=N_stop, step=N_step)
fname1 = joinpath(datadir, "N_tpRog_tpmine.dat")
writedlm(fname1, [N_values t_p_Rog_array t_min_array])

plotdir = joinpath(@__DIR__, "..","misc","plots","Rog", "N_start_"*string(N_start), "N_stop_"*string(N_stop))
    
# check if a directory exists, and if it doesn't, create it using mkpath
isdir(plotdir) || mkpath(plotdir)

# Create the plot
plot(N_values, t_p_Rog_array, label="Rog_tp", xlabel="N_sites", ylabel="Minimum Time (t_p)", title = "Table I Rogerro(2021) \N expanded for a unsymmetric δω=0.25", legend=:topleft, aspect_ratio=:auto,margin= 10mm)
plot!(N_values, t_min_array, label="Our_tp")
savefig(joinpath(plotdir,"t_p_vs_N_for_unsymmetric_del_w.pdf"))