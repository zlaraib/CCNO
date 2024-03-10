using ITensors
using Plots
using Measures
include("../src/evolution.jl")
include("../src/constants.jl")

# This file evolves the system under the vaccum oscillations + self-interaction
# Hamiltonian and then plots survival probability for 
# multiple system sizes through loops. 

function main(N_sites)
    cutoff = 1E-14
    τ = 0.05
    ttotal = 5
    tolerance  = 5E-1
    Δx = 1E-3
    maxdim = 1000 #bond dimension
    s = siteinds("S=1/2", N_sites; conserve_qns=false)
    a_t = 0.965
    b_t = 0
    c_t = 0
    mu = ones(N_sites)
    N = mu .* fill((Δx)^3/(sqrt(2) * G_F), N_sites)
    B = [0, 0, -1]
    ω_a = fill(2, div(N_sites, 2))
    ω_b = fill(0, div(N_sites, 2))
    Δω = (ω_a - ω_b)/2
    ω = vcat(ω_a, ω_b)
    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites]
    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "..","misc","datafiles","Rog_N_loop", "par_"*string(N_sites), "tt_"*string(ttotal))
    #extract output from the expect.jl file where the survival probability values were computed at each timestep
    Sz_array, prob_surv_array = evolve(s, τ, N, ω, B, N_sites, Δx, ψ, energy_sign, cutoff, maxdim, datadir,ttotal)

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

    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","Rog_N_loop", "par_"*string(N_sites), "tt_"*string(ttotal))

    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)
    plot!(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probability p(t)", title = "Running Rog_particle_loop script", aspect_ratio=:auto, margin= 10mm, legend= true, label= ["My_plot_for_N$(N_sites)"])
    scatter!([t_p_Rog],[prob_surv_array[i_first_local_min]], label= ["t_p_Rog_for_N$(N_sites)"])
    scatter!([t_min],[prob_surv_array[i_first_local_min]], label= ["My_t_min__for_N$(N_sites)"], legendfontsize=5, legend=:topright)
    
end

# Loop from 4 to 8 particles with an increment of 4 particles each time
for N_sites in 4:4:8
    @time main(N_sites)
end

savefig(joinpath(plotdir,"Survival probability vs t (Rog)_loop.pdf"))
