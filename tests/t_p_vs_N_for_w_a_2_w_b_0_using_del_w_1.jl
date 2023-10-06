using ITensors
using Plots
using Measures
include("../src/evolution.jl")
include("../src/constants.jl")

function main(N)
    cutoff = 1E-14
    τ = 0.05
    ttotal = 5
    tolerance  = 5E-1
    Δx = 1E-3
    s = siteinds("S=1/2", N; conserve_qns=false)
    a_t = 0.965
    b_t = 0
    c_t = 0
    mu = ones(N)
    n = mu .* fill((Δx)^3/(sqrt(2) * G_F), N)
    B = [0, 0, -1]
    ω_a = fill(2, div(N, 2))
    ω_b = fill(0, div(N, 2))
    Δω = (ω_a - ω_b)/2
    ω = vcat(ω_a, ω_b)
    ψ = productMPS(s, n -> n <= N/2 ? "Dn" : "Up")
    Sz_array, prob_surv_array = evolve(s, τ, n, ω, B, N, Δx, ψ, cutoff, ttotal)
    function find_first_local_minima_index(arr)
        n = length(arr)
        for i in 2:(n-1)
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
    t_p_Rog = a_t*log(N) + b_t * sqrt(N) + c_t
    println("t_p_Rog= ",t_p_Rog)

    return t_p_Rog, t_min
end

# Arrays to store t_p_Rog and t_min for each N
t_p_Rog_array = Float64[]
t_min_array = Float64[]

# Loop from 4 to 8 particles with an increment of 4 particles each time
for N in 4:4:20
    t_p_Rog, t_min = main(N)
    push!(t_p_Rog_array, t_p_Rog)
    push!(t_min_array, t_min)
end

# Create the plot
plot(4:4:20, t_p_Rog_array, label="t_p_Rog", xlabel="N", ylabel="Minimum Time (t_p)", title = "Running t_p_vs_N_for_w_a_2_w_b_0_using_del_w_1 script", aspect_ratio=:auto,margin= 10mm)
plot!(4:4:20, t_min_array, label="t_min")
savefig("t_p_vs_N_for_w_a_2_w_b_0_using_del_w_1.pdf")