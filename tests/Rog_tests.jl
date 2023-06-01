using ITensors
using Plots

N = 96 # number of sites
cutoff = 1E-4 # specifies a truncation threshold for the SVD in MPS representation
tau = 0.05 # time step 
ttotal = 10 # total time of evolution 

s = siteinds("S=1/2", N; conserve_qns=true)  

gates = ITensor[] 

for j in 1:(N - 1)
    s1 = s[j] 
    s2 = s[j + 1]
    
    # hj = pi * (
    #     op("Sz", s1) * op("Id", s2)  + op("Sz", s2) * op("Id", s1)
    # ) + 
    # op("Sz", s1) * op("Sz", s2) +
    # 1 / 2 * op("S+", s1) * op("S-", s2) +    # 1 / 2 * op("S-", s1) * op("S+", s2)
    
    hj =
    op("Sz", s1) * op("Sz", s2) +
    1 / 2 * op("S+", s1) * op("S-", s2) +
    1 / 2 * op("S-", s1) * op("S+", s2)

    Gj = exp(-im * tau / 2 * hj)  

    push!(gates, Gj) 
end

append!(gates, reverse(gates)) 

psi = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

c = div(N, 2) # c = N/2 

Sz_array = Float64[] 
P_elec_array = Float64[]  # Array to store P_elec values
prob_surv_array = Float64[] 

for t in 0.0:tau:ttotal
    sz = expect(psi, "Sz"; sites=1)
    push!(Sz_array, sz) 
    #println("$t $sz")
    t ≈ ttotal && break  
    psi = apply(gates, psi; cutoff)
    normalize!(psi) 

    # P_z = 2 * sz
    # P_elec = 0.5 * (1 + P_z)
    # push!(P_elec_array, P_elec)  # Store P_elec values

    prob_surv= 0.5 * (1- 2*sz)
    push!(prob_surv_array, prob_surv)
    println("$t $prob_surv")
    
end

@assert prob_surv_array[1] == 1.0
# Plotting P_elec vs t
#plot(0.0:tau:tau*(length(P_elec_array)-1), P_elec_array, xlabel = "t", ylabel = "P_elec", legend = false)
plot(0.0:tau:tau*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "prob_surv", legend = false)