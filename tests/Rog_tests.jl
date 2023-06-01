using ITensors
using Plots

N = 10 # number of sites
cutoff = 1E-8 # specifies a truncation threshold for the SVD in MPS representation
tau = 0.05 # time step 
ttotal = 20 # total time of evolution 

s = siteinds("S=1/2", N; conserve_qns=true)  

gates = ITensor[] 

for i in 1:(N-1) 
    # s1 = s[j] 
    # s2 = s[j + 1]
    
    for j in i+1:N
        s1 = s[i]
        s2 = s[j]
        # println(i)
        # println(j)
    # Rog sigma = 2 * sz of our code here 
        hj = 2.0/N *
        (op("Sz", s1) * op("Sz", s2) +
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2))
        #println(2/N)
        Gj = exp(-im * tau / 2 * hj)  
        
        push!(gates, Gj) 
        #println(gates)
    end

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
    global psi = apply(gates, psi; cutoff)
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