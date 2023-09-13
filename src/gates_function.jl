using ITensors
include("constants.jl")

# Expected units of the quantities used in gates                                                                    
# del_m2 = ergs squared
# hbar =  erg s
# c = cm/s
# G_F = erg cm^3           
# ∇x = cm
# E = erg

function create_gates(s, n, ω, B, N, Δx, τ)
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[]                                                              
    
    for i in 1:(N-1)
        for j in i+1:N
            #s1, s2 are non-neighbouring spin site/indices from the s array
            s_i = s[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
            s_j = s[j]
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
            
            hj = 
            ((2.0* √2 * G_F * (n[i]+ n[j])/(2*((Δx)^3)) * 1/N) * 
            (op("Sz", s_i) * op("Sz", s_j) +
             1/2 * op("S+", s_i) * op("S-", s_j) +
             1/2 * op("S-", s_i) * op("S+", s_j)))
             
             if ω[i] != 0 && ω[j] != 0
                hj += (1/(N-1))* 
                ((ω[i] * B[1] * op("Sx", s_i)* op("Id", s_j))  + (ω[j] * op("Sx", s_j) * op("Id", s_i))) + 
                ((ω[i] * B[2] * op("Sy", s_i)* op("Id", s_j))  + (ω[j] * op("Sy", s_j) * op("Id", s_i))) +
                ((ω[i] * B[3] * op("Sz", s_i)* op("Id", s_j))  + (ω[j] * op("Sz", s_j) * op("Id", s_i))) 
                     
             end
            
            
            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * τ/2 * hj)

            # The push! function adds (appends) an element to the end of an array;
            # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
            # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
            push!(gates, Gj)
        end
    end

    # append! adds all the elements of a gates in reverse order (i.e. (N,N-1),(N-1,N-2),...) to the end of gates array.
    # appending reverse gates to create a second-order Trotter-Suzuki integration
    append!(gates, reverse(gates))
    return gates
end


