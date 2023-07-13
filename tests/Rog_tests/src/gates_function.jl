using ITensors

function create_gates(s, n, del_m2, B, E, mu, N, del_x, G_F, tau)
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[]
    
    for i in 1:(N-1)
        for j in i+1:N
            #s1, s2 are non-neighbouring spin site/indices from the s array
            s1 = s[i]
            s2 = s[j]
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B[i]) == 1
            
            #B[i] = normalize(B[i])
            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.
            # mu pairs divided by 2 to avoid double counting
            hj = -((del_m2/(2*E[i])) * B[i][3] * op("Sz", s1)* op("Id", s2)) + 
                    ((2.0* √2 * G_F * (n[i]+ n[j])/(2*((del_x)^3)) * 1/N) * 
                    (op("Sz", s1) * op("Sz", s2) +
                     1/2 * op("S+", s1) * op("S-", s2) +
                     1/2 * op("S-", s1) * op("S+", s2)))

            # writing in terms of mu 
            # hj = (2.0* (mu[i]+ mu[j])/(2) * 1/N) * ((op("Sz", s1) * op("Sz", s2) +
            # 1/2 * op("S+", s1) * op("S-", s2) +
            # 1/2 * op("S-", s1) * op("S+", s2)))
            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * tau/2 * hj)

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
