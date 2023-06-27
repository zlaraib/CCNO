using ITensors

function create_gates(s, N, tau)
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[]

    for i in 1:(N-1)
        for j in i+1:N
            #s1, s2 are non-neighbouring spin site/indices from the s array
            s1 = s[i]
            s2 = s[j]
            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.
            hj = 2.0/N * (op("Sz", s1) * op("Sz", s2) +
                         1/2 * op("S+", s1) * op("S-", s2) +
                         1/2 * op("S-", s1) * op("S+", s2))

            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * tau/2 * hj)

            # The push! function adds (appends) an element to the end of an array;
            # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
            # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
            push!(gates, Gj)
        end
    end

    # append! adds all the elements of a gates in reverse order (i.e. (N,N-1),(N-1,N-2),...) to the end of gates array.
    append!(gates, reverse(gates))
    return gates
end
