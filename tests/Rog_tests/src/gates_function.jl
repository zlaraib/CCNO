using ITensors

function create_gates(s, N, tau)
    gates = ITensor[]

    for i in 1:(N-1)
        for j in i+1:N
            s1 = s[i]
            s2 = s[j]
            hj = 2.0/N * (op("Sz", s1) * op("Sz", s2) +
                         1/2 * op("S+", s1) * op("S-", s2) +
                         1/2 * op("S-", s1) * op("S+", s2))
            Gj = exp(-im * tau/2 * hj)
            push!(gates, Gj)
        end
    end

    append!(gates, reverse(gates))
    return gates
end
