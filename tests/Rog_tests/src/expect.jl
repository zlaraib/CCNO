include("gates_function.jl")  # Include the gates_functions.jl file

function calc_expect(s, tau, N, cutoff, ttotal)
    psi = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    c = div(N, 2)
    Sz_array = Float64[]
    prob_surv_array = Float64[]

    gates = create_gates(s, N, tau)

    for t in 0.0:tau:ttotal
        sz = expect(psi, "Sz"; sites=1)
        push!(Sz_array, sz)

        t ≈ ttotal && break
        psi = apply(gates, psi; cutoff)
        normalize!(psi)

        prob_surv = 0.5 * (1 - 2 * sz)
        push!(prob_surv_array, prob_surv)
        println("$t $prob_surv")

        # Write the values to the datafile
        # println(datafile, "$t $prob_surv")
        # flush(datafile)
    end

    return Sz_array, prob_surv_array
end