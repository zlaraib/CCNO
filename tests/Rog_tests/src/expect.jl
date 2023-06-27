include("gates_function.jl")  # Include the gates_functions.jl file

function calc_expect(s, tau, N, cutoff, ttotal)
    # Initialize psi to be a product state (alternating up and down)
    psi = productMPS(s, n -> isodd(n) ? "Dn" : "Up")

    # center site because the Sz operator can be arbitralily chosen to be placed at this site at each time step
    # to match up with Rog results we changed the placement of Sz operator on first site instead of middle.
    c = div(N, 2)

    # Create empty array to store sz values 
    Sz_array = Float64[]
    # Create empty array to store survival probability values 
    prob_surv_array = Float64[]

    # extract the gates array generated in the gates_function file
    gates = create_gates(s, N, tau)

    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
    for t in 0.0:tau:ttotal
        # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site on the chain
        sz = expect(psi, "Sz"; sites=1)
        # add an element sz to the end of Sz array 
        push!(Sz_array, sz)

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
        psi = apply(gates, psi; cutoff)
        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
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