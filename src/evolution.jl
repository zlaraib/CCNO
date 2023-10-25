include("gates_function.jl")  # Include the gates_functions.jl file
include("momentum.jl")
"""
    Expected units of the quantities defined in the files in tests directory that are being used in the evolve function                                                                   
    s = site index array (dimensionless and unitless) 
    τ = time step (sec)      
    n = no.of neutrinos (dimensionless and unitless)
    ω = vacuum oscillation angular frequency (rad/s)
    B = Normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
    N = Total no.of sites (dimensionless and unitless)
    Δx = length of the box of interacting neutrinos at a site (cm) 
    cutoff = truncation threshold for the SVD in MPS (unitless, number)
    ttotal = ttotal time (sec)
"""

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, n, B, N, Δx, del_m2, p, x, Δp, ψ, shape_name, cutoff, tolerance, ttotal)
    
    # Create empty array to store sz values 
    Sz_array = Float64[]
    # Create empty array to store survival probability values 
    prob_surv_array = Float64[]

    # extract the gates array generated in the gates_function file
    gates = create_gates(s, n, B, N, Δx,del_m2, p, x, Δp, shape_name, τ)

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p_hat = momentum(p,N) 
    p_x_hat = [sub_array[1] for sub_array in p_hat]

     # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
     for t in 0.0:τ:ttotal
        
        x .+=  (p_x_hat .* t)  # displacing particle's position at each timestep #Discuss with sherwood if this should be t ot tau or maybe initialize the x array just before this line. b/c x2-x1 = pxt means that theres uneven time diff between successive displacements of particles

        # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site on the chain
        sz = expect(ψ, "Sz"; sites=1)
        # add an element sz to the end of Sz array 
        push!(Sz_array, sz)
        
        # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
        prob_surv = 0.5 * (1 - 2 * sz)

        # add an element prob_surv to the end of  prob_surv_array 
        push!(prob_surv_array, prob_surv)

        if n == fill(0, N)
            println("$t $sz")
        else println("$t $prob_surv")
        end

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It automatically handles truncating the MPS and handles the non-nearest-neighbor gates in this example.
        ψ = apply(gates, ψ; cutoff)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end
    t_array = 0.0:τ:ttotal

    # Writing data to files with corresponding headers
    fname2 = joinpath(datadir, "t_probsurv.dat")
    writedlm(fname2, [t_array prob_surv_array])
    return Sz_array, prob_surv_array
end


