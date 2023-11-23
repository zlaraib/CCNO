include("gates_function.jl")  # Include the gates_functions.jl file
include("momentum.jl")
"""
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the gates function.                                                                   
    s = site index array (dimensionless and unitless)          
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    B = array of normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
    N_sites = Total no.of sites (dimensionless and unitless)
    Δx = length of the box of interacting neutrinos at a site (cm)
    del_m2 = difference in mass squared (erg^2)
    p = array of momentum vectors (erg)
    x = array of positions of sites (cm)
    Δp = width of shape function (cm)
    shape_name = name of the shape function (string) ["none","triangular","flat_top"]
    τ = time step (sec)
    energy_sign = array of sign of the energy (1 or -1): 1 for neutrinos and -1 for anti-neutrinos
    maxdim = max bond dimension in MPS truncation (unitless and dimensionless)
    cutoff = truncation threshold for the SVD in MPS representation (unitless and dimensionless)
"""

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, N, B, N_sites, Δx, del_m2, p, x, Δp, ψ, shape_name, energy_sign, cutoff, maxdim, datafile, ttotal)
    
    # Create empty array to store sz values 
    Sz_array = Float64[]
    Sy_array = Float64[]
    Sx_array = Float64[]
    # Create empty array to store survival probability values 
    prob_surv_array = Float64[]

    x_values = []  # Store x values for plotting

    # extract the gates array generated in the gates_function file
    gates = create_gates(s, N, B, N_sites, Δx,del_m2, p, x, Δp, shape_name, τ, energy_sign)

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p_hat = momentum(p,N_sites) 
    p_x_hat = [sub_array[1] for sub_array in p_hat]

     # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
     for t in 0.0:τ:ttotal
        push!(x_values, copy(x))  # Record x values at each time step
        x .+=  (p_x_hat .* τ)  # displacing particle's position at each timestep 


        # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site on the chain

        sz = expect(ψ, "Sz"; sites=1)

        sy = expect(complex(ψ), "Sy"; sites=1)
        #sy =  inner(complex(ψ), apply(op("Sy",s,1), complex(ψ)))
        sx = expect(ψ, "Sx"; sites=1)

        # add an element sz to the end of Sz array 
        push!(Sz_array, sz)
        push!(Sy_array, sy)
        push!(Sx_array, sx)
        
        # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
        prob_surv = 0.5 * (1 - 2 * sz)

        # add an element prob_surv to the end of  prob_surv_array 
        push!(prob_surv_array, prob_surv)

        if B[1] != 0
            # Write the values to the data file
            println(datafile, "$t $sz $sy $sx")
            flush(datafile) # Flush the buffer to write immediately
            # println("$t $sz")
        else 
            println(datafile, "$t $prob_surv")
            flush(datafile) # Flush the buffer to write immediately
        end

        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(gates, ψ; cutoff, maxdim)
        

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end
    t_array = 0.0:τ:ttotal

    return Sz_array, Sy_array, Sx_array, prob_surv_array, x_values
end


