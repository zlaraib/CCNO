using DelimitedFiles
include("gates_function.jl")  # Include the gates_functions.jl file
include("H_MPO.jl") 
include("momentum.jl")
include("constants.jl")
"""
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the evolve function.                                                                   
    s = site index array (dimensionless and unitless)          
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    B = array of normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
    N_sites = Total no.of sites (dimensionless and unitless)
    Δx = length of the box of interacting neutrinos at a site (cm)
    Δm² = difference in mass squared (erg^2)
    p = array of momentum vectors (erg)
    x = array of positions of sites (cm)
    Δp = width of shape function (cm)
    shape_name = name of the shape function (string) ["none","triangular","flat_top"]
    τ = time step (sec)
    energy_sign = array of sign of the energy (1 or -1): 1 for neutrinos and -1 for anti-neutrinos
    maxdim = max bond dimension in MPS truncation (unitless and dimensionless)
    cutoff = truncation threshold for the SVD in MPS representation (unitless and dimensionless)
    periodic = boolean indicating whether boundary conditions should be periodic
"""

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 
function evolve(s, τ, N, B,L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, periodic= true)

    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(datadir) || mkpath(datadir)

    # Create empty array s to...
    Sz_array = Float64[] # to store sz values 
    Sy_array = Float64[] # to store sy values 
    Sx_array = Float64[] # to store sx values 
    t_array = [] # to store t values 
    prob_surv_array = Float64[]   # to store survival probability values 
    x_values = []  # to store x values for all sites 
    pₓ_values = [] # to store px vector values for all sites
    ρₑₑ_array = Float64[] # to store ρₑₑ values
    ρ_μμ_array = Float64[] # to store ρ_μμ values
    ρₑμ_array = Float64[] # to store ρₑμ values
    
    ρₑμ_at_t1 = nothing  # Initialize a variable to store ρₑμ at t1
    ρₑμ_at_t2 = nothing  # Initialize a variable to store ρₑμ at t2
    Δt = t2 - t1 #time difference between growth rates 

    # extract the gates array generated in the gates_function file
    gates = create_gates(s, N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name,L, τ, energy_sign, periodic)
    # H = Hamiltonian_mpo(s, N, B, N_sites, Δx, Δm², p, x, Δp, shape_name,L, τ, energy_sign, periodic)
    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(p,N_sites) 
    p̂ₓ= [sub_array[1] for sub_array in p̂]
    
     # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
     for t in 0.0:τ:ttotal
        push!(x_values, copy(x))  # Record x values at each time step
        px = p[:, 1]  # Extracting the first column (which corresponds to px values)
        push!(pₓ_values, copy(px)) # Record px values at each time step

        for i in 1:N_sites
            x[i] += p̂ₓ[i] * c * τ
            if periodic
                # wrap around position from 0 to domain size L
                x[i] = mod(x[i],L)

                # Checking if the updated x[i] satisfies the boundary conditions
                @assert (x[i] >= 0 && x[i] <= L)
            end
        end
        
        if shape_name!=="none" #specfic to the inhomogenous case Test4
            # compute the avg expectation value of Sz at all sites
            sz_tot = expect(ψ, "Sz")
            sz = mean(sz_tot)
        else 
            # compute expectation value of Sz (inbuilt operator in ITensors library) at the first site on the chain
            sz = expect(ψ, "Sz"; sites=1)
        end 

        # compute expectation value of sy and sx using S+ and S- (inbuilt operator in ITensors library) at the first site on the chain
        if p == zeros(N_sites, 3) #for rogerro's case only (b/c S+ S- needed to keep conservation of QN number)
            sy = -0.5 *im * (expect(complex(ψ), "S+"; sites=1) - expect(complex(ψ), "S-"; sites=1)) #re-check
            sx = 0.5 * (expect(ψ, "S+"; sites=1) + expect(ψ, "S-"; sites=1)) #recheck
        else 
            sy = expect(complex(ψ), "Sy"; sites=1)
            sx = expect(ψ, "Sx"; sites=1)
        end

        # add an element sz to ... 
        push!(Sz_array, sz) # .. the end of Sz array 
        push!(Sy_array, sy) # .. the end of Sy array 
        push!(Sx_array, sx) # .. the end of Sx array 
        
        # survival probability for a (we took first) neutrino to be found in its initial flavor state (in this case a spin down)
        prob_surv = 0.5 * (1 - 2 * sz)

        # add an element prob_surv to the end of  prob_surv_array 
        push!(prob_surv_array, prob_surv)

        if B[1] == 1
            println("$t $sz")
        else 
            println("$t $prob_surv")
        end
        
        # recall that in our code sigma_z = 2*Sz so make sure these expressions are consistent with "Sz in ITensors" 
        ρₑₑ = ( (2 * sz) + 1)/2 
        push!(ρₑₑ_array,abs(ρₑₑ))
        ρ_μμ = ( (-2 * sz) + 1)/2 
        push!(ρ_μμ_array,abs(ρ_μμ))
        ρₑμ = sqrt(sx^2 + sy^2)
        push!(ρₑμ_array, ρₑμ)

        # Check if the current time is approximately t1
        if abs(t - t1) < τ / 2
            ρₑμ_at_t1 = ρₑμ
        end
        
        # Check if the current time is approximately t2
        if abs(t - t2) < τ / 2
            ρₑμ_at_t2 = ρₑμ
        end
        # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
        # and if so, it executes the break statement, which causes the loop to terminate early.
        t ≈ ttotal && break

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        ψ = apply(gates, ψ; cutoff, maxdim)
        # ψ = tdvp(H, ψ,  -im *τ;   nsweeps=1,
        # reverse_step=true,outputlevel=1)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(ψ)
    end
    t_array = 0.0:τ:ttotal

    # After the time evolution loop, calculate and print the growth rate from the ratio if both values have been captured
    if ρₑμ_at_t1 !== nothing && ρₑμ_at_t2 !== nothing
        global Im_Ω = (1/Δt ) * log(ρₑμ_at_t2 / ρₑμ_at_t1)
        println("Growth rate of flavor coherence of ρₑμ at t2 to ρₑμ at t1: $Im_Ω")
    else
        println("ρₑμ was not captured at both t1 and t2.")
    end

    # Writing data to files with corresponding headers
    fname1 = joinpath(datadir, "t_<Sz>_<Sy>_<Sx>.dat")
    writedlm(fname1, [t_array Sz_array Sy_array Sx_array])
    fname2 = joinpath(datadir, "t_probsurv.dat")
    writedlm(fname2, [t_array prob_surv_array])
    fname3 = joinpath(datadir, "t_xsiteval.dat")
    writedlm(fname3, [t_array x_values])
    fname4 = joinpath(datadir, "t_pxsiteval.dat")
    writedlm(fname4, [t_array pₓ_values])
    fname5 = joinpath(datadir, "t_ρₑₑ.dat")
    writedlm(fname5, [t_array ρₑₑ_array])
    fname6 = joinpath(datadir, "t_ρ_μμ.dat")
    writedlm(fname6, [t_array ρ_μμ_array])
    fname7 = joinpath(datadir, "t_ρₑμ.dat")
    writedlm(fname7, [t_array ρₑμ_array])
    return Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, Im_Ω
end


