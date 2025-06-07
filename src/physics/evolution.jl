using DelimitedFiles
using ITensors
using ITensorMPS

# This file generates the evolve function which evolves the ψ state in time and computes the expectation values of Sz at each time step, along 
# with their survival probabilities. The time evolution utilizes the unitary operators created as gates from the create_gates function.
# The <Sz> and Survival probabilities output from this function are unitless. 

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
function evolve(params::CCNO.Parameters, state::CCNO.SimulationState)

    t_initial::Float64 = 0.0
    iteration::Int64 = 0

    if params.do_recover
        if isfile(params.recover_file)
            println("Recovering from checkpoint: $params.recover_file")
            state, t_initial, iteration = recover_checkpoint_hdf5(params)
        else
            error("Checkpoint file not found")
        end
    else
        # make sure the data files don't already exist if starting from scratch
        if ispath(params.chkptdir)
            error("Path already exists: "*params.chkptdir)
        end
        if ispath(params.datadir)
            error("Path already exists: "*params.datadir)
        end
        mkpath(params.chkptdir)
        mkpath(params.datadir)
    end    

    # get sites in the correct order
    sort_sites!(state)

    # Compute and print survival probability (found from <Sz>) at each time step then apply the gates to go to the next time
    for t in t_initial:params.τ:params.ttotal
        println("iteration= $iteration time= $t")

        store_data(params.datadir, t, state)

        # extract output of p_hat and p_mod for the p vector defined above for all sites. 
        p_mod, p̂ = momentum(state.p)

        # move particles
        state.xyz += p̂ * c * params.τ
        if params.periodic
            state.xyz = mod.(state.xyz, params.L)
            @assert all(state.xyz .>= 0 .&& state.xyz .<= params.L)
        end

        sort_sites!(state)

        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.
        # apply each gate in gates(ITensors array) successively to the wavefunction ψ (MPS)(it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
        # The apply function is a matrix-vector multiplication operation that is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
        # It truncates the MPS according to the set cutoff and maxdim for all the non-nearest-neighbor gates.


        # extract the gates array generated in the gates_function file
        gates_1site, gates_2site_even, gates_2site_odd, gates_2site_other = create_gates(params, state)
        

        state.ψ = apply(        gates_2site_other , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(        gates_2site_even  , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(        gates_2site_odd   , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(        gates_1site       , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(reverse(gates_2site_odd)  , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(reverse(gates_2site_even) , state.ψ; params.cutoff, params.maxdim)
        state.ψ = apply(reverse(gates_2site_other), state.ψ; params.cutoff, params.maxdim)

        # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
        # This is necessary to ensure that the MPS represents a valid quantum state.
        normalize!(state.ψ)

        iteration += 1

        if iteration % params.checkpoint_every == 0 
            checkpoint_filename = joinpath(params.chkptdir, "checkpoint.chkpt.it" * lpad(iteration, 6, "0") * ".h5")
            checkpoint_simulation_hdf5(params, checkpoint_filename, state, t, iteration)
        end
        
    end

end


