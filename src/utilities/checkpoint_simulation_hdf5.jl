using HDF5
using ITensors
using ITensorMPS

# This file generates functions that read and write into an hdf5 File
# for the values needed for recovery in the evolution process.
function checkpoint_simulation_hdf5(params::CCNO.Parameters, checkpoint_filename::String, state::CCNO.SimulationState, t::Float64, iteration::Int64)
    # Open an HDF5 file for writing (or create it if it doesn't exist)
    f = h5open(checkpoint_filename, "w")
    
    # Write scalar values
    #write(f, "τ(sec)", τ)
    write(f, "L(cm)", params.L)
    write(f, "N_sites(unitless const)", params.N_sites)
    write(f, "Δx(cm)", params.Δx)
    write(f, "m1(erg)", params.m1)
    write(f, "m2(erg)", params.m2)
    write(f, "Δp(cm)", params.Δp)
    write(f, "theta_nu(rad)", params.theta_nu)
    write(f, "cutoff(unitless const)", params.cutoff)
    write(f, "maxdim(unitless const)", params.maxdim)
    write(f, "t(sec)", t)
    write(f, "iteration(unitless const)", iteration)
    
    # Write arrays and vectors
    write(f, "s(unitless array)", state.s)
    write(f, "N(unitless array)", state.N)
    write(f, "p(erg)", state.p)
    write(f, "xyz(cm)", state.xyz)
    write(f, "energy_sign(unitless array)", state.energy_sign)
    
    # Write string values
    write(f, "shape_name(unitless string)", params.shape_name)

    # Write MPS ψ (ITensor type)
    write(f, "ψ(unitless ITensor)", state.ψ)

    # Close the HDF5 file after writing
    close(f)
    
    println("Checkpoint created at $checkpoint_filename, iteration $iteration, time $t")
end

function recover_checkpoint_hdf5(checkpoint_filename::String)
    # Open the HDF5 file for reading
    f = h5open(checkpoint_filename, "r")
    
    # Read scalar values with units
    τ = read(f, "τ(sec)")
    L = read(f, "L(cm)")
    Δx = read(f, "Δx(cm)")
    m1 = read(f, "m1(erg)")
    m2 = read(f, "m2(erg)")
    t_initial = read(f, "t(sec)")
    iteration = read(f, "iteration(unitless const)")
    
    # Read arrays and vectors with units
    s = read(f, "s(unitless array)")
    N = read(f, "N(unitless array)")
    theta_nu = read(f, "theta_nu(rad)")
    p = read(f, "p(erg)")
    xyz = read(f, "xyz(cm)")
    energy_sign = read(f, "energy_sign(unitless array)")
    
    # Read string values
    shape_name = read(f, "shape_name(unitless string)")
    
    # Read MPS ψ (ITensor type)
    ψ = read(f, "ψ(unitless ITensor)", MPS)
    
    # Close the file after reading
    close(f)
    
    println("Recovered checkpoint from $checkpoint_filename, iteration $iteration, time $t_initial")

    s = siteinds(ψ)
    state = simulation_state(ψ=ψ, s=s, energy_sign=energy_sign, N=N, p=p, xyz=xyz)
    
    return state, t_initial, iteration
end
