using HDF5


# This file generates functions that read and write into an hdf5 File
# for the values needed for recovery in the evolution process.
function checkpoint_simulation_hdf5(params::CCNO.parameters, checkpoint_filename, s, N, B, L, Δx, Δm², p, x, ψ, energy_sign, t, iteration)
    # Open an HDF5 file for writing (or create it if it doesn't exist)
    f = h5open(checkpoint_filename, "w")
    
    # Write scalar values
    #write(f, "τ(sec)", τ)
    write(f, "L(cm)", L)
    write(f, "N_sites(unitless const)", params.N_sites)
    write(f, "Δx(cm)", Δx)
    write(f, "Δm²(erg^2)", Δm²)
    write(f, "Δp(cm)", params.Δp)
    write(f, "theta_nu(rad)", params.theta_nu)
    write(f, "cutoff(unitless const)", params.cutoff)
    write(f, "maxdim(unitless const)", params.maxdim)
    write(f, "t(sec)", t+params.τ)
    write(f, "iteration(unitless const)", iteration)
    
    # Write arrays and vectors
    write(f, "s(unitless array)", s)
    write(f, "N(unitless array)", N)
    write(f, "B(unitless array)", B)
    write(f, "p(erg)", p)
    write(f, "x(cm)", x)
    write(f, "energy_sign(unitless array)", energy_sign)
    
    # Write string values
    write(f, "shape_name(unitless string)", params.shape_name)

    # Write MPS ψ (ITensor type)
    write(f, "ψ(unitless ITensor)", ψ)

    # Close the HDF5 file after writing
    close(f)
    
    println("Checkpoint created at $checkpoint_filename, iteration $iteration, time $t")
end

function recover_checkpoint_hdf5(checkpoint_filename)
    # Open the HDF5 file for reading
    f = h5open(checkpoint_filename, "r")
    
    # Read scalar values with units
    τ = read(f, "τ(sec)")
    L = read(f, "L(cm)")
    Δx = read(f, "Δx(cm)")
    Δm² = read(f, "Δm²(erg^2)")
    t_initial = read(f, "t(sec)")
    iteration = read(f, "iteration(unitless const)")
    
    # Read arrays and vectors with units
    s = read(f, "s(unitless array)")
    N = read(f, "N(unitless array)")
    B = read(f, "B(unitless array)")
    p = read(f, "p(erg)")
    x = read(f, "x(cm)")
    energy_sign = read(f, "energy_sign(unitless array)")
    
    # Read string values
    shape_name = read(f, "shape_name(unitless string)")
    
    # Read MPS ψ (ITensor type)
    ψ = read(f, "ψ(unitless ITensor)", MPS)
    
    # Close the file after reading
    close(f)
    
    println("Recovered checkpoint from $checkpoint_filename, iteration $iteration, time $t_initial")

    return s, τ, N, B, L, Δx, Δm², p, x, ψ, energy_sign, t_initial, iteration
end
