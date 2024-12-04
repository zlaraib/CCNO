
using HDF5


# This file generates functions that read and write into an hdf5 File
# for the values needed for recovery in the evolution process.
function checkpoint_simulation_hdf5(checkpoint_filename, s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, gates, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, ttotal, t, iteration)
    # Open an HDF5 file for writing (or create it if it doesn't exist)
    f = h5open(checkpoint_filename, "w")
    
    # Write scalar values
    write(f, "τ", τ)
    write(f, "L", L)
    write(f, "N_sites", N_sites)
    write(f, "Δx", Δx)
    write(f, "Δm²", Δm²)
    write(f, "Δp", Δp)
    write(f, "theta_nu", theta_nu)
    write(f, "cutoff", cutoff)
    write(f, "maxdim", maxdim)
    write(f, "t1", t1)
    write(f, "t2", t2)
    write(f, "t", t)
    write(f, "iteration", iteration)
    
    # Write arrays and vectors
    write(f, "s", s)
    write(f, "N", N)
    write(f, "B", B)
    write(f, "p", p)
    write(f, "x", x)
    write(f, "energy_sign", energy_sign)
    
    # Write string values
    write(f, "shape_name", shape_name)

    # Write MPS ψ (ITensor type)
    write(f, "ψ", ψ)

    # Close the HDF5 file after writing
    close(f)
    
    println("Checkpoint created at $checkpoint_filename, iteration $iteration, time $t")
end

function recover_checkpoint_hdf5(checkpoint_filename)
    # Open the HDF5 file for reading
    f = h5open(checkpoint_filename, "r")
    
    # Read scalar values
    τ = read(f, "τ")
    L = read(f, "L")
    N_sites = read(f, "N_sites")
    Δx = read(f, "Δx")
    Δm² = read(f, "Δm²")
    Δp = read(f, "Δp")
    theta_nu = read(f, "theta_nu")
    cutoff = read(f, "cutoff")
    maxdim = read(f, "maxdim")
    t1 = read(f, "t1")
    t2 = read(f, "t2")

    t_initial = read(f, "t")
    iteration = read(f, "iteration")
    
    # Read arrays and vectors
    s = read(f, "s")
    N = read(f, "N")
    B = read(f, "B")
    p = read(f, "p")
    x = read(f, "x")
    energy_sign = read(f, "energy_sign")
    
    # Read string values
    shape_name = read(f, "shape_name")
    
    # Read MPS ψ (ITensor type)
    ψ = read(f, "ψ", MPS)
    
    # Close the file after reading
    close(f)
    
    println("Recovered checkpoint from $checkpoint_filename, iteration $iteration, time $t_initial")

    return s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, t1, t2, t_initial, iteration
end
