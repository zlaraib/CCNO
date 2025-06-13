# Struct of parameters. Allow specifying arguments by keywords
Base.@kwdef struct Parameters
    shape_name::String  # Specifies the shape function used in the nu-nu Hamiltonian ["flat_top", "triangular", "none"]
    geometric_name::String # Specifies the name of the geometric function used in the nu-nu Hamiltonian ["physical", "none"]
    recover_file::String # path to the file to recover from
    datadir::String # path to write data to
    chkptdir::String # path to write checkpoints to
    plotdir::String # path to write plots to

    do_recover::Bool # False: start anew at t=0. True: recover from recover_file
    checkpoint_every::Int64 # number of timesteps between writing checkpoints
    save_plots_flag::Bool # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    periodic::Bool  # true = imposes periodic boundary conditions while false doesn't

    N_sites::Int64 # total particles/sites for all neutrino and anti neutrino electron flavored
    maxdim::Int64 # max bond dimension in MPS truncation
    cutoff::Float64 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)

    τ::Float64 # [seconds] size of time step
    ttotal::Float64 # [seconds] total time of evolution
    tolerance::Float64 # acceptable level of error or deviation from the exact value or solution. TODO: remove from parameters struct
    m1::Float64 # [eV]  1st mass eigenstate of neutrino
    m2::Float64 # [eV]  2nd mass eigenstate of neutrino
    theta_nu::Float64 # mixing_angle
    Δx::Float64 # [cm] assumed size of each particle cloud, used for determining neutrino density
    Δp::Float64 # [cm] width of shape function. Numerical means of making interactions more or less local without changing physical phase space volume represented by site
    α::Float64 # perturbation strength as mentioned in the paper for the inhomogenous Richers test
    L::Float64 # [cm] total size of the domain if using periodic boundary conditions
end
