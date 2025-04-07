# Struct of parameters. Allow specifying arguments by keywords
Base.@kwdef struct Parameters
    shape_name::String  # Change this to the desired shape name
    recover_file::String
    datadir::String
    chkptdir::String
    plotdir::String

    do_recover::Bool    
    save_plots_flag::Bool # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    periodic::Bool  # true = imposes periodic boundary conditions while false doesn't

    N_sites::Int # total particles/sites for all neutrino and anti neutrino electron flavored
    maxdim::Int # max bond dimension in MPS truncation
    checkpoint_every::Int

    τ::Float64 # time step to include 50 steps every 10 picoseconds # sec
    ttotal::Float64 # total time of evolution # sec
    tolerance::Float64 # acceptable level of error or deviation from the exact value or solution
    m1::Float64 #eV  1st mass eigenstate of neutrino
    m2::Float64 #eV  2nd mass eigenstate of neutrino
    theta_nu::Float64 # mixing_angle
    cutoff::Float64 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT)
    Δp::Float64 # width of shape function  # cm
    α::Float64 # perturbation strength as mentioned in the paper for the inhomogenous Richers test 

end
