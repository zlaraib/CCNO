
include("main_self_interaction.jl")

    """ My initial conditions""" 
    N_sites_eachflavor= 3 # total sites/particles that evenly spaced "for each (electron) flavor" 
    N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
    τ = 5e-12 # time step # sec # variable
    ttotal = 1e-9 # total time of evolution # sec #variable
    tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
    Δm² = 0 # mass square difference # fixed for 'only' self interactions # (erg^2)
    maxdim = 1 # max bond dimension in MPS truncation
    cutoff = 0 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
    L = 1 # cm # domain size # (aka big box length)
    n_νₑ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
    n_νₑ̄ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor antineutrino
    Eνₑ =  50.0e6 # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
    Eνₑ̄ = -1* Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
    #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
    shape_name = "triangular"  # Change this to the desired shape name #variable 
    Δp = 1/N_sites_eachflavor # width of shape function  # cm #variable
    periodic = true  # true = imposes periodic boundary conditions while false doesn't


    @time main(N_sites_eachflavor,τ,ttotal,tolerance,Δm²,maxdim,cutoff,L,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,shape_name,periodic)