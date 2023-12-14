using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics

include("main_self_interaction.jl")
include("../src/constants.jl")

""" Richers(2021) Test 4 initial conditions: """
N_sites_eachflavor= 50 # total sites/particles that evenly spaced "for each (electron) flavor" 
N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
τ = 2E-13 # time step to include 50 steps every 10 picoseconds # sec # variable
ttotal = 9E-11 # total time of evolution # sec #variable
tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
m1 = 8.6e3 #eV  1st mass eigenstate of neutrino
m2 = 0 #eV  2nd mass eigenstate of neutrino
Δm² = 0 # mass square difference # fixed for 'only' self interactions # (erg^2)
maxdim = 1 # max bond dimension in MPS truncation
cutoff = 0 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
L = 1 # cm # domain size # (aka big box length)
n_νₑ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
n_νₑ̄ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor antineutrino
Eνₑ =  50.0e6 # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
B_pert = [0.02, -0.02, 1] # Create a B vector that allows for perturbation to inital state in different directions #variable 
theta_nu = 10e-6 # mixing_angle
B = [sin(2*theta_nu), 0, -cos(2*theta_nu)] # actual b vector that activates the vacuum oscillation term in Hamiltonian
#Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
shape_name = "triangular"  # Change this to the desired shape name #variable 
Δp = 1/N_sites_eachflavor # width of shape function  # cm #variable
periodic = true  # true = imposes periodic boundary conditions while false doesn't
# generate x_array such that the first particle is at position L/(2*N_sites) while subsequent particles are at a position incremental by L/N_sites. # grid style
function generate_x_array(N_sites, L)
    return [(i - 0.5) * L / N_sites for i in 1:N_sites]
end

x = generate_x_array(N_sites, L)
y = fill(rand(), N_sites) #variable
z = fill(rand(), N_sites) #variable

#generate a momentum array that depicts the energy of neutrinos and anti-neutrinos in opposing beams
function generate_p_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(Eνₑ, half_N_sites); fill(Eνₑ̄, half_N_sites)]
end

# p matrix with numbers generated from the p_array for all components (x, y, z)
p = hcat(generate_p_array(N_sites), generate_p_array(N_sites), generate_p_array(N_sites))
# Create an array with the first half as 1 and the rest as -1
energy_sign = [i <= N_sites ÷ 2 ? 1 : -1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry
s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function
ψ = productMPS(s, N_sites -> "Up") # inital state all in spin up direction i.e. all electron flavor

@time main(N_sites_eachflavor,τ,ttotal,tolerance,Δm²,maxdim,cutoff,x, p, ψ, L,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,B_pert,B,shape_name,periodic)
