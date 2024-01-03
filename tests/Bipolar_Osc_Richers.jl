using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics

#Needs revision

include("main_self_interaction.jl")
include("../src/constants")
# dk what bipolar osc time is in the test
""" Richers(2021) Test 2 initial conditions: """
N_sites_eachflavor= 1 # total sites/particles that evenly spaced "for each (electron) flavor" 
N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
τ = 2E-13 # time step to include 50 steps every 10 picoseconds # sec # variable
ttotal = 9E-11 # total time of evolution # sec #variable
tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
m1 = 8.6e3 #eV  1st mass eigenstate of neutrino
m2 = 0 #eV  2nd mass eigenstate of neutrino
#Δm² = (m2^2-m1^2) # mass square difference # (erg^2) #not used in the test
#add other parameters like charcteristic bipolar time and H neu/ H vac =10  ...?
maxdim = 1 # max bond dimension in MPS truncation
cutoff = 0 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
L = 1 # cm # domain size # (aka big box length)
n_νₑ =  (10* (m2-m1)^2 *c^4)/(2*√2*G_F*Eνₑ) # cm^-3 # number density of electron flavor neutrino
n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
Eνₑ =  50.0e6 # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
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

# s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
# We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
# conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
# i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
# conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

# Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up)
ψ = productMPS(s, N_sites -> "Up") 

@time main(N_sites_eachflavor,τ,ttotal,tolerance,Δm²,maxdim,cutoff,x, p, ψ, L,n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,B_pert,B,shape_name,periodic)
