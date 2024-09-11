using ITensors
using Plots
using Measures
using LinearAlgebra
using DelimitedFiles
using Statistics
using Random

include("main_Hamiltonian.jl")
include("../src/constants.jl")

""" Richers(2021) Test 4 initial conditions: """
N_sites_eachflavor= 5 # total sites/particles that evenly spaced "for each (electron) flavor" 
N_sites = 2* (N_sites_eachflavor) # total particles/sites for all neutrino and anti neutrino electron flavored
τ = 5E-13 # time step to include 50 steps every 10 picoseconds # sec # variable
ttotal = 9.0E-11 # total time of evolution # sec #variable
tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution #variable
m1 = -0.008596511*eV #eV  1st mass eigenstate of neutrino
m2 = 0*eV #eV  2nd mass eigenstate of neutrino
Δm² = (m2^2-m1^2) # mass square difference # (erg^2)
maxdim = 1 # max bond dimension in MPS truncation
cutoff = 1e-100 # specifies a truncation threshold for the SVD in MPS representation (SMALL CUTOFF = MORE ENTANGLEMENT) #variable
L = 1 # cm # domain size # (aka big box length)
n_νₑ =  4.891290848285061e+32 # cm^-3 # number density of electron flavor neutrino
n_νₑ̄ =  n_νₑ # cm^-3 # number density of electron flavor antineutrino
Eνₑ =  50.0*MeV # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
Eνₑ̄ = -1 * Eνₑ # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign
Δx = L/N_sites # length of the box of interacting neutrinos at a site in cm  #variable
# Δx = L
theta_nu = 1.74532925E-8  #1e-6 degrees # mixing_angle # = 1.74532925E-8 radians 
B = [-sin(2 * theta_nu), 0, cos(2 * theta_nu)]  # actual b vector that activates the vacuum oscillation term in Hamiltonian
B = B / norm(B) 
#Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
shape_name = "triangular"  # Change this to the desired shape name #variable 
Δp = Δx # width of shape function  # cm #variable
# Δp = L/N_sites
t1 = 33e-15 #choose initial time for growth rate calculation
t2 = 53e-15 #choose final time for growth rate calculation
periodic = true  # true = imposes periodic boundary conditions while false doesn't
# analytic_growth_rate=  (abs(m2^2 - m1^2)/ (2*hbar* Eνₑ)) # analytic growth rate #fix it for inhomo from paper

# generate x_array such that the first particle is at position L/(2*N_sites) while subsequent particles are at a position incremental by L/N_sites. # grid style
function generate_x_array(N_sites, L)
    return [(i - 0.5) * L / N_sites for i in 1:N_sites]
end

x = generate_x_array(N_sites, L)
y = generate_x_array(N_sites, L)
z = generate_x_array(N_sites, L)

#generate a momentum array in px direction that depicts the energy of neutrinos and anti-neutrinos in opposing beams
function generate_px_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(Eνₑ, half_N_sites); fill(Eνₑ̄, half_N_sites)]
end

function generate_py_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end

function generate_pz_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end

# p matrix with numbers generated from the p_array for all components (x, y, z) #sherood has 
p = hcat(generate_px_array(N_sites), generate_py_array(N_sites), generate_pz_array(N_sites))

# Create an array with the first half as 1 and the rest as -1
energy_sign = [i <= N_sites ÷ 2 ? -1 : 1 for i in 1:N_sites] # half sites are (e) neutrinos with positive 1 entry while other half is anti (e) neutrinos with negative 1 entry

# s is an array of spin 1/2 tensor indices (Index objects) which will be the site or physical indices of the MPS.
# We overload siteinds function, which generates custom Index array with Index objects having the tag of total spin quantum number for all N.
# conserve_qns=true conserves the total spin quantum number "Sz" in the system as it evolves,
# i.e. examples of conservation of quantum numbers are the total number of neutrino particles, or the total of all S_z components of this system of spins
# conserving total Sz requires Sx and Sy in terms of S+ and S- by design choice.

s = siteinds("S=1/2", N_sites; conserve_qns=false) #fixed #switched conserve_qns to false to avoid fluxes error in expect function

# Initialize psi to be a product state (Of all electron flavor neutrino i.e. spin up in Richers notation which is equivalently half spin up and half chain spin down in my TN notation)
ψ = productMPS(s, n -> n <= N_sites/2 ? "Up" : "Dn")

function generate_B_pert(α)
    # Generate two random perturbations for x and y
    # x_pert = α *   (2 * rand() - 1)  # Random number between -α and α
    # y_pert = α *   (2 * rand() - 1)  # Random number between -α and α
    x_pert = rand() 
    y_pert = rand() 
    # Calculate the z component to maintain normalization
    z_pert = sqrt(max(0, 1 - x_pert^2 - y_pert^2))

    # Return the B_pert vector
    return [x_pert, y_pert, z_pert]
end
# 
α = 1e-6 # perturbation strength as mentioned in the paper
# Generate the perturbed B vector scaled by α  as mentioned in the paper
B_pert = generate_B_pert(α)

# Since the perturbation is small, B_pert should already be normalized, but you can normalize again for precision
B_pert = B_pert / norm(B_pert)
# println("B_pert= ", B_pert)  

# Perturb the state via one-body Hamiltonian
ψ₀= evolve_perturbation(s, τ, B_pert, α, x, L, N_sites, ψ, cutoff, maxdim, energy_sign, ttotal)

@time main(s, τ, B,L, N_sites, N_sites_eachflavor, tolerance,
                n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, ttotal,periodic)

# @assert abs((Im_Ω - analytic_growth_rate)/ analytic_growth_rate) < tolerance #fix it for inhomo from paper