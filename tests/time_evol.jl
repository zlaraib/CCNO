using ITensors
using Plots 
using Measures 
using LinearAlgebra
include("Rog_tests/src/gates_function.jl")
include("Rog_tests/src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.
let
  N = 2 # number of sites
  cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation
  tau = 0.1 # time step 
  ttotal = 5.0 # total time of evolution 
  del_x = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
  del_m2= 1.93E-23 # mass difference between the second and first neutrino mass eigenstates (associated with solar neutrino oscillations) in ergs^2
  tolerance  = 5E-1 # acceptable level of error or deviation from the exact value or solution

  # Make an array of 'site' indices and label as s 
  # conserve_qns=true conserves the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N; conserve_qns=false)  

  # Initialize an array of ones for all N particles
  mu = zeros(N)
                                
  # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F)
  n = mu.* fill((del_x)^3/(sqrt(2) * G_F), N)
      
  # Create an array B with N elements. Each element of the array is a vector [1, 0, 0]
  B = fill([1, 0, 0], N)              
      
  # Create an array of neutrino vaccum energy
  E = fill(4/(del_m2),N)
  
  gates = create_gates(s, n, del_m2, B, E, N, del_x, tau)
  
  
  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  # center site because the Sz operator is arbitralily chosen to be placed at this site at each time step 
  c = div(N, 2) # c = N/2 
  
  # Create empty array to store sz values 
  Sz_array = Float64[] 

  # Compute and print <Sz> at each time step then apply the gates to go to the next time
  for t in 0.0:tau:ttotal
    # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the first site of the chain
    sz = expect(psi, "Sz"; sites=c)
    # add an element sz to the end of Sz array  
    push!(Sz_array, sz) 
    println("$t $sz")

    # Compute the expected value based on the derived analytic formula
    expected_sz = 0.5 * cos(pi * t)

    # Checking that the value of Sz at the center spin oscillates between 0.5 and -0.5 
    # Compare the actual value with the expected value using a tolerance
    @assert abs(sz - expected_sz) < tolerance

    # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
    # and if so, it executes the break statement, which causes the loop to terminate early.
    t ≈ ttotal && break  

    # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
    # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
    # It automatically handles truncating the MPS (and can even handle non-nearest-neighbor gates but this feature is not used in this example).
    psi = apply(gates, psi; cutoff)

    # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
    # This is necessary to ensure that the MPS represents a valid quantum state.
    normalize!(psi) 
  end

end
