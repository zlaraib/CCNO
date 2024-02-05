using ITensors
using Plots 
using Measures 
using LinearAlgebra
include("../src/evolution.jl")
include("../src/constants.jl")

# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

function main()
  N_sites = 4 # number of sites
  cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation
  τ = 0.1 # time step 
  ttotal = 5.0 # total time of evolution 
  Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm 
  tolerance  = 1E-5 # acceptable level of error or deviation from the exact value or solution

  # Make an array of 'site' indices and label as s 
  # conserve_qns=false doesnt conserve the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N_sites; conserve_qns=false)  

  # Initialize an array of zeros for all N_sites particles
  mu = zeros(N_sites)
                                
  # Create an array of dimension N_sites and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos.
  n = mu.* fill((Δx)^3/(sqrt(2) * G_F), N_sites)
      
  # Create a B vector which would be same for all N_sites particles 
  B = [1, 0, 0]          
  
  # Create arrays ω_a and ω_b
  ω_a = fill(π, div(N_sites, 2))
  ω_b = fill(π, div(N_sites, 2))

  # Concatenate ω_a and ω_b to form ω with N_sites elements. Each element of the array is a const pi.
  ω = vcat(ω_a, ω_b)

  gates = create_gates(s, n, ω, B, N_sites, Δx, τ)
  
  # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
  ψ = productMPS(s, n -> n <= N_sites/2 ? "Dn" : "Up")

  #extract output from the expect.jl file where the survival probability values were computed at each timestep
  Sz_array, prob_surv_array = evolve(s, τ, n, ω, B, N_sites, Δx, ψ, cutoff, ttotal)
  
  expected_sz_array = Float64[]
  expected_sz= Float64[]
  
  for t in 0.0:τ:ttotal

    for i in 1:(N_sites-1)

        if ω[i] != 0

          if B[1] == 1

            # Compute the expected value based on the derived analytic formula
            expected_sz = -0.5 * cos(ω[i] * t)

          end
          if B[3] == -1

            # Compute the expected value based on the derived analytic formula
            expected_sz = -0.5

          end

        end

    end

    push!(expected_sz_array, expected_sz)

  end

  # Check if every element in Sz_array is less than tolerance away from the corresponding element in expected_sz_array
  # for B vector in x, it checks that the value of Sz at the first spin site oscillates between -0.5 and 0.5 
  # for B vector in -z, it checks that the value of Sz at the firstspin site never oscillates from -0.5 
  @assert all(abs.(Sz_array .- expected_sz_array) .< tolerance)

  # Plotting P_surv vs t
  plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script",legend = true, size=(700, 600), aspect_ratio=:auto,margin= 10mm, label = "My_sz") 
  plot!(0.0:τ:τ*(length(Sz_array)-1), expected_sz_array, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script", legendfontsize=8, legend=:topright, label = "Expected_sz from Sakurai") 
  # Save the plot as a PDF file
  savefig("<Sz> vs t (only vacuum oscillation term plot).pdf")
end

@time main()


