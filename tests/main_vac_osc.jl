using ITensors
using Plots 
using Measures 
using LinearAlgebra
using DelimitedFiles
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")
include("../src/momentum.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.


function main()
  N = 4 # number of sites, #variable
  cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation #variable
  τ = 0.1 # time step #variable
  ttotal = 5.0 # total time of evolution #variable
  Δx = 1E-3 # length of the box of interacting neutrinos at a site/shape function width of neutrinos in cm #variable
  tolerance  = 1E-5 # acceptable level of error or deviation from the exact value or solution #variable
  Δp = 1 # shape function width #variable
  del_m2 = 2*π # Fixed for this vacuum oscillation case for omega =pi. dont change it to keep consistent results. 


  # Make an array of 'site' indices and label as s 
  # conserve_qns=false doesnt conserve the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N; conserve_qns=false)  #Fixed

  # Initialize an array of zeros for all N particles
  mu = zeros(N) #Fixed
                                
  # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos.
  n = mu.* fill((Δx)^3/(sqrt(2) * G_F), N) 
      
  # Create a B vector which would be same for all N particles 
  B = [1, 0, 0] # variable. But only other case that can be tested from this file is B = [0,0,-1].

  x = zeros(N) # variable.
  y = fill(rand(), N) # variable.
  z = zeros(N) # variable.
 
  # Generate an Nx3 array for p with random values
  p = ones(N, 3) # variable, but will need to make sure that p_vector.jl file if statment stays constsnet 

  # extract output of p_hat and p_mod for the p vector defined above for all sites. 
  p_mod, p_hat = momentum(p,N)

  #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
  shape_name = "none"  # variable.

  # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
  ψ = productMPS(s, n -> n <= N/2 ? "Dn" : "Up") # Fixed to produce consistent results for the test assert conditions 

  #extract output from the expect.jl file where the survival probability values were computed at each timestep
  Sz_array, prob_surv_array = evolve(s, τ, n, B, N, Δx,del_m2, p, p_mod, p_hat, x, Δp, ψ, shape_name, cutoff, tolerance, ttotal)

  expected_sz_array = Float64[]
  expected_sz= Float64[]
  
  for t in 0.0:τ:ttotal

    for i in 1:(N-1)

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
  savefig(joinpath(plotdir,"<Sz> vs t (only vacuum oscillation term plot).pdf"))
end

@time main()