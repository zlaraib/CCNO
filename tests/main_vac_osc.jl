using ITensors
using Plots 
using Measures 
using LinearAlgebra
using DelimitedFiles
using DelimitedFiles
include("../src/evolution.jl")
include("../src/constants.jl")
include("../src/shape_func.jl")
include("../src/momentum.jl")

# We are simulating the time evolution of a 1D spin chain with N_sites sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.


function main()
  N_sites = 6 # number of sites, #variable
  cutoff = 1E-14 # specifies a truncation threshold for the SVD in MPS representation #variable
  τ = 0.1 # time step # sec #variable
  ttotal = 5.0 # total time of evolution #sec #variable
  Δx = 1E-3 # length of the box of interacting neutrinos at a site in cm #variable
  tolerance  = 1E-5 # acceptable level of error or deviation from the exact value or solution #variable
  Δm² = 2*π # Fixed for this vacuum oscillation case for omega =pi. dont change it to keep consistent results. 
  maxdim = 1000 # max bond dimension in MPS truncation
  L = 1 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
  Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent.
  t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
  t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test
  
  # Make an array of 'site' indices and label as s 
  # conserve_qns=false doesnt conserve the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N_sites; conserve_qns=false)  #Fixed

  # Initialize an array of zeros for all N_sites particles
  mu = zeros(N_sites) #Fixed
                                
  # Create an array of dimension N and fill it with the value 1/(sqrt(2) * G_F). This is the number of neutrinos.
  N = mu.* fill((Δx)^3/(sqrt(2) * G_F), N_sites) 
  
  # Create a B vector which would be same for all N particles 
  theta_nu = π/4 # mixing_angle #rad 
  B = [sin(2*theta_nu), 0, -cos(2*theta_nu)] # is equivalent to B = [1, 0, 0] # variable. But only other case that can be tested from this file is B = [0,0,-1] for which theta_nu = π/4.
  B = B / norm(B)

  x = fill(rand(), N_sites) # variable.
  y = fill(rand(), N_sites) # variable.
  z = zeros(N_sites) # variable.
 
  # Generate an Nx3 array for p with random values
  p = ones(N_sites, 3) # variable, but will need to make sure that p_vector.jl file if statment stays constsnet 
  energy_sign = fill(1, N_sites) # all of the sites are neutrinos

  #Select a shape function based on the shape_name variable form the list defined in dictionary in shape_func file
  shape_name = "none"  # variable.

  # Initialize psi to be a product state (First half to be spin down and other half to be spin up)
  ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up") # Fixed to produce consistent results for the test assert conditions 

  # Specify the relative directory path
  datadir = joinpath(@__DIR__, "..","misc","datafiles","vac_osc", "par_"*string(N_sites), "tt_"*string(ttotal))

  #extract output for the survival probability values at each timestep
  Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, Im_Ω = evolve(s, τ, N, B,L, N_sites, 
  Δx,Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal)
                      
  expected_sz_array = Float64[]
  expected_sz= Float64[]
  
  for t in 0.0:τ:ttotal

    i = 1 # change it according to the corresponding site number in the expect function 
    if B[1] == 1

      # Compute the expected value based on the derived analytic formula
      expected_sz = -0.5 * cos(ω[i] * t)

    end
    if B[3] == -1

      # Compute the expected value based on the derived analytic formula
      expected_sz = -0.5

    end

    push!(expected_sz_array, expected_sz)

  end

  # Check if every element in Sz_array is less than tolerance away from the corresponding element in expected_sz_array
  # for B vector in x, it checks that the value of Sz at the first spin site oscillates between -0.5 and 0.5 
  # for B vector in -z, it checks that the value of Sz at the firstspin site never oscillates from -0.5 
  @assert all(abs.(Sz_array .- expected_sz_array) .< tolerance)
  
  # Specify the relative directory path
  plotdir = joinpath(@__DIR__, "..","misc","plots","vac_osc", "par_"*string(N_sites), "tt_"*string(ttotal))
  
  # check if a directory exists, and if it doesn't, create it using mkpath
  isdir(plotdir) || mkpath(plotdir)
  plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script",
  legend = true, size=(700, 600), aspect_ratio=:auto,left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, 
  bottom_margin = 10mm, label = "My_sz") 
  plot!(0.0:τ:τ*(length(Sz_array)-1), expected_sz_array, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script", 
  legendfontsize=8, legend=:topright, label = "Expected_sz from Sakurai") 
  # Save the plot as a PDF file
  savefig(joinpath(plotdir,"<Sz> vs t(vac_osc).pdf"))

  plot(0.0:τ:τ*(length(Sy_array)-1), Sy_array, xlabel = "t", ylabel = "<Sy>", legend = false, size=(800, 600),left_margin = 20mm,
   right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, aspect_ratio=:auto,margin= 10mm) 
  #Save the plot as a PDF file
  savefig(joinpath(plotdir,"<Sy> vs t(vac_osc).pdf"))

  plot(0.0:τ:τ*(length(Sx_array)-1), Sx_array, xlabel = "t", ylabel = "<Sx>", legend = false, size=(800, 600),left_margin = 20mm,
   right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, aspect_ratio=:auto,margin= 10mm) 
  #Save the plot as a PDF file
  savefig(joinpath(plotdir,"<Sx> vs t(vac_osc).pdf"))
end

@time main()