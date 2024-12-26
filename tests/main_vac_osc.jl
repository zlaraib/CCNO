using ITensors
using Plots 
using Measures 
using LinearAlgebra
using DelimitedFiles
using HDF5

"""
For github unit tests runs: 
src_dir = ../
save_data and save_plots_flag should be false to run test files. 
"""

src_dir = "../"
save_data = false  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs


"""
For science runs: 
src_dir = /home/zohalaraib/Oscillatrino/ # should be changed to users PATH
save_data and save_plots_flag should be true to run test files. 

"""
# src_dir= "/home/zohalaraib/Oscillatrino/" # should be changed to users PATH
# save_data = true  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
# save_plots_flag = true # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs
    

include(src_dir * "src/evolution.jl")
include(src_dir * "src/constants.jl")
include(src_dir * "src/momentum.jl")
include(src_dir * "Utilities/save_plots.jl")
include(src_dir * "src/chkpt_hdf5.jl")
include(src_dir * "Utilities/save_datafiles.jl")

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
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
  periodic = false  # true = imposes periodic boundary conditions while false doesn't
      
  checkpoint_every = 4
  do_recover = false
  recover_type = "auto" 
  recover_iteration = 80 # change it to the iteration you want to recover from, for manual iteration. Currently auto recovery already recovers from last iteration (i.e. recover_iteration = -1 for auto recovery). 
  
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
  datadir = joinpath(@__DIR__, "datafiles", "par_"*string(N_sites))
  chkptdir = joinpath(@__DIR__, "checkpoints", "par_"*string(N_sites))

  #extract output for the survival probability values at each timestep
  Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover = evolve(s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal, chkptdir, checkpoint_every, do_recover, recover_type, recover_iteration, save_data::Bool, periodic)
    
  # extract the Sz on the first site 
  Sz_array_site1= [row[1] for row in Sz_array]
  expected_sz_array_site1 = Float64[]
  expected_sz_site1= Float64[]
  
  if !do_recover 
    for t in 0.0:τ:ttotal

      i = 1 # change it according to the corresponding site number in the expect function 
      if B[1] == 1

        # Compute the expected value based on the derived analytic formula
        expected_sz_site1 = -0.5 * cos(ω[i] * t)

      end
      if B[3] == -1

        # Compute the expected value based on the derived analytic formula
        expected_sz_site1 = -0.5

      end

      push!(expected_sz_array_site1, expected_sz_site1)

    end
    # Check if every element in Sz_array is less than tolerance away from the corresponding element in expected_sz_array
    # for B vector in x, it checks that the value of Sz at the first spin site oscillates between -0.5 and 0.5 
    # for B vector in -z, it checks that the value of Sz at the firstspin site never oscillates from -0.5 
    @assert all(abs.(Sz_array_site1 .- expected_sz_array_site1) .< tolerance)
    
  elseif do_recover 
    for t in t_recover:τ:ttotal

      i = 1 # change it according to the corresponding site number in the expect function 
      if B[1] == 1

        # Compute the expected value based on the derived analytic formula
        expected_sz_site1 = -0.5 * cos(ω[i] * t)

      end
      if B[3] == -1

        # Compute the expected value based on the derived analytic formula
        expected_sz_site1 = -0.5

      end

      push!(expected_sz_array_site1, expected_sz_site1)

    end
    # Check if every element in Sz_array is less than tolerance away from the corresponding element in expected_sz_array
    # for B vector in x, it checks that the value of Sz at the first spin site oscillates between -0.5 and 0.5 
    # for B vector in -z, it checks that the value of Sz at the firstspin site never oscillates from -0.5 
    @assert all(abs.(Sz_array_site1 .- expected_sz_array_site1) .< tolerance)
    
  end

  if save_plots_flag
    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "plots", "par_"*string(N_sites))
    # Read the data files
    t_Sz_tot = readdlm(joinpath(datadir, "t_<Sz>.dat"))
    t_Sy_tot = readdlm(joinpath(datadir, "t_<Sy>.dat"))
    t_Sx_tot = readdlm(joinpath(datadir, "t_<Sx>.dat"))
    t_probsurv_tot = readdlm(joinpath(datadir, "t_probsurv.dat"))
    t_xsiteval = readdlm(joinpath(datadir, "t_xsiteval.dat"))
    t_pxsiteval = readdlm(joinpath(datadir, "t_pxsiteval.dat"))
    t_ρₑₑ_tot = readdlm(joinpath(datadir, "t_ρₑₑ.dat"))
    t_ρ_μμ_tot = readdlm(joinpath(datadir, "t_ρ_μμ.dat"))
    t_ρₑμ_tot = readdlm(joinpath(datadir, "t_ρₑμ.dat"))

    # Extract time array and corresponding values for plotting
    t_array = t_Sz_tot[:, 1]  

    #Extract the array for first site only
    Sz_array = t_Sz_tot[:,2]
    Sy_array = t_Sy_tot[:,2]
    Sx_array= t_Sx_tot[:,2] 
    prob_surv_array = t_probsurv_tot[:, 2]
    ρₑₑ_array = t_ρₑₑ_tot[:, 2]
    ρ_μμ_array = t_ρ_μμ_tot[:, 2]
    ρₑμ_array = t_ρₑμ_tot[:, 2]
    
    x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
    pₓ_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first

    # # Parsing arrays containing strings like "[1.0,", into a numeric array suitable for plotting
    Sz_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sz_array]
    Sy_array = [parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sy_array]
    Sx_array =[parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in Sx_array]
    prob_surv_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in prob_surv_array]
    ρₑₑ_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρₑₑ_array]
    ρ_μμ_array = [ parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρ_μμ_array]
    ρₑμ_array =[parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in ρₑμ_array]

    save_plots(τ, N_sites,L,t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array,datadir, plotdir, save_plots_flag)
  end 
  if !save_plots_flag
    plot(t_array, Sz_array_site1, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script",
    legend = true, size=(700, 600), aspect_ratio=:auto,left_margin = 20mm, right_margin = 5mm, top_margin = 5mm, 
    bottom_margin = 10mm, label = "My_sz_site1") 
    plot!(t_array, expected_sz_array_site1, xlabel = "t", ylabel = "<Sz>", title = "Running main_vac_osc script", 
    legendfontsize=8, legend=:topright, label = "Expected_sz from Sakurai_site1") 
    # Save the plot as a PDF file # for jenkins archive 
    savefig("<Sz>_site1 vs t.pdf")
  end
end

@time main()